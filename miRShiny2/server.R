library(ggplot2)
library(limma)
library(RColorBrewer)
library(DESeq2)
library(edgeR)
library(shiny)
library(reshape)
library(viridis)
library(sva)
library(grid)
library(reader)
library(NMF)
library(plyr)
library(RnaSeqSampleSize)
library(circlize)
#library(pwr)

#######
shinyServer(function(input, output) {
  ##########################
  
  storageValues <- reactiveValues()
  
  theme_set(theme_bw())
  
  theme_update(axis.title = element_text(face = "italic"), plot.title = element_text(hjust = 0.5, size = 15, face = "plain"), legend.position = "bottom",
               legend.direction = "horizontal", panel.border = element_rect(fill = NA, colour = "grey50"))
  
  #############################
  
  #make Conditions File selecton, add voom to normalize matrix function
  
  rawElist <- eventReactive(input$upload, ignoreNULL = FALSE, {
    if ((!is.null(input$files) && !is.null(input$condition)) || input$sampleData) {
      withProgress(message = 'Reading files', value = 0, {
        #read in matrix file
        if(input$sampleData){
          matrix <- read.csv(file = "MCW_CRC_Exo_Matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
        } else{
          sep <- get.delim(win(input$files$datapath), n = 20, skip = 1)
          matrix <-
            read.table(
              win(input$files$datapath),
              row.names = 1,
              header = TRUE,
              sep = sep,              
              check.names = FALSE
            )
        }
        incProgress(1 / 6, detail = 'Reading matrix')
        
        
        #clean matrix data
        incProgress(1 / 6, detail = 'Cleaning matrix data')
        if(length(unique(colnames(matrix))) != ncol(matrix)){
          showModal(modalDialog(
            title = "Error Processing Input",
            "Sample Labels (Columns) Duplicated or Missing."
          ))
          return(NULL)
        }
        matrix[is.na(matrix)] <- 0
        matrix[is.null(matrix)] <- 0
        matrix[matrix < 0] <- 0
        matrix <- floor(matrix)
        
        #read in Conditions File
        incProgress(1 / 6, detail = 'Getting conditions')
        if(input$sampleData){
          condition <-
            read.csv(
              file = "MCW_CRC_Exo_Conditions.csv", 
              header = TRUE,
              check.names = FALSE
            )
        } else{
          condSep <- get.delim(win(input$condition$datapath), n = 10, skip = 1)
          if(is.na(condSep)){
            condSep <- NULL
          }
          condition <-
            read.table(
              win(toString(input$condition$datapath)), 
              header = TRUE,
              sep = condSep,
              check.names = FALSE
            )
        }
        
        #checking objects
        incProgress(1 / 6, detail = 'Checking objects')
        if (ncol(matrix) != nrow(condition)) {
          showModal(modalDialog(
            title = "Error Processing Input",
            "Possible missing data. Dimensions of Expression Data and Conditions File do not correspond."
          ))
          return(NULL)
        }
        if (is.null(condition$condition)) {
          showModal(modalDialog(
            title = "Error Processing Input",
            "\"condition\" column is missing or incorrectly named."
          ))
          return(NULL)
        }
        if (ncol(matrix) <= 1) {
          showModal(modalDialog(
            title = "Error Processing Input",
            "Expression file has too few rows or is incorrectly formatted. "
          ))
          return(NULL)
        }
        
        
        #build elist objects
        incProgress(1 / 6, detail = 'Building expression list objects')
        ##make genes list
        elistraw.genes <-
          data.frame(factor(row.names(matrix)), factor(row.names(matrix)))
        colnames(elistraw.genes) <-
          c("ProbeName", "SystematicName")
        
        ##make targets file
        elistraw.targets <-
          data.frame(colnames(matrix), make.names(condition$condition))
        cn <- c("FileName", "Condition")
        
        if (!is.null(condition$group)) {
          elistraw.targets$group <- condition$group
          cn <- c(cn, "group")
        }
        if (!is.null(condition$normalizer)) {
          elistraw.targets$normalizer <- condition$normalizer
          cn <- c(cn, "normalizer")
        }
        if (!is.null(condition$batch)) {
          elistraw.targets$batch <- condition$batch
          cn <- c(cn, "batch")
        }
        
        colnames(elistraw.targets) <- cn
        
        ##make E
        elistraw.e <- data.matrix(matrix)
        elistraw.list <- list()
        
        incProgress(1 / 3, detail = 'Generating new expression list')
        elistraw <- new("EListRaw", elistraw.list)
        elistraw$E <- elistraw.e
        elistraw$genes <- elistraw.genes
        elistraw$targets <- elistraw.targets
        
        showModal(modalDialog(
          title = "Data Read Successfully",
          "Uploaded files were successfully read. Proceed to Normalization or reupload a different set of data."
        ))
        return(elistraw)
      })
    } else{
      return(NULL)
    }
  })
  
  
  processedElist <-
    eventReactive(input$processData, ignoreNULL = FALSE, {
      if (is.null(rawElist()) || input$processData == 0) {
        return(NULL)
      } else {
        withProgress(message = 'Processing data', value = 0, {
          #get cleaned raw data
          rawEL <- rawElist()
          
          #subset
          if (!is.null(rawEL$targets$group)) {
            rawEL <- rawEL[, subsetVector()]
          }
          
          #mark low values
          incProgress(0.25, detail = "Marking low value genes")
          if (input$filtThresh == "Global Mean") {
            filterVector <-
              filterProbes(rawEL,
                           threshVal = "Global Mean",
                           percentSamples = (input$filtPercent %% 100))
          } else {
            filterVector <-
              filterProbes(
                rawEL,
                threshVal = input$filtThreshVal,
                percentSamples = (input$filtPercent %% 100)
              )
          }
          
          #normalize / check voom
          incProgress(0.25, detail = "Applying normalization and evaluating batches")
          
          if (input$norm == "DESeq") {
            row.names(rawEL$E) <- rawEL$genes$SystematicName
            countData <- ceiling(rawEL$E)
            colData <-
              data.frame(as.character(rawEL$targets$Condition))
            
            colnames(colData) <- c("condition")
            row.names(colData) <- rawEL$targets$FileName
            
            dds <-
              DESeqDataSetFromMatrix(
                countData = countData,
                colData = colData,
                design = ~ condition
              )
            dds <- estimateSizeFactors(dds)
            
            normEL <-
              normalizeBetweenArrays(rawEL, method = "none")
            if (input$voom == FALSE) {
              normE <- counts(dds, normalized = TRUE)
              normE <- log2(normE + 1) #avoid log(0) error
              normEL$E <- normE
            } else {
              normE <- counts(dds, normalized = TRUE)
              v <- voom(normE+0.5, plot = FALSE)
              normEL$E <- v$E
            }
            
          } else if (input$norm == "TMM") {
            row.names(rawEL$E) <- rawEL$genes$SystematicName
            countData <- ceiling(rawEL$E)
            
            colData <- as.character(rawEL$targets$Condition)
            
            er <- DGEList(counts = countData, group = colData)
            er <- calcNormFactors(er)
            
            normEL <-
              normalizeBetweenArrays(rawEL, method = "none")
            normE <- cpm(er, normalized.lib.sizes = TRUE)
            
            if (input$voom == FALSE) {
              normE <- log2(normE + 1)
              normEL$E <- normE
            }
            else{
              v <- voom(normE+0.5, plot = FALSE)
              normEL$E <- v$E
            }
            
          } else if (input$norm == "Quantile" || input$norm == "CyclicLoess" || input$norm == "Scale" || input$norm == "None") {
            if (input$voom == FALSE) {
              rawEL$E <- rawEL$E + 1 #avoid log(0) error
              normEL <-
                normalizeBetweenArrays(rawEL, method = tolower(input$norm))
            } else{
              v <-
                voom(rawEL$E+0.5,
                     plot = FALSE,
                     normalize.method = tolower(input$norm))
              #create elist object from elistraw
              normEL <-
                normalizeBetweenArrays(rawEL, method = "none")
              #replace E with voom-generated E
              normEL$E <- v$E
            }
          } else if(input$norm == "VSN"){
            normEL <-
              normalizeBetweenArrays(rawEL, method = "none")
            normEL$E <- justvsn(rawEL$E)
          } else {
            #CPM, total count, UQ, housekeeping, median
            normEL <-
              normalizeBetweenArrays(rawEL, method = "none")
            #housekeeping gene
            if (input$norm == "Housekeeping Gene") {
              if (input$housekeep == "Custom vector (\"normalizer\")" && !is.null(rawEL$targets$normalizer)) {
                normEL$E <-
                  normalizeMatrix(input$norm,
                                  rawEL$E,
                                  rawEL$targets$normalizer,
                                  voom = input$voom)
              }
              else{
                normEL$E <-
                  normalizeMatrix(input$norm,
                                  rawEL$E,
                                  rawEL$E[input$housekeep,],
                                  voom = input$voom)
              }
            } else{
              normEL$E <-
                normalizeMatrix(input$norm, rawEL$E, voom = input$voom)
            }
          }
          
          storageValues$noFiltElist <- normEL
          
          #clear out prior instances of sva objects
          storageValues$sv <- NULL
          storageValues$svaE <- NULL
          
          #remove negative counts
          normEL$E[normEL$E < 0] <- 0
          
          if(!is.null(normEL$targets$batch)){
            if(input$batchCorrect == "sva"){
              mod1 <- model.matrix(~as.factor(normEL$targets$Condition));
              mod0 <- model.matrix(~1,data=normEL$targets); 
              normEL$E <- (2 ^ normEL$E) - 1
              if(input$svaN == "Auto"){
                svaobj <- svaseq(normEL$E, mod1, mod0);
              } else{
                svaobj <- svaseq(normEL$E, mod1, mod0, n.sv = as.integer(input$svaN));
              }
              
              
              storageValues$sv <- svaobj$sv
              storageValues$svaE <- cleanE(normEL$E, mod1, svaobj$sv)
              normEL$E <- log2(normEL$E+1)
            } else if(input$batchCorrect == "comBat"){
              normEL$E <- (2 ^ normEL$E) - 1
              mod.combat <- model.matrix(~1, data = normEL$targets)
              combatE <- ComBat(dat = normEL$E, batch = normEL$targets$batch, mod = mod.combat, par.prior = TRUE, prior.plots = FALSE)
              
              normEL$E <- log2(combatE+1)
            } else if(input$batchCorrect == "Targetted Recentering"){
              if(!is.null(input$trGrps) && !setequal(unique(rawElist()$targets$group), input$trGrps)){
                normEL$E <- (2 ^ normEL$E) - 1
                groups <- input$trGrps
                groupVector <- rawElist()$targets$group
                E <- normEL$E
                for(i in 1:length(groups)){
                  target <- groups[i]
                  inGrp <- groupVector %in% target
                  outGrp <- !(groupVector %in% groups)
                  for(j in 1:10){
                    diff <- mean(E[,inGrp]) - mean(E[,outGrp]) 
                    E[,inGrp] <- E[,inGrp] - diff
                    E[E < 0] <- 0
                  }
                }
                normEL$E <- log2(E+1)
              }
            }
          }
          
          #filter
          incProgress(0.25, detail = "Filtering out low count genes")
          normEL <- normEL[filterVector,]
          
          incProgress(0.25, detail = "Reformatting data")
          Sys.sleep(0.1)
          showModal(modalDialog(
            title = "Data Processed Successfully",
            "Uploaded raw data was successfully normalized and filtered."
          ))
          return(normEL)
        })
      }
    })
  
  
  topTableList <- reactive({
    if(!(is.null(processedElist()) || is.null(input$condition1) || is.null(input$condition2))){
      if(length(intersect(input$condition1, input$condition2)) == 0){
        withProgress(message = 'Performing DE Analysis', value = 0, {
          incProgress(0.2, detail = "Getting expression values")
          normEL <- processedElist()
          
          
          
          #create contrast string
          incProgress(0.2, detail = "Making contrasts")
          
          condition2 <- input$condition2
          if(length(condition2)>1){
            l2 <- length(condition2)
            condition2 <- paste(condition2, collapse = "+")
            condition2 <- paste0("(", condition2, ")/", l2)
          }
          condition1 <- input$condition1
          if(length(condition1)>1){
            l1 <- length(condition1)
            condition1 <- paste(condition1, collapse = "+")
            condition1 <- paste0("(", condition1, ")/", l1)
          }
          
          contrast <-
            paste(condition2, "-", condition1, sep = "")
          
          #design model matrix
          incProgress(0.25, detail = "Building model matrix")
          f <- factor(normEL$targets$Condition, levels = unique(normEL$targets$Condition))
          design <- model.matrix(~0 + f)
          colnames(design) <- levels(f)
          if(!is.null(storageValues$sv)){
            design <- cbind(design, storageValues$sv)
            a <- c(levels(f), paste("Surrogate",seq_along(1:ncol(storageValues$sv)),sep = ""))
            colnames(design) <- a
          }
          
          #make lmfit
          incProgress(0.25, detail = "Applying linear model")
          fit <- lmFit(normEL, design)
          cont.matrix <- makeContrasts(contrasts = contrast, levels = design)
          secondfit <- contrasts.fit(fit, cont.matrix)
          secondfit <- eBayes(secondfit)
          
          #create top table for mir expression
          incProgress(0.1, detail = "Labelling DE features")
          tt <-
            topTable(
              secondfit,
              adjust = "BH",
              coef = 1,
              number = 3000,
              sort.by = "logFC"
            )
          return(tt)
        })
      }
      else{
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })
  
  #  #  #  #  #  #  #  #  #  #  #  #  #  #
  
  #List all uploaded files, with condition
  output$uploadStatusTable <- renderTable({
    input$upload
    input$condition
    input$files
    input$sampleData
    isolate({
      table <-
        cbind(c("Expression File", "Conditions File"),
              c("-", "-"),
              c("No", "No"),
              c("No", "No"),
              c("-", "-"),
              c("-", "-"))
      colnames(table) <-
        c("File Type", "File Name", "File Uploaded?", "File Read?", "Samples", "Features")
      if (!is.null(input$files$datapath)) {
        files <- input$files
        table[1, 2] <- toString(files$name[1])
        table[1, 3] <- "Yes"
        table[1, 4] <- "No"
      }
      if (!is.null(input$condition$datapath)) {
        condition <- input$condition
        table[2, 2] <- toString(condition$name[1])
        table[2, 3] <- "Yes"
        table[2, 4] <- "No"
      }
      if(input$sampleData){
        table[1, 2] <- "MCW_CRC_Exo_Matrix.csv"
        table[1, 3] <- "Yes"
        table[1, 4] <- "No"
        table[2, 2] <- "MCW_CRC_Exo_Conditions.csv"
        table[2, 3] <- "Yes"
        table[2, 4] <- "No"
      }
      if (!is.null(rawElist())) {
        table[1, 4] <- "Yes"
        table[2, 4] <- "Yes"
        table[1, 5] <- ncol(rawElist()$E)
        table[2, 5] <- nrow(rawElist()$targets)
        table[1, 6] <- nrow(rawElist()$E)
      }
      return(table)
    })
  })
  
  output$normStatusTable <-
    renderTable({
      input$processData
      vec <- rep("-", 7)
      m <- matrix(vec, 1, 7)
      colnames(m) <- c("Subgroups Processed", "Normalization Method", "Batch Correction", "Filtering Cutoff", "Filtering Sample Threshold", "miRNAs Remaining", "Used voom?")
      isolate(
        if (!is.null(processedElist())) {
          rEL <- rawElist()$E
          pEL <- processedElist()$E
          
          geneRemaining <-
            paste(nrow(pEL),
                  "remaining genes from",
                  nrow(rEL),
                  "initial genes.")
          
          if (input$filtThresh == "Global Mean") {
            filtCutoff <-
              paste0("Global Mean", " (", toString(sum(rEL) / (nrow(rEL) * ncol(rEL))), ")")
          } else{
            filtCutoff <- input$filtThreshVal
          }
          
          filtSamples <-
            paste0(toString(ceiling(input$filtPercent * ncol(pEL) / 100)), " (", input$filtPercent, "%)")
          
          if(input$batchCorrect == "sva"){
            batchCorrect <- paste("sva, ", ncol(storageValues$sv), "surrogates")
          }else if(input$batchCorrect == "Targetted Recentering"){
            if(!is.null(input$trGrps) && !setequal(unique(rawElist()$targets$group), input$trGrps)){
              batchCorrect <- paste("Targetted recentering: ", paste(input$trGrps, collapse = ", "))
            } else{
              batchCorrect <- "None"
            }
            
          } else{
            batchCorrect <- input$batchCorrect
          }
          
          m[1,1] <- toString(input$subsetGrps)
          m[1,2] <- input$norm
          m[1,3] <- batchCorrect
          m[1,4] <- filtCutoff
          m[1,5] <- filtSamples
          m[1,6] <- geneRemaining
          m[1,7] <- input$voom
        })
      return(m)
    }, include.colnames = TRUE, include.rownames = FALSE)
  
  output$upperQCPlot <- renderPlot({
    #dont plot until button has been pushed
    if (input$plotQC==0) {
      return(NULL)
    }
    isolate(
      if(!is.null(processedElist())) {
        pEL <- processedElist()
        rEL <- rawElist()
        
        #subset
        if(!is.null(pEL$targets$group)) {
          rEL <- rEL[, subsetVector()]
        }
        
        if(input$plotType == "PCA Plot") {
          withProgress(message = 'Generating PCA Plot', value = 0, {
            incProgress(0.25, detail = "Calculating PCs")
            matrix <- pEL$E
            
            pc <- prcomp(t(matrix))
            pc1 <- pc$x[, 1]
            pc2 <- pc$x[, 2]
            pc3 <- pc$x[, 3]
            
            #fix the labels/legends
            pVar <- round(100 * (pc$sdev) ^ 2 / sum(pc$sdev ^ 2))
            
            incProgress(0.25, detail = "Sorting sample information")
            if(input$sortPCA == "Condition"){
              colour <- factor(pEL$targets$Condition)
            } else if(input$sortPCA == "Group"){
              colour <- factor(pEL$targets$group)
            } else if(input$sortPCA == "Batch"){
              colour <- factor(pEL$targets$batch)
            }
            labels <- colnames(matrix)
    
            condition <- pEL$targets$Condition
            comps <- data.frame(labels, pc1, pc2, pc3, colour, condition)
            if(!is.null(pEL$targets$group)){
              comps$group <- pEL$targets$group
            }
            if(!is.null(pEL$targets$batch)){
              comps$batch <- pEL$targets$batch
            }
            
            incProgress(0.25, detail = "Generating base plot")
            p <- ggplot(data = comps, aes(x = pc1, y = pc2)) + geom_point(aes(colour = colour), size = 3, shape = 1) + geom_text(aes(label = labels), size = 3, hjust = 1.5) + labs(
              title = "Principal Component Analysis",
              x = paste("PC1 (", pVar[1], "%)"),
              y = paste("PC2 (", pVar[2], "%)"),
              colour = input$sortPCA
            ) 
            
            incProgress(0.25, detail = "Applying facets")
            if(input$pcaType == "PC1/PC2"){
              plotColor(p, fill = FALSE, input = input$plotColors, rev = input$colorsRev, n = length(unique(colour)))
            } else if(input$pcaType == "PC1/PC2 with K-means"){
              nClust <- switch(input$sortPCA, Batch = length(unique(pEL$targets$batch)), Group = length(unique(pEL$targets$Group)), Condition = length(unique(pEL$targets$Condition)))
              pcCluster <- kmeans(t(matrix), nClust, nstart = 100)
              comps$Cluster <- factor(pcCluster$cluster)
              
              p <- ggplot(data = comps, aes(x = pc1, y = pc2, shape = Cluster)) + geom_point(aes(colour = colour), size = 3) + geom_text(aes(label = labels), size = 3, hjust = 1.5) + labs(
                title = "Principal Component Analysis",
                x = paste("PC1 (", pVar[1], "%)"),
                y = paste("PC2 (", pVar[2], "%)"),
                colour = input$sortPCA
              ) 
              plotColor(p, fill = FALSE, input = input$plotColors, rev = input$colorsRev, n = length(unique(colour)))
            } else if(input$pcaType == "PC1/PC2/PC3"){
              p1 <- ggplot(data = comps, aes(x = pc1, y = pc2)) + geom_point(aes(colour = colour), size = 3, shape = 1) + labs(
                title = "PC1/PC2",
                x = paste("PC1 (", pVar[1], "%)"),
                y = paste("PC2 (", pVar[2], "%)"),
                colour = input$sortPCA
              ) 
              p2 <- ggplot(data = comps, aes(x = pc2, y = pc3)) + geom_point(aes(colour = colour), size = 3, shape = 1) + labs(
                title = "PC2/PC3",
                x = paste("PC2 (", pVar[2], "%)"),
                y = paste("PC3 (", pVar[3], "%)"),
                colour = input$sortPCA
              ) 
              p3 <- ggplot(data = comps, aes(x = pc1, y = pc3)) + geom_point(aes(colour = colour), size = 3, shape = 1) + labs(
                title = "PC1/PC3",
                x = paste("PC1 (", pVar[1], "%)"),
                y = paste("PC3 (", pVar[3], "%)"),
                colour = input$sortPCA
              ) 
              multiplot(p1, p2, p3, cols = 2)
            } else if(input$pcaType == "Facet by Condition"){
              plotColor(p + facet_wrap(~condition, scales = "free", nrow = 2), fill = FALSE, input = input$plotColors, rev = input$colorsRev, n = length(unique(colour)))
            } else if(input$pcaType == "Facet by Group"){
              plotColor(p + facet_wrap(~group, scales = "free", nrow = 2), fill = FALSE, input = input$plotColors, rev = input$colorsRev, n = length(unique(colour)))
            } else if(input$pcaType == "Facet by Batch"){
              plotColor(p + facet_wrap(~batch, scales = "free", nrow = 2), fill = FALSE, input = input$plotColors, rev = input$colorsRev, n = length(unique(colour)))
            }
          })
        } else if (input$plotType == "Correlation Coefficient Matrix") {
          withProgress(message = 'Generating Correlation Coefficient Matrix', value = 0, {
            matrix <- pEL$E
            
            incProgress(0.25, detail = "Sorting sample information")
            if(input$cctype == "By Group"){
              group <- pEL$targets$group
              
              u <- unique(group)
              for(i in 1:length(u)){
                if(i ==1){
                  means <- data.frame(rowMeans(matrix[,group==u[1]]))
                }
                else{
                  means <- cbind(means, rowMeans(matrix[,group==u[i]]))
                }
              }
              matrix <- as.matrix(means)
              colnames(matrix) <- u
            } else if(input$cctype == "By Batch"){
              batch <- pEL$targets$batch
              
              u <- unique(batch)
              for(i in 1:length(u)){
                if(i==1){
                  means <- data.frame(rowMeans(matrix[,batch==u[1]]))
                }
                else{
                  means <- cbind(means, rowMeans(matrix[,batch==u[i]]))
                }
              }
              matrix <- as.matrix(means)
              colnames(matrix) <- u
            }
            
            incProgress(0.25, detail = "Calculating correlation coefficients")
            corA <- cor(matrix, method = tolower(input$CCcalcType))
            cor.m <- melt(corA)
            
            incProgress(0.25, detail = "Adjusting matrix")
            if(input$cctype == "By Group"){
              cor.m$X1 <- rep(unique(group), ncol(matrix))
              
              X2 <- numeric()
              
              for(i in 1:ncol(matrix)){
                X2 <- c(X2, rep(unique(group)[i], ncol(matrix)))
              }
              cor.m$X2 <- X2
            }
            else if(input$cctype == "By Batch"){
              cor.m$X1 <- rep(unique(batch), ncol(matrix))
              
              X2 <- character()
              
              for(i in 1:ncol(matrix)){
                X2 <- c(X2, rep(unique(batch)[i], ncol(matrix)))
              }
              cor.m$X2 <- X2
            }
            else{
              cor.m$X1 <- rep(1:ncol(matrix), ncol(matrix))
              
              X2 <- numeric()
              
              for(i in 1:ncol(matrix)){
                X2 <- c(X2, rep(i, ncol(matrix)))
              }
              cor.m$X2 <- X2
            }
            incProgress(0.25, detail = "Building plot")
            if(input$cctype == "By Group"){
              p <- ggplot(cor.m, aes(x = factor(X2), y = factor(X1))) + geom_tile(aes(fill = value), colour = "white") + labs(x = "Group", y = "Group") 
            } else if(input$cctype == "By Batch"){
              p <- ggplot(cor.m, aes(x = factor(X2), y = factor(X1))) + geom_tile(aes(fill = value), colour = "white") + labs(x = "Batch", y = "Batch") 
            } else{
              p <- ggplot(cor.m, aes(x = X2, y = X1)) + geom_tile(aes(fill = value), colour = "white") + labs(x = "Sample Number", y = "Sample Number") 
            }
            plotColor(p, discrete = FALSE, input = input$plotColors, rev = input$colorsRev)
          })
        } else if (input$plotType == "Boxplot") {
          E <- rEL$E
          withProgress(message = 'Generating raw boxplots', value = 0, {
            t <- "Raw Expression Distribution"
            if(input$sortBox == 'Sample Number'){
              lab <- "Sample Number"
              p <- plotOrderBoxplot(E = E, title = t, lab = lab, cont = TRUE)
              print(p)
            } else{
              if(input$sortBox == 'Condition'){
                lab <- "Condition"
                ov <- rEL$targets$Condition
              } else if(input$sortBox == 'Group'){
                lab <- "Group"
                ov <- rEL$targets$group
              } else if(input$sortBox == 'Batch'){
                lab <- "Batch"
                ov <- rEL$targets$batch
              }
              p <- plotOrderBoxplot(orderVector = ov, E = E, title = t, lab = lab)
              plotColor(p, n = length(unique(ov)), input = input$plotColors, rev = input$colorsRev)
            }
          })
        } else{
          return(NULL)
        }
      } else{
        return(NULL)
      })
  })
  
  #render plot of mirnas post processing
  output$lowerQCPlot <- renderPlot({
    #dont plot until button has been pushed
    if (input$plotQC==0) {
      return(NULL)
    }
    isolate(if (!is.null(processedElist())) {
      pEL <- processedElist()
      
      if(!is.null(storageValues$svaE)){
        pEL$E <- storageValues$svaE
      }
    
      if (input$plotType == "PCA Plot") {
        return(NULL)
      } else if (input$plotType == "Correlation Coefficient Matrix") {
        return(NULL)
      } else if (input$plotType == "Boxplot") {
        withProgress(message = 'Generating normalized boxplots', value = 0, {
          t <- "Normalized Expression Distribution"
          E <- pEL$E
          if(input$sortBox == 'Sample Number'){
            lab <- "Sample Number"
            p <- plotOrderBoxplot(E = E, title = t, lab = lab, cont = TRUE)
            print(p)
          } else {
            if(input$sortBox == 'Condition'){
              lab <- "Condition"
              ov <- pEL$targets$Condition
            } else if(input$sortBox == 'Group'){
              lab <- "Group"
              ov <- pEL$targets$group
            } else if(input$sortBox == 'Batch'){
              lab <- "Batch"
              ov <- pEL$targets$batch
            }
            p <- plotOrderBoxplot(orderVector = ov, E = E, title = t, lab = lab)
            plotColor(p, n = length(unique(ov)), input = input$plotColors,  rev = input$colorsRev)
          }  
        })
      } else{
        return(NULL)
      }
    } else{
      return(NULL)
    })
  })
  
  output$MAPlot <- renderPlot({
    if(!is.null(topTableList())){
      if(is.null(input$condition1)){
        return(NULL)
      }
      if(is.null(input$condition2)){
        return(NULL)
      }
      withProgress(message = 'Generating MA plot', value = 0, {
        incProgress(0.333, detail = "Getting DE values")
        mirlist <- topTableList()
        incProgress(0.333, detail = "Generating thresholds")
        mirlist$threshold = as.factor(mirlist$AveExpr > input$AveEcut)
        
        incProgress(0.333, detail = "Plotting MA Values")
        MAplot <-
          ggplot(mirlist, aes(
            x = AveExpr,
            y = logFC,
            colour = threshold
          )) + geom_point(alpha = 0.75, size = 2, shape = 1) +
          theme(legend.position = "none") +
          ylab("M (Log2FC)") + xlab("A (Average Expression)") +
          ggtitle("M vs A (Log2FC vs. Average Expression)") +
          scale_color_manual(values = c("#8c8c8c", "#ff0000")) +
          geom_vline(xintercept = input$AveEcut, color = "Red")
      })
      return(MAplot)
    }
    else{
      return(NULL)
    }
  })
  
  #create volcano plot
  output$volcanoPlot <- renderPlot({
    if(!is.null(topTableList())){
      if(is.null(input$condition1)){
        return(NULL)
      }
      if(is.null(input$condition2)){
        return(NULL)
      }
      withProgress(message = 'Generating volcano plot', value = 0, {
        incProgress(0.333, detail = "Getting DE values")
        mirlist <- topTableList()
        incProgress(0.333, detail = "Generating thresholds")
        mirlist$threshold = as.factor(abs(mirlist$logFC) > log2(input$logFCcut) & mirlist$P.Value < 0.05 & mirlist$AveExpr > input$AveEcut)
        
        incProgress(0.333, detail = "Plotting features")
        vplot <-
          ggplot(mirlist, aes(
            x = logFC,
            y = -log10(P.Value),
            colour = threshold
          )) + geom_point(alpha = 0.75, size = 2, shape = 1) +
          theme(legend.position = "none") +
          xlim(c(-5, 5)) + ylim(c(0, 15)) +
          xlab("log2(fold change)") + ylab("-log10(p-value)") +
          ggtitle("P-Value vs Log2FC") +
          scale_color_manual(values = c("#8c8c8c", "#ff0000")) +
          geom_vline(xintercept = log2(input$logFCcut), color = "Red") +
          geom_vline(xintercept = -log2(input$logFCcut), color = "Red") +
          geom_abline(intercept = -log10(input$pValCut), slope = 0, color = "Red") 
      })
      return(vplot)
    }
    else{
      return(NULL)
    }
  })
  
  output$heatmap <- renderPlot({
    if(!is.null(topTableList())){
      
      if(is.null(input$condition1)){
        return(NULL)
      }
      if(is.null(input$condition2)){
        return(NULL)
      }
      
      withProgress(message = 'Generating heatmap', value = 0, {
        incProgress(0.25, detail = "Getting DE values")
        mirlist <- topTableList()
        elist <- processedElist()
        incProgress(0.2, detail = "Removing low values")
        #sig represents which elements of the list are eligible
        sig <-
          (mirlist$P.Value < input$pValCut &
             abs(mirlist$logFC) > log2(input$logFCcut) & mirlist$AveExpr > input$AveEcut)
        
        storageValues$sig <- sig
        
        
        sig_mir_list <- mirlist[sig,]
        storageValues$sigMirs <- sig_mir_list$SystematicName
        
        
        heatmap_match <-
          match(elist$genes$ProbeName, sig_mir_list$ProbeName, nomatch = NA)
        heatmap_elist <- elist[!is.na(heatmap_match),]
        heatmap_elist = heatmap_elist[order(heatmap_match[!is.na(heatmap_match)]),] #sorts by the table
        heatmap_elist = heatmap_elist[,order(heatmap_elist$targets$Condition)] #arranges columns by condition
        incProgress(0.2, detail = "Applying row transform")
        topMatrix <- t(scale(t(heatmap_elist$E)))
        topMatrix[topMatrix < -3] <- -3
        topMatrix[topMatrix > 3] <- 3
        storageValues$topMatrix <- topMatrix
        
        #insufficient DE features
        if(nrow(topMatrix) * ncol(topMatrix) < 4){
          return(NULL)
        }
        
        rownames(topMatrix) <- heatmap_elist$genes$SystematicName
        colnames(topMatrix) <- heatmap_elist$targets$FileName
        
        storageValues$hmMatrix <- topMatrix
        
        Condition <- heatmap_elist$targets$Condition
        annotation <- data.frame(Condition)
        
        storageValues$annotation <- Condition
        
        incProgress(0.2, detail = "Building heatmap")
        if(input$dendClust==0){
          rv <- NA
          cv <- NA
        } else if(input$dendClust==1){
          rv <- TRUE
          cv <- NA
        } else{
          rv <- TRUE
          cv <- TRUE
        }
        a <-
          aheatmap(
            topMatrix,
            col = heatmapPal(input$hmColors, input$hmRev),
            Rowv = rv,
            Colv = cv,
            scale = "none",
            annCol = annotation,
            annColors = list(rainbow(length(unique(Condition))))
          )
        incProgress(0.2, detail = "Applying annotation")
        storageValues$rowInd <- a$rowInd
        storageValues$colInd <- a$colInd
      })
      return(a)
    }else{
      return(NULL)
    }
    
  })
  
  output$sigMirTable <- renderTable({
    if(!is.null(topTableList())){
      mirlist <- topTableList()
      sig <-
        (mirlist$P.Value < (input$pValCut) &
           abs(mirlist$logFC) > log2(input$logFCcut) & mirlist$AveExpr > input$AveEcut)
      reg <-
        data.frame(mirlist[sig, c("SystematicName","logFC","AveExpr","t","P.Value","adj.P.Val","B")])
      colnames(reg) <-
        c(            
          "Systematic Name",
          "Log2 Fold Change",
          "Average Expression",
          "t-value",
          "P-value",
          "Adjusted P-value (q)",            
          "Log Odds (B)"
        )
      return(reg)
    } else{
      return(NULL)
    }
  }, digits = 5, caption = "<b><h4>Differentially Expressed miRNAs</h4></b>", caption.placement = getOption("xtable.caption.placement", "top"),
  caption.width = getOption("xtable.caption.width", NULL))
  
  heatmapFile <- reactive({
    E <- storageValues$hmMatrix
    Condition <- storageValues$annotation
    anno <- data.frame(Condition)
    
    
    if(is.null(E) || is.null(anno)){
      return(NULL)
    }
    if(nrow(E) * ncol(E) < 4){
      return(NULL)
    }
    
    if(input$dendClust==0){
      rv <- NA
      cv <- NA
    } else if(input$dendClust==1){
      rv <- NA
      cv <- TRUE
    } else{
      rv <- TRUE
      cv <- TRUE
    }
    a <-
      aheatmap(
        E,
        col = heatmapPal(input$hmColors, input$hmRev),
        Rowv = rv,
        Colv = cv,
        scale = "none",
        annCol = anno,
        annColors = list(rainbow(length(unique(Condition))))
      )
    return(a)
  })
  
  heatmapMatrixFile <- reactive({
    E <- storageValues$topMatrix
    rowInd <- storageValues$rowInd
    colInd <- storageValues$colInd
    if(is.null(rowInd)){
      rowInd <- TRUE
    }
    if(is.null(colInd)){
      colInd <- TRUE
    }
    heatmapMatrixFile <- t(scale(t(E[rev(rowInd),colInd])))
    return(heatmapMatrixFile)
  })
  
  DeTableFile <- reactive({
    mirlist <- topTableList()
    sig <-
      (mirlist$P.Value < (input$pValCut) &
         abs(mirlist$logFC) > log2(input$logFCcut) & mirlist$AveExpr > input$AveEcut)
    reg <-
      data.frame(
        mirlist[sig, c("SystematicName","logFC","AveExpr","t","P.Value","adj.P.Val","B")])
    colnames(reg) <-
      c(
        "Systematic Name",
        "Log2 Fold Change",
        "Average Expression",
        "t-value",
        "P-value",
        "Adjusted P-value (q)",
        "Log Odds (B)"
      )
    row.names(reg) <- NULL
    return(reg)
  })
  FullTableFile <- reactive({
    normEL <- storageValues$noFiltElist
    
    #create contrast string
    contrast <-
      paste(input$condition2, "-", input$condition1, sep = "")
    
    #design model matrix
    f <- factor(normEL$targets$Condition, levels = unique(normEL$targets$Condition))
    design <- model.matrix(~0 + f)
    colnames(design) <- levels(f)
    
    #make lmfit
    fit <- lmFit(normEL, design)
    cont.matrix <- makeContrasts(contrasts = contrast, levels = design)
    secondfit <- contrasts.fit(fit, cont.matrix)
    secondfit <- eBayes(secondfit)
    
    #create top table for mir expression
    mirlist <-
      topTable(
        secondfit,
        adjust = "BH",
        coef = 1,
        number = 10000,
        sort.by = "logFC"
      )
    reg <-
      data.frame(
        mirlist[, c("SystematicName","logFC","AveExpr","t","P.Value","adj.P.Val","B")])
    
    colnames(reg) <-
      c(
        "Systematic Name",
        "Log2 Fold Change",
        "Average Expression",
        "t-value",
        "P-value",
        "Adjusted P-value (q)",
        "Log Odds (B)"
      )
    row.names(reg) <- NULL
    return(reg)
  })
  
  #  #  #  #  #  #  #  #  #  #  #  #  #  #
  
  output$subsetUI <- renderUI({
    if (!is.null(rawElist()$targets$group)) {
      groups <- unique(rawElist()$targets$group)
      ordGroups <- groups[order(groups)]
      selectInput(
        inputId = "subsetGrps",
        label = NULL,
        multiple = TRUE,
        choices = ordGroups,
        selected = ordGroups
      )
    }
    else{
      return(NULL)
    }
  })
  
  subsetVector <- reactive({
    if (!is.null(rawElist()$targets$group)) {
      groupVector <- rawElist()$targets$group
      subVector <- groupVector %in% input$subsetGrps
      if(sum(subVector) == 0) subVector <- rep(TRUE, length(groupVector))
      return(subVector)
    }
    else{
      return(NULL)
    }
  })
  
  output$trGrpsUI <- renderUI({
    if (!is.null(rawElist()$targets$group)) {
      groups <- unique(rawElist()$targets$group)
      ordGroups <- groups[order(groups)]
      selectInput(
        inputId = "trGrps",
        label = "Groups To Recenter",
        multiple = TRUE,
        choices = ordGroups,
        selected = NULL
      )
    }
    else{
      return(NULL)
    }
  })
  
  output$housekeepUI <- renderUI({
    if (!is.null(rawElist())) {
      hk <- rowSums(rawElist()$E)
      names(hk) <- row.names(rawElist()$E)
      hk <-
        hk[order(hk, decreasing = TRUE)]
      
      #pick top 20 high expressed microRNAs
      hk <- hk[1:20]
      hk <- names(hk)
      
      if (!is.null(rawElist()$targets$normalizer)){
        hk <-
          c(hk, "Custom vector (\"normalizer\")")
      }
      selectInput(inputId = "housekeep",
                  label = "Housekeeping Gene (ranked by expression)",
                  choices = hk)
    }
    else{
      return(NULL)
    }
  })
  
  output$filtPercent <- renderUI({
    rEL <- rawElist()
    if (!is.null(rEL$targets$group)) {
      rEL <- rEL[, subsetVector()]
    }
    val <- 50
    if (!is.null(rawElist())) {
      #find 1 sample value
      val <- round(100 / (ncol(rEL$E) + 1), 2)
    }
    fp <- numericInput(
      inputId = "filtPercent",
      label = "Sample % Threshold",
      min = 0,
      max = 100,
      value = val,
      step = 0.1
    )
    
    return(fp)
  })
  
  output$filtSampleNumber <- renderText({
    if (!is.null(rawElist())) {
      rEL <- rawElist()
      if (!is.null(rEL$targets$group)) {
        rEL <- rEL[, subsetVector()]
      }
      print(ceiling(ncol(rEL$E) * (input$filtPercent %% 100) / 100)) #this print statement is important for some reason
    }
    else{
      print("-")
    }
  })
  
  output$upperPlotUI <- renderUI({
    input$plotQC
    #defaults
    w <- 1000
    h <- 400
    isolate(
      if(input$plotType == 'PCA Plot'){
        w <- 950
        h <- 1000
      }
      else if(input$plotType == 'Correlation Coefficient Matrix'){
        w <- 950
        h <- 950
      }
      else if(input$plotType == 'Boxplot'){
        w <- '100%'
        h <- 450
      }
    )
    plotOutput("upperQCPlot", width = w, height = h)
  })
  
  output$lowerPlotUI <- renderUI({
    input$plotQC
    w <- 1000
    h <- 400
    isolate(
      if(input$plotType == 'PCA'){
        return(NULL)
      }
      else if(input$plotType == 'Correlation Coefficient Matrix'){
        return(NULL)
      }
      else if(input$plotType == 'Boxplot'){
        w <- '100%'
        h <- 450
      }
    )
    plotOutput("lowerQCPlot", width = w, height = h)
  })
  
  output$sortBoxplotUI <- renderUI({
    rEL <- rawElist()
    choices <- c("Sample Number", "Condition")
    if(!is.null(rEL$targets$group)){
      choices <- c(choices, "Group")
    }
    if(!is.null(rEL$targets$batch)){
      choices <- c(choices, "Batch")
    }
    selectInput(
      inputId = "sortBox",
      label = "Order Samples by:",
      choices = choices,
      selected = choices[1]
    )
  })
  
  output$PCAtypeUI <- renderUI({
    rEL <- rawElist()
    choices <- c("PC1/PC2", "PC1/PC2 with K-means", "PC1/PC2/PC3", "Facet by Condition")
    if(!is.null(rEL$targets$group)){
      choices <- c(choices, "Facet by Group")
    }
    if(!is.null(rEL$targets$batch)){
      choices <- c(choices, "Facet by Batch")
    }
    selectInput(
      inputId = "pcaType",
      label = "PCA Plot Display Options",
      choices = choices,
      selected = choices[1]
    )
  })
  
  output$sortPCAUI <- renderUI({
    rEL <- rawElist()
    choices <- c("Condition")
    if(!is.null(rEL$targets$group)){
      choices <- c(choices, "Group")
    }
    if(!is.null(rEL$targets$batch)){
      choices <- c(choices, "Batch")
    }
    selectInput(
      inputId = "sortPCA",
      label = "PCA Plot Grouping Options",
      choices = choices,
      selected = choices[1]
    )
  })
  
  output$CCtypeUI <- renderUI({
    rEL <- rawElist()
    choices <- c("By Sample")
    if(!is.null(rEL$targets$group)){
      choices <- c(choices, "By Group")
    }
    if(!is.null(rEL$targets$batch)){
      choices <- c(choices, "By Batch")
    }
    selectInput(
      inputId = "cctype",
      label = "CC Matrix Display Options",
      choices = choices,
      selected = choices[1]
    )
  })
  
  output$condition1UI <- renderUI({
    rEL <- rawElist()
    c <- unique(sapply(X = rEL$targets$Condition, FUN = toString))
    selectInput(
      inputId = "condition1",
      multiple = TRUE,
      label = "Control group",
      choices = c,
      selected = c[1]
    )
  })
  
  output$condition2UI <- renderUI({
    rEL <- rawElist()
    c <- unique(sapply(X = rEL$targets$Condition, FUN = toString))
    selectInput(
      inputId = "condition2",
      multiple = TRUE,
      label = "Patient group",
      choices = c,
      selected = c[2]
    )
  })
  
  output$heatmapUI <- renderUI({
    if(!is.null(processedElist())&&!is.null(topTableList())){
      
      mirlist <- topTableList()
      elist <- processedElist()
      sig <-
        (mirlist$P.Value < input$pValCut &
           abs(mirlist$logFC) > log2(input$logFCcut) & mirlist$AveExpr > input$AveEcut)
      
      sig_mirlist <- mirlist[sig,]
      
      col <- ncol(elist$E)
      row <- nrow(sig_mirlist)
      
      w <- '100%'
      h <- 700
      
      ratio <- col / row
      if(ratio < 4) ratio <- 4
      if(ratio > 16) ratio <- 16
      ratio <- log2(ratio)
      
      h <- 700 - 175 * (ratio-2)
      
      if(col <30){
        scale <- 0.4 + 0.02 * col
        h <- scale * h
        w <- paste0(100 * scale, '%')
      }
      storageValues$hmW <- w
      storageValues$hmH <- h
      return(plotOutput("heatmap", width = w, height = h))
    }
    else{
      return(NULL)
    }
  })
  output$downloadHM <- downloadHandler(
    filename = function() {
      paste(input$hmImageName, '.png', sep = '')
    },
    content = function(file) {
      png(file, width = 3.5*storageValues$hmH, height = 2*storageValues$hmH)
      print(heatmapFile())
      dev.off()
    }
  )
  output$downloadMatrix <- downloadHandler(
    filename = function() {
      paste(input$hmName, '.csv', sep = '')
    },
    content = function(file) {
      write.csv(heatmapMatrixFile(), file)
    }
  )
  output$downloadDETab <- downloadHandler(
    filename = function() {
      paste(input$deName, '.csv', sep = '')
    },
    content = function(file) {
      write.csv(DeTableFile(), file)
    }
  )
  output$downloadFullTab <- downloadHandler(
    filename = function() {
      paste(input$fullName, '.csv', sep = '')
    },
    content = function(file) {
      write.csv(FullTableFile(), file)
    }
  )
  output$dataExists <- reactive({
    return(!is.null(rawElist()))
  })
  outputOptions(output, 'dataExists', suspendWhenHidden = FALSE)
  output$groupExists <- reactive({
    return(!is.null(rawElist()$targets$group))
  })
  outputOptions(output, 'groupExists', suspendWhenHidden = FALSE)
  
  output$batchExists <- reactive({
    return(!is.null(rawElist()$targets$batch))
  })
  outputOptions(output, 'batchExists', suspendWhenHidden = FALSE)
  output$sigMirsExist <- reactive({
    return(!is.null(storageValues$sigMirs) && length(storageValues$sigMirs) > 0)
  })
  outputOptions(output, 'sigMirsExist', suspendWhenHidden = FALSE)
  
  output$deSelectUI <- renderUI({
    sigMirs <- storageValues$sigMirs
    selectInput(
      inputId = "sigMir",
      label = "Select miRNAs to view",
      choices = sigMirs,
      selected = sigMirs[1],
      multiple = TRUE
    )
  })
  output$singleMirBoxplot <- renderPlot({
    input$plotDE
    isolate(
      if(!is.null(input$sigMir)){
        withProgress(message = 'Generating boxplots', value = 0, {
          selected <- input$sigMir
          pEL <- processedElist()
          incProgress(0.1, detail = "Getting expression values")
          
          condition <- rep(pEL$targets$Condition, length(selected))
          
          incProgress(0.1, detail = "Subsetting dataset")
          values <- vector()
          for(i in 1: length(selected)){
            values <- c(values, pEL$E[selected[i],])
          }
          
          sample <- vector()
          for(i in 1:length(selected)){
            sample <- c(sample, rep(selected[i], ncol(pEL)))
          }
          
          incProgress(0.4, detail = "Building plot")
          bpDF <- data.frame(condition, values, sample)
          if(input$dePlotType == "Boxplot"){
            p <- ggplot(bpDF, aes(x = condition, y = values, fill = factor(condition))) + geom_boxplot(outlier.shape = 3, outlier.color = '#ff0000', width = 0.65, outlier.size = 1, outlier.alpha = 0.2, notch = input$deNotch) + xlab("Condition") + ylab("Normalized Expression Value") + labs(fill = "Condition")
            if(input$deDot){
              p <- p + geom_dotplot(binaxis = "y", stackdir = "centerwhole", fill = "White", dotsize = 0.65)
            }
          } else {
            #Violin Plot
            p <- ggplot(bpDF, aes(x = condition, y = values, fill = factor(condition))) + geom_violin() + xlab("Condition") + ylab("Normalized Expression Value") + labs(fill = "Condition")
            if(input$deBox){
              p <- p + geom_boxplot(width = 0.1, fill = "White")
            }
          }
          incProgress(0.4, detail = "Applying facets ")
          plotColor(p+facet_wrap(~sample, scales = "free", ncol = ceiling(sqrt(length(selected)))) + theme(axis.title=element_text(size=14), axis.text=element_text(size=12)), input = input$miRplotColors, rev = input$miRcolorsRev, n = length(unique(condition)))
        })
      }
      else{
        return(NULL)
      }
    )
  })
  
  output$idealProfileTable <- renderTable({
    input$idealButton
    isolate(
      if(!is.null(input$idealProfile) && !is.null(rawElist())){
        withProgress(message = "Calculating correlation", value = 0, {
          incProgress(0.25, detail = "Reading values")
          idealProfile <- read.csv(file = win(input$idealProfile$datapath), header = TRUE, check.names = FALSE)
          if(ncol(idealProfile > 5)){
            idealProfile <- idealProfile[,1:5]
          }
          incProgress(0.25, detail = "Cleaning data")
          idealProfile[is.na(idealProfile)] <- 0
          idealProfile[is.null(idealProfile)] <- 0
          idealProfile[idealProfile < 0] <- 0
          fullE <- rawElist()$E
          if(nrow(fullE) == nrow(idealProfile)){
            incProgress(0.25, detail = "Calculating values")
            cor.m <- cor(fullE, idealProfile, method = tolower(input$idealCorrType))
            avg <- colMeans(cor.m)
            rms <- colSums(((cor.m - 1) ^ 2))
            m <- rbind(avg, rms, cor.m)
            incProgress(0.25, detail = "Binding data")
            sNames <- c("Average", "RMS.Error", colnames(fullE))
            colnames(m) <- colnames(idealProfile)
            row.names(m) <- sNames
            return(m)
          } else{
            return(NULL)
          }
        })
      } else {
        return(NULL)
      }
    )
  }, rownames = TRUE)
  
  observeEvent(input$powerPlotButton, {
    
    if(is.null(storageValues$powerCurves)) {
      storageValues$powerCurves = list(5)
      storageValues$numPowerCurves = 0
    }
    if(is.null(storageValues$powerTables)) {
      storageValues$powerTables = list(5)
    }
    
    withProgress(message = "Setting parameters", value = 0, {
      #default values
      samplesize = input$powerSampleSize
      foldchange = input$powerFoldChange
      fdr = input$powerFDR
      averagereadcounts = input$powerAverageReadCounts
      maxdispersion = input$powerMaxDispersion
      featurecount = input$powerFeatureCount
      prognosticfeaturecount = input$powerPrognosticFeatureCount
      dataMatrixDistribution = NULL
      repNumber = 1
      #if user wants, we set above variables based on calculations of uploaded data
      if(input$powerUseData) {
        setProgress(0.2, message = "Calculating other parameters")
        featurecount = nrow(rawElist())
        if(input$powerUseCutoff) { #use differentially expressed mirna
          mirList = topTableList()
          sigmirs = mirList[storageValues$sig,]
          prognosticfeaturecount = nrow(sigmirs)
          foldchange = 2 ^ min(abs(sigmirs$logFC))
          fdr = max(sigmirs$adj.P.Val)
        }
        #get the distribution object from the raw elist, even if user uses DE mirna
        #produces an NA if there is no replication
        dataMatrixDistribution = est_count_dispersion(rawElist()$E,
          group = rawElist()$targets$group, minAveCount = 0, subSampleNum = 100)
        if(!is.na(dataMatrixDistribution$common.dispersion)) {
          #only used to print out parameters in the table
          maxdispersion = "NA"
        } else {
          showModal(modalDialog(
            title = "Error Calculating Parameters",
            "Uploaded data has no replications, cannot make statistical inferences from data."
          ))
          return(NULL)
        }
        #only used to print out parameters in the table
        averagereadcounts = "NA"
        #only use data points with read counts >= 5
        if(length(which(dataMatrixDistribution$pseudo.counts.mean >= 5)) < input$powerRepNumber) {
          showModal(modalDialog(
            title = "Warning Calculating Parameters",
            paste("Estimation repetition count too high for this data, lowering to", length(which(dataMatrixDistribution$pseudo.counts.mean >= 5)))
          ))
          repNumber = length(which(dataMatrixDistribution$pseudo.counts.mean >= 5))
        } else {
          repNumber = input$powerRepNumber
        }
      }
      if(input$powerPlotType == "Sample size"){
        setProgress(0.5, message = "Calculating power")
        #can only display 5
        storageValues$powerCurves[[(storageValues$numPowerCurves %% 5) + 1]] = est_power_curve_better(
          n = samplesize, f = fdr, rho = foldchange, lambda0 = averagereadcounts, phi0 = maxdispersion, m = featurecount, m1 = prognosticfeaturecount,
          stepsize = input$powerStepInterval, squarestep = (input$powerPlotChoice == "Squared intervals"), usegradient = (input$powerPlotChoice == "Gradient-sensitive intervals"), gradientdetail = input$powerGradientDetail,
          distributionObject = dataMatrixDistribution,repNumber = repNumber)
        #store these paramters to display later
        storageValues$powerCurves[[(storageValues$numPowerCurves %% 5) + 1]]$parameters["m"] = featurecount
        storageValues$powerCurves[[(storageValues$numPowerCurves %% 5) + 1]]$parameters["m1"] = prognosticfeaturecount
      } else if(input$powerPlotType == "Power"){
        setProgress(0.5, message = "Calculating sample size")
        desiredpower = input$powerMaxPower
        if(desiredpower >= 1) { #power can never be 1
          desiredpower = 0.999
        }
        if(foldchange == 1) { #quik n dirty error handling
          showModal(modalDialog(
            title = "Error Calculating Sample Size",
            "No non-zero power for a fold change of 1"
          ))
          return(NULL)
        }
        if(is.null(dataMatrixDistribution)) {
          storageValues$powerCurves[[(storageValues$numPowerCurves %% 5) + 1]] = sample_size(power = desiredpower, f = fdr, rho = foldchange, lambda0 = averagereadcounts, phi0 = maxdispersion, m = featurecount, m1 = prognosticfeaturecount, storeProcess = T)
        } else {
          storageValues$powerCurves[[(storageValues$numPowerCurves %% 5) + 1]] = sample_size_distribution(power = desiredpower, f = fdr, rho = foldchange, distributionObject = dataMatrixDistribution, repNumber = repNumber, m = featurecount, m1 = prognosticfeaturecount, storeProcess = T)
        }
      }
      #debug print
      #print(storageValues$powerCurves[[(storageValues$numPowerCurves %% 5) + 1]])
      setProgress(1, message = "Plotting")
      #displays a table of #samples vs power
      storageValues$powerTables[[(storageValues$numPowerCurves %% 5) + 1]] = 
      if(!is.null(storageValues$powerCurves[[(storageValues$numPowerCurves %% 5) + 1]])){
        curve = storageValues$powerCurves[[(storageValues$numPowerCurves %% 5) + 1]]
        reg = data.frame(curve$process[, c("N","Power")])
        colnames(reg) =
          c(            
            "Samples",
            "Power"
          )
        reg
      } else{
        NULL
      }
      storageValues$numPowerCurves = storageValues$numPowerCurves + 1
    })
  })
  observeEvent(input$powerClearPlotButton, {
    storageValues$powerCurves = NULL
    storageValues$powerTables = NULL
    storageValues$numPowerCurves = 0
  })
  
  output$powerPlot = renderPlot({
    if(!is.null(storageValues$powerCurves)) {
      plot_power_curve_better(storageValues$powerCurves, las = 2, main = "Sample Size vs. Statistical Power")
    }
  })
  #displays parameters of each graph
  output$powerParameterTable = renderTable({
    if(!is.null(storageValues$powerCurves)){
      reg = data.frame(storageValues$powerCurves[[1]]$parameters[c("fdr","rho","lambda0","phi0", "m", "m1")])
      if(length(storageValues$powerCurves) > 1) {
        for(x in 2:length(storageValues$powerCurves)) {
          reg = data.frame(reg, storageValues$powerCurves[[x]]$parameters[c("fdr","rho","lambda0","phi0", "m", "m1")])
        }
      }
      
      rownames(reg) <-
        c(            
          "False Discovery Rate",
          "Fold Change",
          "Average Read Count",
          "Dispersion",
          "Total Number of miRNAs",
          "Expected Number of DE miRNAs"
        )
      colnames(reg) = paste("Power Graph", 1:length(storageValues$powerCurves))
      return(reg)
    } else{
      return(NULL)
    }
  }, digits = 3, caption = "<b><h4>Parameter List</h4></b>", caption.placement = getOption("xtable.caption.placement", "top"),
  rownames = T)
  #THERE HAS TO BE A BETTER WAY
  output$powerTable1 = renderTable(if(!is.null(storageValues$numPowerCurves) && storageValues$numPowerCurves > 0){storageValues$powerTables[[1]]}else{NULL}, digits = 3, caption = "<b><h4>Power table 1</h4></b>", caption.placement = getOption("xtable.caption.placement", "top"),
                                    caption.width = getOption("xtable.caption.width", NULL))
  output$powerTable2 = renderTable(if(!is.null(storageValues$numPowerCurves) && storageValues$numPowerCurves > 1){storageValues$powerTables[[2]]}else{NULL}, digits = 3, caption = "<b><h4>Power table 2</h4></b>", caption.placement = getOption("xtable.caption.placement", "top"),
                                   caption.width = getOption("xtable.caption.width", NULL))
  output$powerTable3 = renderTable(if(!is.null(storageValues$numPowerCurves) && storageValues$numPowerCurves > 2){storageValues$powerTables[[3]]}else{NULL}, digits = 3, caption = "<b><h4>Power table 3</h4></b>", caption.placement = getOption("xtable.caption.placement", "top"),
                                   caption.width = getOption("xtable.caption.width", NULL))
  output$powerTable4 = renderTable(if(!is.null(storageValues$numPowerCurves) && storageValues$numPowerCurves > 3){storageValues$powerTables[[4]]}else{NULL}, digits = 3, caption = "<b><h4>Power table 4</h4></b>", caption.placement = getOption("xtable.caption.placement", "top"),
                                   caption.width = getOption("xtable.caption.width", NULL))
  output$powerTable5 = renderTable(if(!is.null(storageValues$numPowerCurves) && storageValues$numPowerCurves > 4){storageValues$powerTables[[5]]}else{NULL}, digits = 3, caption = "<b><h4>Power table 5</h4></b>", caption.placement = getOption("xtable.caption.placement", "top"),
                                   caption.width = getOption("xtable.caption.width", NULL))
  
  output$circularplot = renderPlot({
    withProgress(message = "Rendering base", value = 0.1, {
      circos.initializeWithIdeogram()
      circos.clear()
    })
  })
  ################### functions ###################
  
  #function to return list of probes that
  filterProbes <- function(elist, threshVal, percentSamples) {
    if (threshVal == "Global Mean") {
      #get global mean
      eMatrix <- elist$E
      eArea <- nrow(eMatrix) * ncol(eMatrix)
      eAvg <- sum(eMatrix) / eArea
      
      #convert into T/F matrix of greater than global average
      eMatrixTF <- (eMatrix > eAvg)
      
      mirRowSumVector <- rowSums(eMatrixTF)#numeric()
      "i'm pretty sure this is just rowsums
      for (i in 1:nrow(eMatrix)) {
        mirRowSumVector <- c(mirRowSumVector, sum(eMatrixTF[i,]))
      }
      "
      mirKeep <-
        mirRowSumVector > percentSamples * ncol(eMatrix) / 100
      return(mirKeep)
    }
    else{
      eMatrix <- elist$E
      #convert into T/F matrix of greater than global average
      eMatrixTF <- (eMatrix > threshVal)
      
      mirRowSumVector <- rowSums(eMatrixTF)#numeric()
      "
      for (i in 1:nrow(eMatrix)) {
        mirRowSumVector <- c(mirRowSumVector, sum(eMatrixTF[i,]))
      }
      "
      mirKeep <-
        mirRowSumVector > percentSamples * ncol(eMatrix) / 100
      return(mirKeep)
    }
  }
  
  win <- function(filepath){
    is.windows <- FALSE
    
    if(is.windows){
      filepath <- gsub("/", "\\\\",filepath)
    }
    return(filepath)
  }
  
  normalizeMatrix <- function(method, E, housekeepVector, voom) {
    if (method == "Upper Quartile") {
      
      A <- E
      upQ <- 0
      
      for (i in 1:ncol(A)) {
        col <- A[, i]
        denom <- quantile(col[col != 0], names = FALSE)[4]
        A[, i] <- A[, i] / denom
        upQ <- upQ + denom
      }
      upQ <- upQ / ncol(A)
      
      A <- A * upQ
      
    } else if (method == "Total Count") {
      
      A <- E
      tc <- 0 
      
      for (i in 1:ncol(A)) {
        denom <- sum(A[, i])
        A[, i] <- (A[, i] / denom)
        tc <- tc + denom
      }
      tc <- tc / ncol(A)
      
      A <- A * tc
      
    } else if (method == "Median") {
      
      A <- E
      med <- 0
      
      for (i in 1:ncol(A)) {
        col <- A[, i]
        denom <- quantile(col[col != 0], names = FALSE)[3]
        A[, i] <- (A[, i] / denom)
        med <- med + denom
      }
      med <- med / ncol(A)
      
      A <- A * med
      
    } else if (method == "CPM") {
      
      A <- E
      
      for (i in 1:ncol(A)) {
        denom <- sum(A[, i]) / 1000000
        A[, i] <- (A[, i] / denom)
      }
    } else if (method == "Housekeeping Gene") {
      A <- E
      for (i in 1:ncol(A)) {
        denom <- housekeepVector[i]
        A[, i] <- (A[, i] / denom)
      }
      A <- A * 1000000
    }else{
      #if no approved option selected
      return(NULL)
    }
    if(voom){
      A <- voom(A+0.5, plot = FALSE)$E
    }
    else{
      A <- log2(A+1)
    }
    return(A)
  }
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  cleanE = function(E, mod, sv) {
    X = cbind(mod, sv)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(E))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(E - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
  }
  
  plotColor <- function(p, both = FALSE, discrete = TRUE, fill = TRUE, n = NULL, input, rev = FALSE, ...){
    palType <- plotPalType(input)
    pal <- plotPal(input, rev)
    d <- 1
    if(rev){
      d <- -1
    }
    if(discrete){
      if(n < 6){
        n <- 6
      }
      if(fill){
        if(both){
          switch(palType,
                 RColorBrewer = print(p + scale_fill_brewer(palette = pal, direction = d)+ scale_color_brewer(palette = pal, direction = d)),
                 viridis = print(p + scale_fill_viridis(option = pal, discrete = TRUE, direction = d) + scale_color_viridis(option = pal, discrete = TRUE, direction = d)),
                 manual = print(p + scale_fill_manual(values = colorRampPalette(pal)(n))+ scale_color_manual(values = colorRampPalette(pal)(n))),
                 none = print(p))
        } else{
          switch(palType,
                 RColorBrewer = print(p + scale_fill_brewer(palette = pal, direction = d)),
                 viridis = print(p + scale_fill_viridis(option = pal, discrete = TRUE, direction = d)),
                 manual = print(p + scale_fill_manual(values = colorRampPalette(pal)(n))),
                 none = print(p))
        }
      } else{
        switch(palType,
               RColorBrewer = print(p + scale_color_brewer(palette = pal, direction = d)),
               viridis = print(p + scale_color_viridis(option = pal, discrete = TRUE, direction = d)),
               manual = print(p + scale_color_manual(values = colorRampPalette(pal)(n))),
               none  = print(p))
      }
    } else {
      if(fill){
        switch(palType,
               RColorBrewer = print(p + scale_fill_distiller(palette = pal, direction = d)),
               viridis = print(p + scale_fill_viridis(option = pal, direction = d)),
               manual = print(p + scale_fill_gradientn(colors = pal)),
               none = print(p))
      } else{
        switch(palType,
               RColorBrewer = print(p + scale_color_distiller(palette = pal, direction = d)),
               viridis = print(p + scale_color_viridis(option = pal, direction = d)),
               manual = print(p + scale_color_continuous(colors = pal)),
               none = print(p))
      }
    }
  } 
  

  
  plotOrderBoxplot<- function(orderVector = NULL, E, lab, title, cont = FALSE){
    if(is.null(orderVector)){
      orderVector <- 1:ncol(E)
    }
    incProgress(0.333, detail = "Applying sample order")
    
    pOrder <- order(orderVector)
    
    #reorder data matrix
    E <- E[,pOrder]
    
    #add S to sample #'s
    cn <- character()
    for(i in 1:ncol(E)){
      cn <- c(cn, paste0("S", i))
    }
    
    #replace column names 
    colnames(E) <- factor(cn[pOrder])
    
    incProgress(0.333, detail = "Building data frames")
    #form boxplot dataframe
    m <- melt(t(E))
    
    #create group value column and append to df
    m$orderVector <- rep(orderVector[pOrder], nrow(E))
    incProgress(0.333, detail = "Building plot")
    m$X1 <- factor(m$X1, levels = cn[pOrder])
    if(cont){
      p <- ggplot(m, aes(x = X1, y = log2(value))) + geom_boxplot(outlier.shape = 3, outlier.color = '#ff0000', width = 0.6) + ggtitle(toString(title)) + xlab("Sample Number") + ylab ("log2ExpVal") + labs(fill = lab) 
    } else{
      p <- ggplot(m, aes(x = X1, y = log2(value), fill = factor(orderVector))) + geom_boxplot(outlier.shape = 3, outlier.color = '#ff0000', outlier.size = 1, width = 0.6) + ggtitle(toString(title)) + xlab("Sample Number") + ylab ("log2ExpVal") + labs(fill = lab) 
    }
    return(p)
  }
  
  heatmapPal <- function(palType, rev){
    a <- switch(
      palType,
      Default =  '-RdYlBu2:100',
      YellowOrangeRed = colorRampPalette(brewer.pal(9, "YlOrRd"))(255),
      YellowGreenBlue = colorRampPalette(brewer.pal(9, "YlGnBu"))(255),
      RedYellowGreen = colorRampPalette(brewer.pal(11, "RdYlGn"))(255),
      RedYellowBlue = colorRampPalette(brewer.pal(11, "RdYlBu"))(255),
      RedBlue = colorRampPalette(brewer.pal(11, "RdBu"))(255),
      Spectral = colorRampPalette(brewer.pal(11, "Spectral"))(255),
      Rainbow = colorRampPalette(rainbow(6))(255),
      RedGreenBlue = colorRampPalette(c("Red", "Green", "Blue"))(255),
      RedBlackGreen = colorRampPalette(c("#ff0000", "#800000", '#000000', "#008000", "#00ff00"))(255),
      Rainbow2.0 = colorRampPalette(c("#ff0000", "#ffa600", "#ffff00", "#00ff00", "#00ffff", "#0000ff"))(255),
      Viridis = viridis(255),
      Magma = viridis(255, option = "A"),
      Inferno = viridis(255, option = "B", begin = 0),
      Plasma = viridis(255, option = "C", begin = 0)
    )
    if(rev){
      return(rev(a))
    }
    return(a)
  }
  
  plotPal <- function(palType, rev){
    a <- switch(
      palType,
      Default = NULL,
      Greys = "Greys",
      YellowOrangeRed = "YlOrRd",
      YellowGreenBlue = "YlGnBu",
      RedYellowGreen = "RdYlGn",
      RedYellowBlue = "RdYlBu",
      Set1 = "Set1",
      Set2 = "Set2",
      Set3 = "Set3",
      RedBlue = "RdBu",
      Spectral = "Spectral",
      Rainbow = rainbow(6),
      RedGreenBlue = colorRampPalette(c("Red", "Green", "Blue"))(6),
      RedBlackGreen = colorRampPalette(c("#ff0000", "#800000", '#000000', "#008000", "#00ff00"))(6),
      Rainbow2.0 = c("#ff0000", "#ffa600", "#ffff00", "#00ff00", "#00ffff", "#0000ff"),
      Viridis = "D",
      Magma = "A",
      Inferno = "B",
      Plasma = "C"
    )
    if(rev){
      return(rev(a))
    }
    return(a)
  }
  plotPalType <- function(pal){
    a <- switch(
      pal,
      Default = "none",
      Greys = "RColorBrewer",
      YellowOrangeRed = "RColorBrewer",
      YellowGreenBlue = "RColorBrewer",
      RedYellowGreen = "RColorBrewer",
      RedYellowBlue = "RColorBrewer",
      Set1 = "RColorBrewer",
      Set2 = "RColorBrewer",
      Set3 = "RColorBrewer",
      RedBlue = "RColorBrewer",
      Spectral = "RColorBrewer",
      Rainbow = "manual",
      RedGreenBlue = "manual",
      RedBlackGreen = "manual",
      Rainbow2.0 = "manual",
      Viridis = "viridis",
      Magma = "viridis",
      Inferno = "viridis",
      Plasma = "viridis"
    )
    return(a)
  }
  #copy of est_power_curve but improved plotting options
  est_power_curve_better = function (n, w = 1, rho = 2, lambda0 = 5, phi0 = 1, alpha = 0.05, 
                                     f = 0.05, stepsize = 5, squarestep = F, usegradient = T, gradientdetail = 50,
                                     distributionObject, repNumber = 10,...) 
  {
    powerList <- NULL
    sampleSizeList = NULL
    withProgress({
      if(usegradient) { #step based on the gradient, if gradient is tiny, step larger
        i = 1
        sampleSizeList = 1      
        while(sampleSizeList[i] <= n) {
          setProgress(sampleSizeList[i]/n, message = "Generating power values")
          if(is.null(distributionObject) || missing(distributionObject)) {
            if (!missing(f)) {
              powerList[i] <- est_power(n = sampleSizeList[i], 
                                        w = w, rho = rho, lambda0 = lambda0, phi0 = phi0, 
                                        f = f, ...)
            }
            else {
              powerList[i] <- est_power(n = sampleSizeList[i], 
                                        w = w, rho = rho, lambda0 = lambda0, phi0 = phi0, 
                                        alpha = alpha, ...)
            }
          } else {
            powerList[i] = est_power_distribution(n = sampleSizeList[i], f=alpha,rho=rho,distributionObject = distributionObject, repNumber = repNumber)
          }
          if(i > 1) {
            gradient = (powerList[i] - powerList[i-1])/sqrt(sampleSizeList[i] - sampleSizeList[i-1]) #sqrt is to discourage large leaps in sample sizes
            sampleSizeList[i+1] = sampleSizeList[i] + stepsize / (gradientdetail * gradient + 1)
          } else {
            sampleSizeList[i+1] = sampleSizeList[i] + stepsize
          }
          i = i + 1
        } 
        #edge case
        sampleSizeList[i] = n
        if(is.null(distributionObject) || missing(distributionObject)) {
          if (!missing(f)) {
            powerList[i] <- est_power(n = sampleSizeList[i], 
                                      w = w, rho = rho, lambda0 = lambda0, phi0 = phi0, 
                                      f = f, ...)
          }
          else {
            powerList[i] <- est_power(n = sampleSizeList[i], 
                                      w = w, rho = rho, lambda0 = lambda0, phi0 = phi0, 
                                      alpha = alpha, ...)
          }
        } else {
          powerList[i] = est_power_distribution(n = sampleSizeList[i], f=alpha,rho=rho,distributionObject = distributionObject, repNumber = repNumber)
        }
      } else {
        if(squarestep) {
          sampleSizeList = (seq(1, n/stepsize) ^ 2) * stepsize^2/n #todo remove the extra n at the end
        } else {
          if (!n%%stepsize) {      
            sampleSizeList <- c(1, seq(stepsize, n, by = stepsize))
          }
          else {
            sampleSizeList <- c(1, seq(stepsize, n, by = stepsize), n)
          }
        }
        for (i in 1:length(sampleSizeList)) {
          setProgress(i/length(sampleSizeList), message = "Generating power values")
          if(is.null(distributionObject) || missing(distributionObject)) {
            if (!missing(f)) {
              powerList[i] <- est_power(n = sampleSizeList[i], 
                                        w = w, rho = rho, lambda0 = lambda0, phi0 = phi0, 
                                        f = f, ...)
            }
            else {
              powerList[i] <- est_power(n = sampleSizeList[i], 
                                        w = w, rho = rho, lambda0 = lambda0, phi0 = phi0, 
                                        alpha = alpha, ...)
            }
          } else {
            powerList[i] = est_power_distribution(n = sampleSizeList[i], f=alpha,rho=rho,distributionObject = distributionObject, repNumber = repNumber)
          }
        }
      }
    })
    process <- cbind(N = sampleSizeList, Power = powerList)
    if (!missing(f)) {
      parameters <- c(n, f, w, rho, lambda0, phi0)
      names(parameters) <- c("n", "fdr", "w", "rho", "lambda0", 
                             "phi0")
    }
    else {
      parameters <- c(n, alpha, w, rho, lambda0, phi0)
      names(parameters) <- c("n", "alpha", "w", "rho", "lambda0", 
                             "phi0")
    }
    return(list(process = process, parameters = parameters, power = process[nrow(process), 
                                                                            2]))
  }
  #copy of function plot_power_curve but with a legend that makes sense
  plot_power_curve_better = function (result, cexLegend = 1, type = "b", xlab = "Sample Size", 
            ylab = "Power", pch = 16, lwd = 3, las = 1, cex = 1.5, main = "Power Curve", 
            col = "red") 
  {
    if (identical(names(result), c("iter", "f.root", "root", 
                                   "process", "parameters")) || identical(names(result), 
                                                                          c("process", "parameters", "power"))) {
      plot(result$process[, 1], result$process[, 2], type = type, 
           xlab = xlab, ylab = ylab, pch = pch, lwd = lwd, las = las, 
           cex = cex, main = main, col = col)
      legend("bottomright", legend = paste(paste(names(result$parameters), 
                                                 result$parameters, sep = "="), collapse = ";"), bty = "n", 
             text.col = "red", cex = 0.9)
      abline(h = result$parameters["power"], lty = 2, col = "grey")
    }
    else {
      if (length(result) > 5) {
        result <- result[(length(result) - 4):length(result)]
        warning("At most 5 curves were allowed in plot_power_curve function, the last 5 in result parameter will be used")
      }
      resultRange <- apply(sapply(result, function(x) (x$process)[nrow(x$process), 
                                                                  ]), 1, max)
      plot(c(0, resultRange[1]), c(0, resultRange[2]), type = "n", 
           xlab = xlab, ylab = ylab, las = las, main = main)
      col <- c("brown1", "steelblue2", "mediumpurple2", "seagreen3", 
               "lightgoldenrod")
      legendEach <- ""
      for (x in 1:length(result)) {
        resultEach <- result[[x]]
        lines(resultEach$process[, 1], resultEach$process[, 
                                                          2], type = "b", pch = 16, lwd = 3, cex = 1.5, 
              col = col[x])
      }
      legendEach = paste("Power Graph",1:length(result))
      abline(h = result[[length(result)]]$parameters["power"], 
             lty = 2, col = "grey")
      legend("bottomright", legend = legendEach, bty = "n", 
             text.col = col[1:length(result)], cex = cexLegend)
    }
  }
})


