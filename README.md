# About miRShiny
miRShiny is an R Shiny-based interactive pipeline for the processing and visualization of microRNA sequencing data.

# Tabs/Features
- Data Upload
- Data Pre-Process
- Quality Control
- Differential Analysis
- Genomic Visualization
- Individual Feature Visualization
- Data Download
- Accuracy Evaluation
- Statistical Power Analysis

# Input Data Format
Data input requires two files: **Expression Data** and a **Conditions File**.

**Expression Data:** A counts matrix or similar numerical matrix, with rows corresponding to features, and columns corresponding to samples.

*Requirements:*
1. Uploaded in tab-seperated, comma-seperated, or similar text format
1. Formatted as a matrix with dimensions [i,j], with multiple rows *i* and columns *j*
1. One row and column of sample and feature name annotation, with NO repeated names

**Conditions File:** A columnar text file including sample conditions and other sample annotations.

*Requirements:*
1. Uploaded in tab-seperated, comma-seperated, or similar text format
1. Number of sample information rows (discounting header) in each column is equal to *j*, the number of columns/samples in the expression data
1. Includes at minimum one complete column titled `condition` to identify the state of each sample
1. Exactly one header row

*Optional columns:* Additional column information is accepted in the Conditions File for various purposes.
- `group`: for subsetting the uploaded data set
- `batch`: for considering results between batches and correcting batch error
- `normalizer`: for a custom scaling vector in normalization

# Screenshots
<img src="https://user-images.githubusercontent.com/5732925/41006669-6e5f03f0-68d7-11e8-8155-6171f6627210.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/5732925/41006676-748af798-68d7-11e8-80e7-1fcaba95b72c.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/5732925/41006678-760bf2f2-68d7-11e8-90a2-59282df11f41.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/5732925/41006679-76f624da-68d7-11e8-94e1-98dbda3199da.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/5732925/41006682-77df9372-68d7-11e8-9210-45d680308c8b.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/5732925/41006684-78b68756-68d7-11e8-8d0f-96721aabeef5.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/5732925/41006685-79819a18-68d7-11e8-906f-e2a38b968024.png" width="45%"></img> <img src="https://user-images.githubusercontent.com/5732925/41006687-7a3c12d0-68d7-11e8-8015-156db6219b31.png" width="45%"></img> 

# Developed With
- R
  - R Core Team (2017). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/
- Bioconductor
  - Orchestrating high-throughput genomic analysis with Bioconductor. W. Huber, V.J. Carey, R. Gentleman,..., M. Morgan Nature Methods, 2015:12, 115. https://www.bioconductor.org/
- shiny
  - Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017). shiny: Web Application Framework for R. R package version 1.0.3. https://CRAN.R-project.org/package=shiny
- Limma
  - Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
- edgeR
  - McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research 40, 4288-4297
- DESeq2
  - Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
- ggplot2
  - H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
- sva
  - Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Storey JD, Zhang Y and Torres LC (2017). sva: Surrogate Variable Analysis. R package version 3.24.4.
- NMF
  - Renaud Gaujoux, Cathal Seoighe (2010). A flexible R package for nonnegative matrix factorization. BMC Bioinformatics 2010, 11:367.
- reader
  - Nicholas Cooper (2017). reader: Suite of Functions to Flexibly Read Data from Files. R package version 1.0.6. https://CRAN.R-project.org/package=reader
- viridis
  - Simon Garnier (2017). viridis: Default Color Maps from 'matplotlib'. R package version 0.4.0. https://CRAN.R-project.org/package=viridis
- RnaSeqSampleSize
  - Zhao S, Li C, Guo Y, Sheng Q and Shyr Y (2017). RnaSeqSampleSize: RnaSeqSampleSize. R package version 1.10.0. http://bioconductor.org/packages/release/bioc/html/RnaSeqSampleSize.html
- circlize
  - Gu, Z. (2014) circlize implements and enhances circular visualization in R. Bioinformatics. 10.1093/bioinformatics/btu393
- ComplexHeatmap
  - Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics. https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
- openxlsx
  - Alexander Walker (2017). openxlsx: Read, Write and Edit XLSX Files. R package version 4.0.17. https://CRAN.R-project.org/package=openxlsx
- heatmaply
  - Tal Galili, Alan O'Callaghan, Jonathan Sidi, Carson Sievert; heatmaply: an R package for creating interactive cluster heatmaps for online publishing, Bioinformatics, , btx657. https://doi.org/10.1093/bioinformatics/btx657
