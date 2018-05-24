require(shiny)
require(ggplot2)
require(limma)
require(RColorBrewer)
require(DESeq2)
require(edgeR)
require(reshape)
require(viridis)
require(sva)
require(grid)
require(reader)
require(NMF)
require(plyr)
require(shinydashboard)
require(plotly)

fontSize = "font-size:150%"
headerFontSize = "font-size:200%"
subFontSize = "font-size:130%"

dashboardPage(skin="yellow",
              dashboardHeader(
                title = "miRShiny",
                titleWidth = 250),
              dashboardSidebar(
                width = 250,
                sidebarMenu(
                  menuItem("About", icon = icon("info-circle"), startExpanded = TRUE,
                           menuSubItem("Input Data Format", tabName = "dataFormat"),
                           menuSubItem("Developed With", tabName = "devWith"),
                           menuSubItem("Change Log", tabName = "chLog")
                  ),
                  menuItem("Data Analysis", icon = icon("bar-chart"), startExpanded = FALSE,
                           menuSubItem("Data Upload", tabName = "Upload", icon = icon("upload")),
                           menuSubItem("Data Pre-Process", tabName = "Pre-Process", icon = icon("gears")),
                           menuSubItem("Quality Control", tabName = "Quality", icon = icon("certificate")),
                           menuSubItem("Differential Analysis", tabName = "Differential", icon = icon("calculator")),
                           menuSubItem("Genomic Visualization", tabName = "Genome", icon = icon("circle-o")),
                           menuSubItem("Individual Feature Visualization", tabName = "Visualization", icon = icon("eye")),
                           menuSubItem("Data Download", tabName = "Download", icon = icon("download")),
                           menuSubItem("Accuracy Evaluation", tabName = "Accuracy", icon = icon("bullseye")),
                           menuSubItem("Statistical Power Analysis", tabName = "Power", icon = icon("power-off"))
                  )
                )
              ),
              dashboardBody(
                tags$head(tags$style(HTML('
                              .main-header .logo {
                              font-family: Geneva, "Verdana", sans-serif;
                              font-weight: bold;
                              font-size: 24px;
                              }
                              .largeFont{
                                font-size: 50px;
                              }
                              '))),
                tabItems(
                  
                  tabItem(tabName = "dataFormat",
                          fluidRow(
                            column(
                              5,
                              offset = 0,
                              HTML(
                                "<h2 style='color:#4d4d4d; font-size: 250%; font-family:helvetica'><strong>About miRShiny</strong></h2>"
                              ),
                              HTML(
                                "<p style='font-size: 150%; font-family:helvetica'>
                                  miRShiny is an R Shiny-based interactive pipeline for the processing and visualization of microRNA sequencing data. </p>"
                              ),
                              HTML(
                                "<p style='font-size: 150%; font-family:helvetica'>
                                    <strong>Last updated:</strong> 1/17/2018
                                    </p>"
                              )
                            ),
                            column(
                              3,
                              align = "center",
                              offset = 3,
                              img(src = 'isbLogo.png', height = 150)
                            )
                          ),
                          HTML("<br>"),
                          wellPanel(
                            HTML(
                              "<p style='font-size: 200%;font-family:helvetica'>
                              <strong><u>Input Data Format</u></strong>
                              </p>"
                            ),
                            HTML("<br>"),
                            HTML(
                              "<p style='font-size: 150%;font-family:helvetica'>
                                          <strong><em>Warning:</em></strong>
                                          Avoid special characters like <code>whitespace</code>, <code>,</code>, <code>#</code> in sample/feature/condition labels to prevent parsing errors.
                                          </p>"
                            ),
                            HTML("<br>"),
                            HTML(
                              "<p style='font-size: 150%;font-family:helvetica'>
                                          <strong>Expression Data:</strong>
                                          A counts matrix or similar numerical matrix, with rows corresponding to features, and columns corresponding to samples.
                                          </p>"
                            ),
                            HTML(
                              "<p style='font-size: 150%;font-family:helvetica'>
                                          <em>Requirements:</em>
                                          <ol style='font-size: 150%;font-family:helvetica'>
                                          <li>Uploaded in tab-seperated, comma-seperated, or similar text format</li>
                                          <li>Formatted as a matrix with dimensions [i,j], with multiple rows (i) and columns (j)</li>
                                          <li>One row and column of sample and feature name annotation, with NO repeated names</li>
                                          </ol>
                                          </p>"
                            ),
                            HTML("<br>"),
                            HTML(
                              "<p style='font-size: 150%;font-family:helvetica'>
                                          <strong>Conditions File:</strong>
                                          A columnar text file including sample conditions and other sample annotations.
                                          </p>"
                            ),
                            HTML(
                              "<p style='font-size: 150%;font-family:helvetica'>
                                          <em>Requirements:</em>
                                          <ol style='font-size: 150%;font-family:helvetica'>
                                          <li >Uploaded in tab-seperated, comma-seperated, or similar text format</li>
                                          <li>Number of sample information rows (discounting header) in each column is equal to <em>j</em>, the number of columns/samples in the expression data.</li>
                                          <li>Includes at minimum one complete column titled <code>condition</code> to identify the state of each sample</li>
                                          <li>Exactly one header row</li>
                                          </ol>
                                          </p>"
                            ),HTML(
                              "<p style='font-size: 150%;font-family:helvetica'>
                                          <em>Optional columns:</em> Additional column information is accepted in the Conditions File for various purposes.
                                          <ul style='font-size: 150%;font-family:helvetica'>
                                          <li><code>group</code>: for subsetting the uploaded data set</li>
                                          <li><code>batch</code>: for considering results between batches and correcting batch error</li>
                                          <li><code>normalizer</code>: for a custom scaling vector in normalization</li>
                                          </ul>
                                          </p>"
                            )
                          )
                  ),
                  tabItem(tabName = "devWith",
                          fluidRow(
                            column(
                              5,
                              offset = 0,
                              HTML(
                                "<h2 style='color:#4d4d4d; font-size: 250%; font-family:helvetica'><strong>About miRShiny</strong></h2>"
                              ),
                              HTML(
                                "<p style='font-size: 150%; font-family:helvetica'>
                                miRShiny is an R Shiny-based interactive pipeline for the processing and visualization of microRNA sequencing data. </p>"
                              ),
                              HTML(
                                "<p style='font-size: 150%; font-family:helvetica'>
                                <strong>Last updated:</strong> 1/17/2018
                                </p>"
                              )
                            ),
                            column(
                              3,
                              align = "center",
                              offset = 3,
                              img(src = 'isbLogo.png', height = 150)
                            )
                          ),
                          HTML("<br>"),
                          wellPanel(
                            HTML(
                              "<p style='font-size: 200%;font-family:helvetica'>
                              <strong><u>Developed With</u></strong>
                              </p>"
                            ),
                            HTML(
                              "<p style='font-size: 150%;font-family:helvetica'>
                              <ul style='font-size: 130%;font-family:helvetica'>
                              <li><em>R</em>:
                              <ul>
                              <li>R Core Team (2017). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.</li>
                              <li><a href = 'https://www.R-project.org/'>https://www.R-project.org/</a></li>
                              </ul>
                              </li>
                              <li><em>Bioconductor</em>:
                              <ul>
                              <li>Orchestrating high-throughput genomic analysis with Bioconductor. W. Huber, V.J. Carey, R. Gentleman,..., M. Morgan Nature Methods, 2015:12, 115.</li>
                              <li><a href = 'https://www.bioconductor.org/'>https://www.bioconductor.org/</a></li>
                              </ul>
                              </li>
                              <li><em>shiny</em>:
                              <ul>
                              <li>Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017). shiny: Web Application Framework for R. R package version 1.0.3.</li>
                              <li><a href = 'https://CRAN.R-project.org/package=shiny'>https://CRAN.R-project.org/package=shiny</a></li>
                              </ul>
                              </li>
                              <li><em>Limma</em>:
                              <ul>
                              <li>Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.</li>
                              </ul>
                              </li>
                              <li><em>edgeR</em>:
                              <ul>
                              <li>McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research 40, 4288-4297</li>
                              </ul>
                              </li>
                              <li><em>DESeq2</em>:
                              <ul>
                              <li>Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)</li>
                              </ul>
                              </li>
                              <li><em>ggplot2</em>:
                              <ul>
                              <li>H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.</li>
                              </ul>
                              </li>
                              <li><em>sva</em>:
                              <ul>
                              <li>Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Storey JD, Zhang Y and Torres LC (2017). sva: Surrogate Variable Analysis. R package version 3.24.4.</li>
                              </ul>
                              </li>
                              <li><em>NMF</em>:
                              <ul>
                              <li>Renaud Gaujoux, Cathal Seoighe (2010). A flexible R package for nonnegative matrix factorization. BMC Bioinformatics 2010, 11:367.</li>
                              </ul>
                              </li>
                              <li><em>reader</em>:
                              <ul>
                              <li>Nicholas Cooper (2017). reader: Suite of Functions to Flexibly Read Data from Files. R package version 1.0.6.</li>
                              <li><a href = 'https://CRAN.R-project.org/package=reader'>https://CRAN.R-project.org/package=reader</a></li>
                              </ul>
                              </li>
                              <li><em>viridis</em>:
                              <ul>
                              <li>Simon Garnier (2017). viridis: Default Color Maps from 'matplotlib'. R package version 0.4.0.</li>
                              <li><a href = 'https://CRAN.R-project.org/package=viridis'>https://CRAN.R-project.org/package=viridis</a></li>
                              </ul>
                              </li>
                              <li><em>RnaSeqSampleSize</em>:
                              <ul>
                              <li>Zhao S, Li C, Guo Y, Sheng Q and Shyr Y (2017). RnaSeqSampleSize: RnaSeqSampleSize. R package version 1.10.0.</li>
                              <li><a href = 'http://bioconductor.org/packages/release/bioc/html/RnaSeqSampleSize.html'>http://bioconductor.org/packages/release/bioc/html/RnaSeqSampleSize.html</a></li>
                              </ul>
                              </li>
                              <li><em>circlize</em>:
                              <ul>
                              <li>Gu, Z. (2014) circlize implements and enhances circular visualization in R. Bioinformatics.</li>
                              <li><a href = '10.1093/bioinformatics/btu393'>10.1093/bioinformatics/btu393</a></li>
                              </ul>
                              </li>
                              </ul>
                              </p>"
                            )
                          )
                  ),
                  tabItem(tabName = "chLog",
                          fluidRow(
                            column(
                              5,
                              offset = 0,
                              HTML(
                                "<h2 style='color:#4d4d4d; font-size: 250%; font-family:helvetica'><strong>About miRShiny</strong></h2>"
                              ),
                              HTML(
                                "<p style='font-size: 150%; font-family:helvetica'>
                                miRShiny is an R Shiny-based interactive pipeline for the processing and visualization of microRNA sequencing data. </p>"
                              ),
                              HTML(
                                "<p style='font-size: 150%; font-family:helvetica'>
                                <strong>Last updated:</strong> 1/17/2018
                                </p>"
                              )
                            ),
                            column(
                              3,
                              align = "center",
                              offset = 3,
                              img(src = 'isbLogo.png', height = 150)
                            )
                          ),
                          HTML("<br>"),
                          wellPanel(
                            HTML(
                              "<p style='font-size: 200%;font-family:helvetica'>
                              <strong><u>Change Log</u></strong>
                              </p>"
                            ),
                            HTML(
                              "<p style='font-family:helvetica'>
                              <ul style='font-size: 150%;font-family:helvetica'>
                              <li><em>11/8/2017</em></li>
                              <ul>
                              <li><em>Feature Changes</em></li>
                              <ul>
                              <li>Added a new tab: <b>Statistical Power Analysis.</b></li>
                              <ul>
                              <li>Given parameters of your data, it will determine the relation between sample size and statistical power.</li>
                              <li>Can automatically set parameters based on your uploaded data.</li>
                              <li>Includes a variety of different methods to plot this relationship.</li>
                              </ul>
                              </ul>
                              </ul>
                              <li><em>8/2/2017</em></li>
                              <ul>
                              <li><em>Usage Changes</em></li>
                              <ul>
                              <li>Changed required label in sample information file to <code>condition</code> from <code>conditions</code> for compatibility purposes.</li>
                              </ul>
                              <li><em>Feature Changes</em></li>
                              <ul>
                              <li>PCA Plots/Correlation Coefficient Matrix now generatable in Quality Control tab.</li>
                              <li>Quality Control plots can now be sorted/grouped/facetted by various factors</li>
                              <li>Batch correction from package 'sva' now enabled for experiments with batch effects.</li>
                              <li>Differentially expressed microRNAs can now be plotted side-by-side for comparison.</li>
                              <li>Expanded information shown in status boxes.</li>
                              </ul>
                              <li><em>Interface Changes</em></li>
                              <ul>
                              <li>Data download links moved to a new tab.</li>
                              <li>Modals are now generated after successful upload and processing of data.</li>
                              <li>Loading bars reimplemented for most processes.</li>
                              <li>Minor label text changes.</li>
                              <li>Plot colors darkened for increased clarity.</li>
                              </ul>
                              <li><em>Bug Fixes</em></li>
                              <ul>
                              <li>Filtering threshold now displays properly when subsetting data.</li>
                              </ul>
                              </ul>
                              </ul>
                              </p>"
                            )
                          )
                  ),
                  tabItem(tabName = "Upload",
                          fluidRow(
                            br(),
                            column(
                              5,
                              wellPanel(
                                checkboxInput(
                                  inputId = "sampleData",
                                  label = "Use Sample Data?",
                                  value = FALSE
                                ),
                                conditionalPanel(
                                  condition = "input.sampleData",
                                  HTML("<hr>"),
                                  HTML(
                                    "<p style='font-size: 150%;font-family:helvetica'>
                                    Sample Data Information
                                    <ul style='font-size: 150%;font-family:helvetica'>
                                    <li><strong>GEO Accession: </strong>GSE71008</li>
                                    <li><strong>Dataset Name: </strong>Plasma extracellular RNA profiles in Healthy and Cancer patients</li>
                                    <li><strong>Source: </strong>Medical College of Wisconsin, Department of Pathology</li>
                                    <li><strong>Platform: </strong>Illumina Genome Analyzer (HS), NEB sRNA library preparation kit</li>
                                    <li><strong>Dataset Description: </strong>Small RNA Sequencing of 50 healthy controls, 100 colorectal cancer patients (25 each of Stage I-IV), 36 prostate cancer patients, and 3 pancreatic cancer patients.</li>
                                    <li><strong>Citation: </strong>Yuan T, Huang X, Woodcock M, Du M et al. Plasma extracellular RNA profiles in healthy and cancer patients. Sci Rep 2016 Jan 20;6:19413.</li>
                                    </ul>
                                    </p>"
                                  )
                                ),
                                conditionalPanel(
                                  condition = "!input.sampleData",
                                  HTML(
                                    "<p style='font-size: 130%; font-family:helvetica'>
                                    <strong><em><font color = '#bd4147'>Warning:</font></em></strong>
                                    Avoid special characters like <code>whitespace</code>, <code>,</code>, <code>#</code> in sample/feature/condition labels to prevent parsing errors.
                                    </p>"
                                    
                                  ),
                                  h4("Expression Data", style = fontSize),
                                  helpText("Counts matrix or similar matrix-form data including whole number read counts. Include exactly one row/column of sample/feature name annotation. Most common delimited text formats accepted.", style = fontSize),
                                  fileInput(
                                    inputId = "files",
                                    label = "Upload expression data:",
                                    multiple = FALSE,
                                    accept = c('.txt', '.csv', '.tsv', '.sum')
                                  ),
                                  HTML('<hr style="border-color: #bfbfbf">'),
                                  h4("Sample Information", style = fontSize),
                                  helpText("File that includes columns of sample annotation. Number of rows must match number of samples. Must include column named \"condition\". Can include optional columns \"group\",\"batch\", and \"normalizer\" (see About page for further information). Most common delimited text formats accepted.", style = fontSize),
                                  fileInput(
                                    inputId = "condition",
                                    label = "Upload conditions file:",
                                    multiple = FALSE,
                                    accept = c('.txt', '.csv', '.tsv', '.sum')
                                  )
                                )
                              ),
                              wellPanel(
                                actionButton(
                                  inputId = "upload",
                                  label = "Read Files",
                                  width = '100%'
                                )
                              )
                            ),
                            column(
                              7,
                              wellPanel(
                                align = "center",
                                h4("File Upload Status", style = fontSize),
                                div(tableOutput("uploadStatusTable"), #table output here
                                    style = "font-size:130%; overflow:scroll")
                              )
                            )
                          )
                  ),
                  tabItem(tabName = "Pre-Process",
                          fluidRow(
                            br(),
                            column(
                              4,
                              conditionalPanel(
                                condition = "output.groupExists",
                                wellPanel(
                                  h4("Subset Data Files"),
                                  helpText(
                                    "Select groups (from \"group\" column in condition file) to include them in data processing. If no groups are selected, all groups will be included."
                                  ),
                                  htmlOutput("subsetUI")
                                )
                              ),
                              wellPanel(
                                h4("Normalization"),
                                selectInput(
                                  inputId = "norm",
                                  label = "Normalization Method",
                                  choices = c(
                                    "Total Count",
                                    "Median",
                                    "Quantile",
                                    "Scale",
                                    "CyclicLoess",
                                    "Upper Quartile",
                                    "DESeq",
                                    "TMM",
                                    "VSN",
                                    "CPM",
                                    "Housekeeping Gene"
                                  ),
                                  selected = "CPM"
                                ),
                                
                                conditionalPanel(
                                  condition = "input.norm == 'Housekeeping Gene'",
                                  htmlOutput("housekeepUI")
                                ),
                                conditionalPanel(
                                  helpText("Voom: uses mean-variance generated precision weights in normalization."),
                                  condition = "input.norm != 'VSN'",
                                  checkboxInput(
                                    inputId = "voom",
                                    label = "Use voom?",
                                    value = FALSE
                                  )
                                )
                              ),
                              conditionalPanel(
                                condition = "output.batchExists",
                                wellPanel(
                                  h4("Correct for batch error?"),
                                  helpText("Use with quality control figures. Methods are often unreliable."),
                                  selectInput(
                                    inputId = "batchCorrect",
                                    label = "Batch Correction Method",
                                    choices = c("None", "comBat", "sva", "Targetted Recentering"),
                                    selected = "None"
                                  ),
                                  conditionalPanel(
                                    condition = "input.batchCorrect == 'sva'",
                                    selectInput(
                                      inputId = "svaN",
                                      label = "Number of surrogate variables to use?",
                                      choices = c("Auto", 1:20),
                                      selected = "Auto"
                                    )
                                  ),
                                  conditionalPanel(
                                    condition = "input.batchCorrect == 'Targetted Recentering'",
                                    helpText("Targetted groups will be scaled to mean value of untargetted groups. An empty or complete selection will result in no correction."),
                                    htmlOutput("trGrpsUI")
                                  )
                                )
                              ),
                              wellPanel(
                                h4("Data Filtering"),
                                helpText("Removes a feature unless expressed above set threshold on designated quantity of samples."),
                                selectInput(
                                  inputId = "filtThresh",
                                  label = "Filtering Threshold Type",
                                  choices = c("Global Mean", "Value"),
                                  selected = "Global Mean"
                                ),
                                conditionalPanel(
                                  condition = "input.filtThresh == 'Value'",
                                  numericInput(
                                    inputId = "filtThreshVal",
                                    label = "Expression Threshold Value",
                                    min = 0,
                                    max = 1000,
                                    value = 5,
                                    step = 1
                                  )
                                ),
                                HTML('<hr style="border-color: #bfbfbf;">'),
                                helpText("Required percentage/number of samples in which a feature must be expressed above the threshold value."),
                                htmlOutput("filtPercent"),
                                strong("Sample # Threshold:"),
                                verbatimTextOutput("filtSampleNumber")
                              ),
                              wellPanel(
                                actionButton(
                                  inputId = "processData",
                                  label = "Process Data",
                                  width = '100%'
                                )
                              )
                            ),
                            column(8,
                                   wellPanel(align = "center", h4("Pre-Processing Status"), div(tableOutput("normStatusTable"), style = "font-size:100%; overflow: scroll")),
                                   wellPanel(
                                     h4(align = "center", "Pre-Processing Details"),
                                     HTML(
                                       "<p style='font-family:helvetica'>
                                       <strong>Normalization Options:</strong>
                                       <ul>
                                       <li><em>Total Count: </em>Linear scaling based on each sample's total counts.</li>
                                       <li><em>Median: </em>Linear scaling based on each sample's median counts.</li>
                                       <li><em>Quantile: </em>Quantile normalization method from 'limma' package.</li>
                                       <li><em>Scale: </em>Scale normalization method from 'limma' package.</li>
                                       <li><em>Cyclic Loess: </em>Loess normalization method from 'Limma' package.</li>
                                       <li><em>Upper Quartile: </em>Linear scaling based on the samples' upper quartile values.</li>
                                       <li><em>DESeq: </em>Normalization method from the 'DESeq2' package.</li>
                                       <li><em>TMM: </em>Trimmed Means of M-Values, implemented in the 'limma' package.</li>
                                       <li><em>VSN: </em>Variance Stabilization Normalization, implemented in the 'vsn' package.</li>
                                       <li><em>CPM: </em>Linear scaling based on library counts per million.</li>
                                       <li><em>Housekeeping Gene: </em>Linear scaling using highly expressed features.</li>
                                       </ul>
                                       </p>"
                                     ),
                                     HTML(
                                       "<p style='font-family:helvetica'>
                                       <strong>Pre-Processing Steps:</strong>
                                       <ol>
                                       <li>Uploaded data is cleaned to replace non-numeric values.</li>
                                       <li>Data is subset by groups, if enabled.</li>
                                       <li>Low count features are marked by user criteria to be filtered out.</li>
                                       <li>Data is normalized.</li>
                                       <li><code>voom</code> function is applied, if enabled.</li>
                                       <li>Batch correction is applied, if enabled.</li>
                                       <li>Data is log-transformed.</li>
                                       <li>Marked low-value features (step 3) are removed.</li>
                                       </ol>
                                       </p>"
                                     )
                                   )
                            )
                          )
                  ),
                  tabItem(tabName = "Quality",
                          fluidRow(
                            br(),
                            column(
                              3,
                              wellPanel(
                                selectInput(
                                  inputId = "plotType",
                                  label = "Select plot to generate:",
                                  choices = c("Boxplot", "PCA Plot", "Correlation Coefficient Matrix")
                                ),
                                conditionalPanel(
                                  condition = "input.plotType == 'Boxplot'",
                                  htmlOutput("sortBoxplotUI")
                                ),
                                conditionalPanel(
                                  condition = "input.plotType == 'PCA Plot'",
                                  htmlOutput("PCAtypeUI"),
                                  htmlOutput("sortPCAUI")
                                ),
                                conditionalPanel(
                                  condition = "input.plotType == 'Correlation Coefficient Matrix'",
                                  selectInput(
                                    inputId = "CCcalcType",
                                    label = "Calculation Method",
                                    choices = c("Pearson", "Spearman", "Kendall"),
                                    selected = "Pearson"
                                  ),
                                  htmlOutput("CCtypeUI")
                                ),
                                actionButton(
                                  inputId = "plotQC",
                                  label = "Plot!",
                                  width = '100%'
                                )
                              ),
                              wellPanel(
                                checkboxInput(
                                  inputId = "QCcolorOptions",
                                  label = "Show color options?",
                                  value = FALSE
                                ),
                                conditionalPanel(
                                  condition = "input.QCcolorOptions",
                                  selectInput(
                                    inputId = "plotColors",
                                    label = "Color Palette",
                                    choices =
                                      c(
                                        'Default',
                                        'YellowOrangeRed',
                                        'YellowGreenBlue',
                                        'RedYellowGreen',
                                        'RedYellowBlue',
                                        'RedBlue',
                                        'Set1',
                                        'Set2',
                                        'Set3',
                                        'Spectral',
                                        'Rainbow',
                                        'RedGreenBlue',
                                        'RedBlackGreen',
                                        'Rainbow2.0',
                                        'Viridis',
                                        'Magma',
                                        'Inferno',
                                        'Plasma'
                                      ),
                                    selected = 'Default'
                                  ),
                                  checkboxInput(
                                    inputId = "colorsRev",
                                    label = "Reverse Palette?",
                                    value = FALSE
                                  )
                                )
                              )
                            ),
                            column(
                              9,
                              align = "center",
                              htmlOutput("upperPlotUI"),
                              #plotlyOutput("upperPlotUI"),
                              HTML('<br>'),
                              htmlOutput("lowerPlotUI")
                            )
                          )
                  ),
                  tabItem(tabName = "Differential",
                          fluidRow(
                            br(),
                            column(
                              3,
                              wellPanel(
                                h4("Select comparison groups"),
                                helpText("Condition group overlaps will not be processed. Invalid condition names may have been slightly modified"),
                                htmlOutput("condition1UI"),
                                htmlOutput("condition2UI")
                              ),
                              wellPanel(
                                h4("MA Plot Threshold"),
                                numericInput(
                                  inputId = "AveEcut",
                                  label = "A-value Cutoff",
                                  value = 7.5,
                                  min = 0,
                                  max = 10000,
                                  step = 1
                                )
                              ),
                              wellPanel(
                                h4("Volcano Plot Thresholds"),
                                numericInput(
                                  inputId = "logFCcut",
                                  label = "Fold Change Cutoff",
                                  value = 1.2,
                                  min = 1.0,
                                  max = 10,
                                  step = 0.05
                                ),
                                numericInput(
                                  inputId = "pValCut",
                                  label = "p-Value Cutoff",
                                  value = 0.05,
                                  min = 0,
                                  max = 0.1,
                                  step = 0.001
                                )
                              ),
                              wellPanel(
                                h4("Heatmap"),
                                checkboxInput(
                                  inputId = "dendClustR",
                                  label = "Cluster features?",
                                  value = F
                                ),
                                checkboxInput(
                                  inputId = "dendClustC",
                                  label = "Cluster samples?",
                                  value = F
                                ),
                                selectInput(
                                  inputId = "hmColors",
                                  label = "Heatmap Color Palette",
                                  choices =
                                    c(
                                      'Default',
                                      'YellowOrangeRed',
                                      'YellowGreenBlue',
                                      'RedYellowGreen',
                                      'RedYellowBlue',
                                      'RedBlue',
                                      'Spectral',
                                      'Rainbow',
                                      'RedGreenBlue',
                                      'RedBlackGreen',
                                      'Rainbow2.0',
                                      'Viridis',
                                      'Magma',
                                      'Inferno',
                                      'Plasma'
                                    ),
                                  selected = "RedBlackGreen"
                                ),
                                checkboxInput(
                                  inputId = "hmRev",
                                  label = "Reverse Color Palette?",
                                  value = FALSE
                                )
                              ),
                              wellPanel(
                                h4("Correlation Matrix Options"),
                                selectInput(
                                  inputId = "CCcalcType",
                                  label = "Calculation Method",
                                  choices = c("Pearson", "Spearman", "Kendall"),
                                )
                              ),
                              wellPanel(
                                checkboxInput(
                                  inputId = "CCcolorOptions",
                                  label = "Show color options?",
                                  value = FALSE
                                ),
                                conditionalPanel(
                                  condition = "input.CCcolorOptions",
                                  selectInput(
                                    inputId = "CCplotColors",
                                    label = "Color Palette",
                                    choices =
                                      c(
                                        'Default',
                                        'YellowOrangeRed',
                                        'YellowGreenBlue',
                                        'RedYellowGreen',
                                        'RedYellowBlue',
                                        'RedBlue',
                                        'Set1',
                                        'Set2',
                                        'Set3',
                                        'Spectral',
                                        'Rainbow',
                                        'RedGreenBlue',
                                        'RedBlackGreen',
                                        'Rainbow2.0',
                                        'Viridis',
                                        'Magma',
                                        'Inferno',
                                        'Plasma'
                                      ),
                                    selected = 'Default'
                                  ),
                                  checkboxInput(
                                    inputId = "CCcolorsRev",
                                    label = "Reverse Palette?",
                                    value = FALSE
                                  )
                                )
                              )
                            ),
                            column(
                              9,
                              align = "center",
                              plotlyOutput("MAPlot", height = 400, width = '85%'),
                              HTML('<br>'),
                              HTML('<br>'),
                              plotlyOutput("volcanoPlot", height = 400, width = '85%'),
                              HTML('<br>'),
                              HTML('<br>'),
                              div(tableOutput("sigMirTable"), style = "font-size:100%"),
                              HTML('<br>'),
                              HTML('<br>'),
                              #htmlOutput("heatmapUI"),
                              #plotlyOutput(renderPlotly({heatmaply(mtcars)}))
                              plotlyOutput("heatmapUI", height = 800, width = '85%'),
                              HTML('<br>'),
                              HTML('<br>'),
                              htmlOutput("corrmapUI")
                            )
                          )
                  ),
                  tabItem(tabName = "Genome",
                        fluidRow(
                          column(3,
                            wellPanel(
                              h4("Plot Parameters"),
                              selectInput(
                                inputId = "genomicPlotType",
                                label = "Plot Type",
                                choices = c("Circular", "Barchart"),
                                selected = "Circular"
                              ),
                              numericInput(
                                inputId = "circlePlotBandWidthMult",
                                label = "Band width multiplier",
                                value = 1,
                                min = 0.25,
                                max = 5,
                                step = 0.25
                              ),
                              numericInput(
                                inputId = "circlePlotBandColorExp",
                                label = "Band color exponent",
                                value = 2,
                                min = 1,
                                max = 4,
                                step = 1
                              ),
                              conditionalPanel(
                                condition = "input.genomicPlotType == 'Circular'",
                                numericInput(
                                  inputId = "circlePlotLinkWidthMult",
                                  label = "Link width multiplier",
                                  value = 1,
                                  min = 0.25,
                                  max = 5,
                                  step = 0.25
                                ),
                                numericInput(
                                  inputId = "circlePlotLinkColorExp",
                                  label = "Link color exponent",
                                  value = 4,
                                  min = 1,
                                  max = 8,
                                  step = 1
                                )
                              )
                            ),
                            wellPanel(
                              actionButton(
                                inputId = "circlePlotButton",
                                label = "Update Plot",
                                width = '100%'
                              )
                            )
                          ),
                          column(9, 
                                 plotOutput("circularplot", width = "100%", height = "1000px")
                          )
                        )
                  ),
                  tabItem(tabName = "Visualization",
                          fluidRow(
                            br(),
                            column(
                              3,
                              wellPanel(
                                htmlOutput('deSelectUI'),
                                HTML('<hr>'),
                                selectInput(
                                  inputId = 'dePlotType',
                                  label = "Plot Type",
                                  choices = c("Boxplot", "Violin Plot", "Split Violin Plot"),
                                  selected = "Boxplot"
                                ),
                                conditionalPanel(
                                  condition = "input.dePlotType == 'Boxplot'",
                                  checkboxInput(inputId = "deNotch", label = "Notch Boxes?", value = FALSE),
                                  checkboxInput(inputId = "deDot", label = "Overlay Dotplot?", value = FALSE)
                                ),
                                conditionalPanel(
                                  condition = "input.dePlotType == 'Violin Plot' || input.dePlotType == 'Split Violin Plot'",
                                  checkboxInput(inputId = "deBox", label = "Overlay Boxplot", value = TRUE)
                                ),
                                HTML('<hr>'),
                                selectInput(
                                  inputId = "miRplotColors",
                                  label = "Color Palette",
                                  choices =
                                    c(
                                      'Default',
                                      'Greys',
                                      'YellowOrangeRed',
                                      'YellowGreenBlue',
                                      'RedYellowGreen',
                                      'RedYellowBlue',
                                      'RedBlue',
                                      'Set1',
                                      'Set2',
                                      'Set3',
                                      'Spectral',
                                      'Rainbow',
                                      'RedGreenBlue',
                                      'RedBlackGreen',
                                      'Rainbow2.0',
                                      'Viridis',
                                      'Magma',
                                      'Inferno',
                                      'Plasma'
                                    ),
                                  selected = 'Greys'
                                ),
                                checkboxInput(
                                  inputId = "miRcolorsRev",
                                  label = "Reverse Palette?",
                                  value = FALSE
                                ),
                                actionButton(
                                  inputId = "plotDE",
                                  label = "Plot!",
                                  width = '100%'
                                )
                              )
                            ),
                            column(
                              9,
                              align = 'center',
                              plotOutput("singleMirBoxplot", height = 800, width = 950)
                            )
                          )
                  ),
                  tabItem(tabName = "Download",
                          fluidRow(
                            column(
                              8,
                              offset = 2,
                              align = "center",
                              br(),
                              wellPanel(
                                h4("File Downloads"),
                                textInput(
                                  inputId = 'hmImageName',
                                  label = 'Name for Heatmap Matrix File Download',
                                  value = "experimentName_heatmap",
                                  width = '75%'
                                ),
                                downloadButton('downloadHM', 'Download'),
                                HTML('<hr style = "width = 75%">'),
                                textInput(
                                  inputId = 'hmName',
                                  label = 'Name for Heatmap Matrix File Download',
                                  value = "experimentName_heatmapVals",
                                  width = '75%'
                                ),
                                downloadButton('downloadMatrix', 'Download'),
                                HTML('<hr style = "width = 75%">'),
                                textInput(
                                  inputId = 'corrhmImageName',
                                  label = 'Name for Correlation Heatmap Matrix File Download',
                                  value = "experimentName_corrheatmap",
                                  width = '75%'
                                ),
                                downloadButton('downloadCorrHM', 'Download'),
                                HTML('<hr style = "width = 75%">'),
                                textInput(
                                  inputId = 'circularPlotImageName',
                                  label = 'Name for Genomic Plot File Download',
                                  value = "experimentName_genomvis",
                                  width = '75%'
                                ),
                                downloadButton('downloadCircPlot', 'Download'),
                                HTML('<hr style = "width = 75%">'),
                                textInput(
                                  inputId = 'deName',
                                  label = 'Name for Differentially Expressed miRNA Table File Download',
                                  value = "experimentName_miRDETable",
                                  width = '75%'
                                ),
                                downloadButton('downloadDETab', 'Download'),
                                HTML('<hr style = "width = 75%">'),
                                textInput(
                                  inputId = 'fullName',
                                  label = 'Name for Complete miRNA Table File Download',
                                  value = "experimentName_miRFullTable",
                                  width = '75%'
                                ),
                                downloadButton('downloadFullTab', 'Download')
                              )
                            )
                          )
                  ),
                  tabItem(tabName = "Accuracy",
                          fluidRow(
                            HTML("<br>"),
                            column(
                              5,
                              wellPanel(
                                helpText("Ideal microRNA profile file should be a columnar .csv file, with rows corresponding to features. Up to five columns are accepted. The same formatting requirements carry over from the expression/conditions file. Row names should not be included."),
                                fileInput(
                                  inputId = "idealProfile",
                                  label = "Upload Synthetic Profile",
                                  width = '100%'
                                ),
                                selectInput(
                                  inputId = "idealCorrType",
                                  label = "Correlation Coefficient Method",
                                  choices = c("Pearson", "Kendall", "Spearman"),
                                  selected = "Pearson"
                                ),
                                HTML("<hr>"),
                                actionButton(inputId = "idealButton", label = "Process", width = '100%')
                              )
                            ),
                            column(
                              7,
                              align = "center",
                              wellPanel(
                                h4("Summary Statistics"),
                                tableOutput("idealProfileTable")
                              )
                            )
                          )
                  ),
                  tabItem(tabName = "Power",
                          fluidRow(
                            br(),
                            column(3,
                                   wellPanel(
                                     h4("Statistical Power Parameters"),
                                     helpText("Estimates relationship between sample size and statistical power."),
                                     conditionalPanel(
                                       condition = "output.dataExists",
                                       checkboxInput(
                                         inputId = "powerUseData",
                                         label = "Use uploaded data to help estimate power?",
                                         value = F
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.powerUseData && output.sigMirsExist",
                                       checkboxInput(
                                         inputId = "powerUseCutoff",
                                         label = "Use differentially expressed miRNA found in Differential Analysis for estimations?",
                                         value = F
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "!(input.powerUseData && input.powerUseCutoff)",
                                       numericInput(
                                         inputId = "powerFoldChange",
                                         label = "Minimum fold change of DE",
                                         value = 1,
                                         min = 0,
                                         max = 20,
                                         step = 0.05
                                       ),
                                       numericInput(
                                         inputId = "powerFDR",
                                         label = "False discovery rate (q-value)",
                                         value = 0.05,
                                         min = 0,
                                         max = 0.5,
                                         step = 0.001
                                       ),
                                       numericInput(
                                         inputId = "powerPrognosticFeatureCount",
                                         label = "Expected number of DE miRNA",
                                         value = 10,
                                         min = 0,
                                         step = 1
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.powerUseData && input.powerUseCutoff",
                                       helpText("This feature is experimental.")
                                     ),
                                     conditionalPanel(
                                       condition = "input.powerUseData",
                                       numericInput(
                                         inputId = "powerRepNumber",
                                         label = "# of estimation repetitions",
                                         value = 10,
                                         min = 5,
                                         step = 5
                                       ),
                                       helpText("Higher numbers mean better accuracy. For real analysis, it is recommended that this is higher than 100.")
                                     ),
                                     conditionalPanel(
                                       condition = "!input.powerUseData",
                                       numericInput(
                                         inputId = "powerFeatureCount",
                                         label = "Total number of miRNA",
                                         value = 1000,
                                         min = 0,
                                         step = 10
                                       ),
                                       numericInput(
                                         inputId = "powerAverageReadCounts",
                                         label = "Average read counts",
                                         value = 5,
                                         min = 0,
                                         step = 1
                                       ),
                                       numericInput(
                                         inputId = "powerMaxDispersion",
                                         label = "Dispersion of data",
                                         value = 0.5,
                                         min = 0,
                                         step = 0.01
                                       )
                                     )
                                   ),
                                   wellPanel(
                                     h4("Plotting Options"),
                                     conditionalPanel(
                                       condition = "input.powerUseData",
                                       helpText("If uploaded data is used, different plotting methods may produce different results.")
                                     ),
                                     selectInput(
                                       inputId = "powerPlotType",
                                       label = "Plotting Method",
                                       choices = c("Sample size", "Power")
                                     ),
                                     conditionalPanel(
                                       condition = "input.powerPlotType == 'Sample size'",
                                       helpText("Estimates a power for each sample size."),
                                       numericInput(
                                         inputId = "powerSampleSize",
                                         label = "Maximum number of samples",
                                         value = 250,
                                         min = 0,
                                         max = 1500,
                                         step = 10
                                       ),
                                       numericInput(
                                         inputId = "powerStepInterval",
                                         label = "Sample Size Interval",
                                         value = 25,
                                         min = 5,
                                         step = 5
                                       ),
                                       selectInput(
                                         inputId = "powerPlotChoice",
                                         label = "Point plot interval method",
                                         choices = c("Constant intervals", "Squared intervals", "Gradient-sensitive intervals")
                                       ),
                                       helpText("This setting affects which points get plotted."),
                                       conditionalPanel(
                                         condition = "input.powerPlotChoice == 'Gradient-sensitive intervals'",
                                         sliderInput(
                                           inputId = "powerGradientDetail",
                                           label = "Gradient sensitivity",
                                           min = 0,
                                           max = 100,
                                           value = 50
                                         ),
                                         helpText("0 = no effect, 100 = high detail in sloped regions")
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.powerPlotType == 'Power'",
                                       helpText("Plots points up to and near the desired power. Only points near the desired power are accurate."),
                                       numericInput(
                                         inputId = "powerMaxPower",
                                         label = "Desired Power",
                                         value = 0.9,
                                         min = 0.01,
                                         max = 1,
                                         step = 0.01
                                       ),
                                       helpText("Higher powers may take longer to estimate. Desired power cannot be 1.")
                                     )
                                   ),
                                   wellPanel(
                                     actionButton(
                                       inputId = "powerPlotButton",
                                       label = "Plot!",
                                       width = '100%'
                                     ),
                                     actionButton(
                                       inputId = "powerClearPlotButton",
                                       label = "Clear Plot",
                                       width = '100%'
                                     ),
                                     helpText("Up to 5 lines can be plotted at one time.")
                                   )
                            ),
                            column(9,
                                   plotOutput("powerPlot", height = 1000, width = '90%'),
                                   fluidRow(
                                     tableOutput("powerParameterTable"),
                                     tableOutput("powerTable1"),
                                     tableOutput("powerTable2"),
                                     tableOutput("powerTable3"),
                                     tableOutput("powerTable4"),
                                     tableOutput("powerTable5")
                                   )
                            )
                          )
                  )
                )
              )
)