library(shiny)

shinyUI(pageWithSidebar(

  headerPanel("Integrative expression analysis for RNAseq data"),
  sidebarPanel(

conditionalPanel(
  condition = "input.seltab == 'HM'", 
    selectInput("order", "Heatmap column order:", 
                selected="No",
                list("No rearrangement" = "No","Cluster Columns" = "hc","Heatmap PC1" = "HM1","Heatmap centroid" = "HMC")
    ),  
                       
   	selectInput("list2", "RNAs to use:",
                list("List of genes (need ≥ 2)" = "user",
                     "Differential expression analysis" = "Differential_expression",
                     "Key Immune Gene sets" = "Key_Immune"
                )
    ),
                     
    conditionalPanel(
      condition = "input.list2 == 'Key_Immune'",
      selectInput("list", "Immune gene sets to use:",  choices)
    ),
               
    conditionalPanel(
      condition = "input.list2 == 'user'",
      textInput("genes","Gene list",value = "IL1A IL1B IL1F6 IL1R1 IL1RN")
    )
),


conditionalPanel(
      condition = "input.list2 == 'Differential_expression'",
      
selectInput("Anal.Type", "What statistical framework?",list("Limma" = "Limma", "EdgeR" = "EdgeR","Limma Voom" = "LV")),
selectInput("wc", "What type of comparison?", as.list(eg4)),      
numericInput("pSlider", "p value", 0.05),
numericInput("max.number", "Number of genes to show, working from top of list", 100),   
numericInput("lfcSlider", "log fold change", 1),   
selectInput("mtc", "Multiple testing correction:",list("None" = "none", "Benjamini-Hochberg" = "BH")),
         
                conditionalPanel(
      condition = "input.wc == 'bespoke'",                
      textInput("bsc","Enter your own contrast", 
                value = "group1 - (group2 + group3)")                       
                )),
                
 br(),   


conditionalPanel(
  condition = "input.seltab == 'QC'",  
  radioButtons("ctype","Correlation type?", c("Pearson", "Spearman"), selected="Pearson")
),

conditionalPanel(
  
condition = "input.seltab == 'DT'",  
  radioButtons("ctype","Correlation type?", c("Pearson", "Spearman"), selected="Pearson"),
  numericInput("trimmer.min", "Enter MINimum mean read number for a gene:", 10),
  numericInput("trimmer.max", "Enter MAXimum mean read number for a gene:", 20000),
  selectInput("norm", "Normalise?",
                list("None" = "none",
                "Quantile (limma)" = "q",
                "TMM (EdgeR)" = "t",
                "RLE (EdgeR)" = "r",
                "Upper Quartile (EdgeR)" = "uq")),
  selectInput("l2", "Log2 transform the data? (not if using EdgeR or limma with Voom)",list("No" = "No", "Yes" = "Yes")),
radioButtons("expl","Show explanations?", c("No", "Yes"), selected="No")
),

conditionalPanel(
condition = "input.seltab == 'HM'",
         radioButtons("Heatmap_control","Show Heatmap controls?",
                     c("Yes", "No"),
                     selected="No")),
                     
                                           br(),
                     textInput("fn", "enter prefix for saved filename if you wish", value=""),
 actionButton("sv", "Save"),
    h6(p("Save file of selected genes and data")) ,  
                     
conditionalPanel(
      condition = "input.Heatmap_control == 'Yes'", 
    selectInput("colors", "Color scheme:",
                list("Green/Red" = "greenred",
                     "Blue/Red" = "bluered",
                     "Blue/Yellow" = "blueyellow",
                     "Matlab" = "matlab")),
                     
br(), 

        radioButtons("radioCLM","Hierarchical clustering method for heatmap rows (and for columns if selected)",
                     c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                     selected="ward.D"),
      
        radioButtons("radioDist","Distance metric",
                     c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                     selected="euclidean"),
                     
 br(), 
      
        radioButtons("radioScale","Scaleing",
                     c("none", "row", "col"),
                     selected="row")              
                     
                     )
    ,width=2
    ),

   # use an ID in the tabset panel to use in server logic to determine which of the current tabs is active
  mainPanel(
    tabsetPanel(id="seltab",
      # add values to the tabpanels so we can see which one is currently active
      tabPanel(h4("QC"), plotOutput("s.c", height = 600, width = 800),plotOutput("pca", height = 600, width = 800),value="QC"),
      tabPanel(h4("Trimming"), textOutput("diag.expl"),   plotOutput("m.hist", height = 600, width = 800), plotOutput("m.sd", height = 600, width = 800), plotOutput("box.pl", height = 600, width = 800), plotOutput("box.pll", height = 600, width = 800), value="DT"),
      tabPanel(h4("Expression Heatmap"), textOutput("dif.expl"),plotOutput("expMap", height = 1000, width = 1000), value="HM"),
      tabPanel(h4("GO gene sets and pathways"), dataTableOutput("t.go"), value="GO"),
      tabPanel(h4("Transcription Factor sets"), dataTableOutput("t.tf"), value="TF")
    ),width=10
  )
)
)

