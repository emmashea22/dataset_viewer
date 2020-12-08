
library(data.table)
library(shiny)
library(DT)
library(Sigfried)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(shinycssloaders)

library(heatmaply)
library(ggcorrplot)


#Outside of app to run---
#source code with functions: 
source('annotation_getWS_Data_code.R')

dataset_names <- c('TCGA', 'Target', 'GTEx', 'Blueprint', 'CoMMpass', 'CCLE', 'PTX')

#GDSx is the name of our internal database
gdsx_output_names <- c('omicsoft_tcga_gene_fpkm_2.1.0', 'omicsoft_target_gene_fpkm_1.0.0', 'omicsoft_gtex_gene_fpkm_1.0.0', 'omicsoft_blueprint_gene_fpkm_1.0.0', 'MMRF_CoMMpass_TPM_rnaseq_1.0.0', 'ccle.rnaseq.gene.tpm_1.0.0', 'ptx.rnaseq.gene.tpm_5.1.0')
gdsx_input_names <- c('omicsoft_tcga_gene_fpkm', 'omicsoft_target_gene_fpkm', 'omicsoft_gtex_gene_fpkm', 'omicsoft_blueprint_gene_fpkm', 'MMRF_CoMMpass_TPM_rnaseq', 'ccle.rnaseq.gene.tpm', 'ptx.rnaseq.gene.tpm')
datasets <- data.table(dataset_names, gdsx_output_names, gdsx_input_names)

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)

sample_annotation <- readRDS('iobp_subset_noNormals_annotation.RDS')


### To run app ---

ui <- fluidPage(
  
  titlePanel("GDSx Visualizations"),
  
  sidebarLayout(
    sidebarPanel(
      textInput('gene', 'Input gene symbol(s):', value = "", width = NULL, placeholder = 'e.g. TP53'),
      textInput('geneSig', 'Custom Gene Signature:', value = "", width = NULL, placeholder = 'e.g. CD274, CD8A, PDCD1'),
      textInput('sigName', 'Name Your Signature:', value="", width=NULL, placeholder= 'e.g. PDL1_PD1_CD8_Sig'),
      radioButtons(inputId="sigsOnly", label="Only Display Signature Data?", 
                   choices=c("Yes","No"),selected="No"),
      checkboxGroupInput("DataInput", "Select data types:",
                         c("TCGA" = "TCGA",
                           "Target" = "Target",
                           "GTEx" = "GTEx",
                           "Blueprint" = "Blueprint",
                           "CoMMpass" = "CoMMpass",
                           "CCLE" = "CCLE",
                           "PTX" = "PTX")),
      actionButton("goButton", strong("Retrieve Data")),
      h5(em(strong("Warning:"),"Requesting lots of genes across many datasets will affect load time."), style='color:red')),
    mainPanel(
      tabsetPanel(
        type="tabs",
        tabPanel("Data Table", withSpinner(uiOutput("test_table", height = "600px"), type=6)),
        tabPanel("Primary Site Scatter", 
                 br(),
                 h5(em("All gene expression values are Log2(TPM+1)")),
                 br(),
                 withSpinner(plotlyOutput("testplot1", height = "800px"), type=6)),
        tabPanel("Primary Site Boxplot",
                 br(),
                 h5(em("All gene expression values are Log2(TPM+1)")),
                 br(),
                 withSpinner(plotlyOutput("testplot2", height = "800px"), type=6)),
        tabPanel("Correlation Matrix",
                 br(),
                 h5(em("All gene expression values are Log2(TPM+1)")),
                 br(),
                 withSpinner(plotlyOutput("corr", height = "800px"), type=6))
      )
    )
    
  ) )



server <- function(input, output) {
    
    #this creates the gene list that is sent to the GDSx WS to create a custom query
    genesToFind<- reactive({
    if (input$goButton == 0)
      return()
    isolate({
      if (!is.null(input$gene)) {
        genes <- input$gene
        genes <- trimws(strsplit(genes, ",")[[1]])
      }
      else {
        genes <- c("")
      } 
      if (!is.null(input$geneSig)) {
        custGenes <- input$geneSig
        custGenes <- trimws(strsplit(custGenes, ",")[[1]])
        sigName <- input$sigName
      } else {
        custGenes <- c("")
      }
      
      allSigGenes <- append(custGenes)
      geneList <- append(genes, allSigGenes)
      genesToFind <- subset(gencode_annotation, GENE_SYMBOL %in% geneList)
    })
  })
  
  #this creates the dataset list that is sent to the GDSx WS to create a custom query
  datasets1 <-reactive({
    if(input$goButton == 0)
      return()
    datasets1<-datasets[(datasets$dataset_names %in% input$DataInput), ]
  })
  
  test_data<-reactive({
    if (input$goButton == 0)
      return()
    
    test_data <- read.csv(getGdsxWsCall(genesToFind(), datasets1()))
  })
  
  #this is the dataset that is presented in the table view
  final_datset <- reactive({
    if (input$goButton == 0)
      return()
    
    req(input$DataInput)
    
    print(genesToFind())
    datasets <- datasets1()
    dataset <- subset(test_data(), OBSERVATION_NAME %in% datasets$gdsx_output_names)
    print(head(dataset))
    
    # We have datasets that have been normalized differently (i.e. FPKM vs TPM), this function will make sure every dataset is normalized the same way
    normalized_data <- normalizeDatasets(dataset, genesToFind())
    
    print(input$gene)
    print(input$geneSig)
    print(input$sigName)
        
    if (input$geneSig != "") {
      custGenes <- input$geneSig
      custGenes <- trimws(strsplit(custGenes, ",")[[1]])
      sigName <- input$sigName
      cust_signature <- addCustomSignature(normalized_data, custGenes, sigName)
    } else {
      cust_signature <- data.frame()
    }
   
    final_data <- getFinalDataset(normalized_data, cust_sig_calc = cust_signature)
   
  })
  
  
  
  output$test_table<-renderTable({
    
    final_datset()
    
  })
  

  #this output is used for the plots
  geneList<- reactive({
    if (input$goButton == 0)
      return()
    if (!is.null(input$gene)) {
      genes<-input$gene
      genes <- trimws(strsplit(genes, ",")[[1]])
    }
    else {
      genes <- c("")
    } 
    if (!is.null(input$geneSig)) {
      custGenes <- input$geneSig
      custGenes <- trimws(strsplit(custGenes, ",")[[1]])
      sigName <- input$sigName
    } else {
      custGenes <- c("")
      sigName <- c("")
    }
    
    allSigGenes <- append(custGenes)
    allGenes <- append(genes, allSigGenes)
    sigNames <- append(sigName, definedName)
    geneList <- append(allGenes, sigNames)
    geneList <- geneList[geneList != ""]
    geneList <- gsub("-", "_", geneList)
  })
  
  output$testplot1<-renderPlotly({
    if (input$goButton == 0)
      return()
    
    if (!is.null(input$gene)) {
      genes<-input$gene
      genes <- trimws(strsplit(genes, ",")[[1]])
    }
    else {
      genes <- c("")
    } 
    
    if (!is.null(input$geneSig)) {
      custGenes <- input$geneSig
      custGenes <- trimws(strsplit(custGenes, ",")[[1]])
      sigName <- input$sigName
    }
    
    sigNames <- append(definedName)
    sigNames <- sigNames[sigNames != ""]
    
    data <- final_datset()
    names(data) <- gsub("-", "_", names(data))
    vars <- geneList()
    vars <- intersect(names(data), vars)
    
    sigOnlyChoice <- input$sigsOnly
    if(sigOnlyChoice=="Yes") {
      vars <- append(genes, sigNames)
    }
    vars <- unique(vars)
    
    
    plots1 <- lapply(vars, function(var) {
      scatterplot <- plot_ly(data = data, 
                              x = ~Source_Primary_Site_Tumor_Or_Normal,
                              y =  as.formula(paste0("~", var)),
                              type = "scatter",
                              mode = "markers",
                              marker = list(size=5),
                              color = ~Primary_Site_Tumor_Or_Normal,
                              colors = brewer.pal(length(unique(data$Primary_Site_Tumor_Or_Normal)), "Accent"),
                              legendgroup = ~Primary_Site_Tumor_Or_Normal,
                              text = ~paste('<b>Gene:</b>', var, '<br><b>Log2(TPM+1):</b>', data[,var], '<br><b>Sample:</b> ',data$SubjectID, '<br><b>Subtype:</b>', data$Source_Primary_Site_Histology),
                              hoverinfo = "text")
      })
    
    subplot(plots1, nrows=length(plots1), shareX=TRUE, titleX=FALSE, titleY = TRUE) %>% 
      layout(xaxis = list(tickfont = list(size = 10), tickangle=-45))
 
  })
  
  output$testplot2<-renderPlotly({
    if (input$goButton == 0)
      return()
    
    data <- final_datset()
  
    if (!is.null(input$gene)) {
      genes<-input$gene
      genes <- trimws(strsplit(genes, ",")[[1]])
    }
    else {
      genes <- c("")
    } 
    
    if (!is.null(input$geneSig)) {
      custGenes <- input$geneSig
      custGenes <- trimws(strsplit(custGenes, ",")[[1]])
      sigName <- input$sigName
    }
    
    sigNames <- append(sigName)
    sigNames <- sigNames[sigNames != ""]
  
    data <- final_datset()
    names(data) <- gsub("-", "_", names(data))
    
    vars <- geneList()
    vars <- intersect(names(data), vars)
    
    sigOnlyChoice <- input$sigsOnly
    if(sigOnlyChoice=="Yes") {
      vars <- append(genes, sigNames)
    }
    vars <- unique(vars)
    
    plots2 <- lapply(vars, function(var) {
      boxplot_plot <- plot_ly(data = data, 
                              x = ~Source_Primary_Site_Tumor_Or_Normal,
                              y =  as.formula(paste0("~", var)),
                              type = "box",
                              #mode = "markers",
                              marker = list(size=5),
                              color = ~Source_Primary_Site_Tumor_Or_Normal,
                              colors = brewer.pal(length(unique(data$Source_Primary_Site_Tumor_Or_Normal)), "Accent"),
                             legendgroup = ~Primary_Site,
                              text = ~paste('<b>Study Source:</b>', data$Study_Source, '<b>Gene:</b>', var, '<br><b>Log2(TPM+1):</b>', data[,var], '<br><b>Sample:</b> ',data$SampleID, '<br><b>Subtype:</b>', data$Source_Primary_Site_Histology),
                              hoverinfo = "text") 
    })
    
  
    
   subplot(plots2,nrows=length(plots2), shareX=TRUE, titleX=FALSE, titleY = TRUE) %>% 
      layout(
        xaxis = list(tickfont = list(size = 10)))
  })
  
  output$corr<-renderPlotly({
    
    if (input$goButton == 0)
      return()
  
    if (!is.null(input$gene)) {
      genes<-input$gene
      genes <- trimws(strsplit(genes, ",")[[1]])
    }
    else {
      genes <- c("")
    } 
    
    if (!is.null(input$geneSig)) {
      custGenes <- input$geneSig
      custGenes <- trimws(strsplit(custGenes, ",")[[1]])
    }
   
    sigNames <- append(sigName)
    sigNames <- sigNames[sigNames != ""]
    
    data <- final_datset()
    names(data) <- gsub("-", "_", names(data))
    
    vars <- geneList()
    vars <- intersect(names(data), vars)
    
    sigOnlyChoice <- input$sigsOnly
    if(sigOnlyChoice=="Yes") {
      vars <- append(genes, sigNames)
    }
    vars <- unique(vars)
    data <- final_datset()
    data1<- data[vars]

    
    heatmaply_cor(
      cor(data1),
      xlab = "Features", 
      ylab = "Features"
    )
    
  })
  
  
  
}
shinyApp(ui = ui, server = server)
