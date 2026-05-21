library(shiny)
library(DESeq2)
library(DT)
library(shinyjs)
library(bslib)
library(ggplot2)
library(shinycssloaders)
library(dplyr)
library(EnhancedVolcano)
library(shinyWidgets)

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = bs_theme(bootswatch = "slate"), # Example theme
  titlePanel(
    div(
      style = "text-align: center; color: maroon;",  # Center and apply the color
      tags$b("DESeq2 Differential Expression Analysis")
    )
  ),
  useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      fileInput("counts", "Upload Count Matrix", accept = c(".csv", ".tsv", ".txt")),
      fileInput("meta", "Upload Metadata", accept = c(".csv", ".tsv", ".txt")),
      uiOutput("condition_col_ui"), # Condition column will be a dynamic dropdown
      uiOutput("control_selector"),
      uiOutput("comparison_selector"),
      actionButton("run_deseq2", "Run DESeq2"),
      ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("DESeq2 Results", DTOutput("results_table") %>% withSpinner(color = "#0dc5c1"),
                 downloadButton("download_table", "Download Results")),
        tabPanel("PCA Plot", plotOutput("PCA_plot") %>% withSpinner(color = "#0dc5c1"),
                 downloadButton("download_PCA", "Download PCA plot")),
        tabPanel("MA Plot", plotOutput("MA_plot") %>% withSpinner(color = "#0dc5c1"),
                 downloadButton("download_MA", "Download MA plot")),
        tabPanel("Volcano Plot", plotOutput("volcano_plot") %>% withSpinner(color = "#0dc5c1"),
                 downloadButton("download_volcano", "Download Volcano plot"))
      )
      
    )
  )
)
    

# Define server logic 
server <- function(input, output, session) {
  
  # Helper to read CSV, TSV, or TXT with separator detection
  read_delim_auto <- function(file_path) {
    ext <- tools::file_ext(file_path)
    
    if (ext == "csv") {
      read.csv(file_path, row.names = 1, check.names = FALSE)
    } else {
      read.delim(file_path, row.names = 1, check.names = FALSE)
    }
  }
  
  # Reactive values
  meta_data <- reactiveVal(NULL)
  dds_obj <- reactiveVal(NULL)
  conditions <- reactiveVal(NULL)
  
  
  # Step 1: Load metadata and get column names
  observeEvent(input$meta, {
    meta <- read_delim_auto(input$meta$datapath)
    meta_data(meta)
    
    # Dynamically create the condition column dropdown
    output$condition_col_ui <- renderUI({
      selectInput("condition_col", "Select Condition Column:", choices = colnames(meta))
    })
    
    # Reset the conditions in case metadata changes
    conditions(NULL)
  })
  
  # Step 2: Get unique conditions once condition column is selected
  observeEvent(input$condition_col, {
    req(input$condition_col)  # Ensure input is available
    
    meta <- meta_data()
    cond_col <- input$condition_col
    
    if (!cond_col %in% colnames(meta)) {
      showNotification("Condition column not found in metadata.", type = "error")
      return()
    }
    
    # Get unique conditions from the selected column
    unique_conds <- unique(meta[[cond_col]])
    conditions(unique_conds)
    
    # Create control selector UI
    output$control_selector <- renderUI({
      selectInput("control", "Select Control Group:", choices = unique_conds)
    })
  })
  
  # Step 3: Create pairwise comparison selector once control is selected
  observeEvent(input$control, {
    req(conditions())
    comps <- setdiff(conditions(), input$control)
    
    output$comparison_selector <- renderUI({
      selectInput("comparison", "Select Comparison vs Control:",
                  choices = comps)
    })
  })
  
  # Step 4: Run DESeq2
  observeEvent(input$run_deseq2, {
    withProgress(
    message = "Performing DE Analysis...",
    detail = "Initializing...",
    value = 0,
    min = 0,
    max = 100,
    {
      for (i in 1:10) {
        Sys.sleep(0.5) # Simulate work
        incProgress(10, detail = paste("Step", i, "of 10"))
      }
      Sys.sleep(1)
      
     req(input$counts, input$meta, input$comparison, input$control)
    
    counts <- read_delim_auto(input$counts$datapath)
    meta <- meta_data()
    cond_col <- input$condition_col
    
    # Ensure samples match
    common_samples <- intersect(colnames(counts), rownames(meta))
    counts <- counts[, common_samples]
    meta <- meta[common_samples, , drop = FALSE]
    # Step 2
   #incProgress(2/3, message = "Processing data...")
    #Sys.sleep(20) # Simulate work
    
    # Ensure factors
    meta[[cond_col]] <- factor(meta[[cond_col]])
    
    # Set reference/control level
    meta[[cond_col]] <- relevel(meta[[cond_col]], ref = input$control)
    
    # Simulate the long DESeq2 computation
   #Sys.sleep(10)  # Simulate long computation time (remove in real code)
    
    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                                  colData = meta,
                                  design = as.formula(paste("~", cond_col)))
    
    
    dds <- DESeq(dds)
    # Step 3
    #incProgress(3/3, message = "Processing data...")
    #Sys.sleep(10) # Simulate work
    
    dds_obj(dds)
    
    setProgress(100, message = "DE Analysis complete!", detail = "Finished.")
    })
    
 }) 
  
  # Step 5: Display results for selected comparison
  output$results_table <- renderDT({
    req(dds_obj(), input$comparison, input$control)
    
    dds <- dds_obj()
    cond_col <- input$condition_col
    res <- results(dds, contrast = c(cond_col, input$comparison, input$control))
    
    # Filter out rows with NA padj values
    res <- res[!is.na(res$padj), ]
    res <- res[order(res$padj),]
    # Render the results as a datatable
    datatable(as.data.frame(res), options = list(pageLength = 10))

  })
  
  
  #Add a download handler for the results table
  output$download_table <- downloadHandler(
    filename = function() {
      req(input$comparison, input$control)
      cond_col <- input$condition_col
      contrast = paste(cond_col, input$comparison, input$control,sep="_")
      paste("DESeq2_results_",contrast,"_",Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(dds_obj(), input$comparison, input$control)
      dds <- dds_obj()
      cond_col <- input$condition_col
      
      # Generate results for the selected contrast
      res <- results(dds, contrast = c(cond_col, input$comparison, input$control))
      res <- res[!is.na(res$padj), ]
      resOrdered <- res[order(res$padj),]
      write.csv(resOrdered, file, row.names = T)
    }
  )
  ############################################################################################################
  #########################. PCA PLOT ########################################################################
  ############################################################################################################
  
  #Plot PCA in a new tab
  output$PCA_plot <- renderPlot({
    req(dds_obj())
    dds <- dds_obj()
    rld <- rlog(dds, blind = TRUE)
    plotPCA(rld, intgroup = input$condition_col)
  })

  #Add a download handler for the PCA plot
  output$download_PCA <- downloadHandler(
    filename = function() {
      req(input$comparison, input$control)
      cond_col <- input$condition_col
      contrast = paste(cond_col, input$comparison, input$control,sep="_")
      paste("PCA_plot_",contrast,"_",Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      req(dds_obj())
      dds <- dds_obj()
      rld <- rlog(dds, blind = TRUE)
      PCA <- plotPCA(rld, intgroup = input$condition_col)
      print(PCA)
      # Close the device
      dev.off()
    },
    contentType = "application/pdf"
  )
  ############################################################################################################
  #########################. MA PLOT ########################################################################
  ############################################################################################################
  
  output$MA_plot <- renderPlot({
    req(dds_obj(), input$comparison, input$control)
    dds <- dds_obj()
    cond_col <- input$condition_col
    
    # Generate results for the selected contrast
    res <- results(dds, contrast = c(cond_col, input$comparison, input$control))

    contrast = paste(cond_col, input$comparison, input$control,sep="_")
    # Plot MA plot
    plotMA(res, ylim = c(-5, 5), main = paste("MA Plot for",contrast,sep=" "))
  })
  
  
  #Add a download handler for the MA plot
  output$download_MA <- downloadHandler(
    filename = function() {
      req(dds_obj(), input$comparison, input$control)
      dds <- dds_obj()
      cond_col <- input$condition_col
      
      # Generate results for the selected contrast
      contrast = paste(cond_col, input$comparison, input$control,sep="_")
      paste("MA_plot_",contrast,"_",Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      req(dds_obj(), input$comparison, input$control)
      dds <- dds_obj()
      cond_col <- input$condition_col
      contrast = paste(cond_col, input$comparison, input$control,sep="_")
      # Generate results for the selected contrast
      res <- results(dds, contrast = c(cond_col, input$comparison, input$control))
      MA <- plotMA(res, ylim = c(-5, 5), main = paste("MA Plot for",contrast,sep=" "))
      print(MA)
      # Close the device
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  ############################################################################################################
  #########################. Volcano PLOT ####################################################################
  ############################################################################################################
  output$volcano_plot <- renderPlot({
    
    req(dds_obj(), input$comparison, input$control)
    dds <- dds_obj()
    cond_col <- input$condition_col
    contrast = paste(cond_col, input$comparison, input$control,sep="_")
    # Generate results for the selected contrast
    res <- results(dds, contrast = c(cond_col, input$comparison, input$control))
    res <- res[!is.na(res$padj), ]
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'padj',title=paste("Volcano Plot for",contrast,sep=" "),
                    pointSize = 4.0,
                    legendLabels=c('Not sig.','Log (base 2) FC','adjP-value',
                                   'adjP-value & Log (base 2) FC'),
                    legendLabSize = 16,
                    legendIconSize = 5.0,
                    xlab = bquote(~Log[2]~ 'fold change'),
                    ylab = bquote(~-Log[10]~ 'adj P-value'),
                    gridlines.major = FALSE,
                    gridlines.minor = FALSE)
  })
  #   ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
  #     geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1) +
  #     scale_color_manual(values = c("gray", "red")) +
  #     theme_minimal() +
  #     xlim(c(-5, 5)) + 
  #     ylim(c(0, 10)) +
  #     labs(title = paste("Volcano Plot for",contrast,sep=" "), x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")
  # })
  
  #Add a download handler for the Volcano plot
  output$download_volcano <- downloadHandler(
    filename = function() {
      req(dds_obj(), input$comparison, input$control)
      dds <- dds_obj()
      cond_col <- input$condition_col
      
      # Generate results for the selected contrast
      contrast = paste(cond_col, input$comparison, input$control,sep="_")
      paste("Volcano_plot_",contrast,"_",Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      req(dds_obj(), input$comparison, input$control)
      dds <- dds_obj()
      cond_col <- input$condition_col
      contrast = paste(cond_col, input$comparison, input$control,sep="_")
      # Generate results for the selected contrast
      res <- results(dds, contrast = c(cond_col, input$comparison, input$control))
      volcano <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
                        geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1) +
                        scale_color_manual(values = c("gray", "red")) +
                        theme_minimal() +
                        xlim(c(-5, 5)) + 
                        ylim(c(0, 10)) +
                        labs(title = paste("Volcano Plot for",contrast,sep=" "), x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")
      print(volcano)
      # Close the device
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
