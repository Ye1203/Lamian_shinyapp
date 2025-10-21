dataVisualizationUI <- function() {
  fluidPage(
    shinyjs::useShinyjs(),
    tags$head(
      tags$style(HTML("
        .description-panel {
          background: #f8f9fa; 
          padding: 15px; 
          border-radius: 5px; 
          border: 1px solid #ddd;
          margin-bottom: 20px;
          transition: opacity 0.3s ease;
        }
      "))
    ),
    fluidRow(
      column(12, tags$label("Result folder"))
    ),
    fluidRow(
      column(6, textInput("result_folder", NULL, value = "", width = "100%")),
      column(2, actionButton("read_xde", "READ", style = "width: 100%")),
      column(4, downloadButton("download_stat", "Download statistical result", style = "width: 100%"))
    ),
    hr(),
    uiOutput("dynamic_content")
  )
}

dataVisualization_server <- function(input, output, session) {
  data_loaded <- reactiveVal(FALSE)
  xde_result <- reactiveVal(NULL)
  stat <- reactiveVal(NULL)
  sce <- reactiveVal(NULL)
  gene_analysis_selected <- reactiveVal(NA)
  gene_plot1_NA <- reactiveVal(NULL)
  gene_plot2_NA <- reactiveVal(NULL)
  gene_plot1 <- reactiveVal(NULL)
  gene_plot2 <- reactiveVal(NULL)
  gene_plot3 <- reactiveVal(NULL)
  gene_plot4 <- reactiveVal(NULL)
  show_points <- reactiveVal(FALSE)
  included_zero <- reactiveVal(FALSE)
  gene_list <- reactiveVal()
  heatmap_plot <- reactiveVal(NULL)
  heatmap_data <- reactiveVal(NULL)
  heatmap_name <- reactiveVal(NULL)
  result_visulization_folder <- reactiveVal(NULL)
  observe({
    if (nchar(input$result_folder) > 0) {
      shinyjs::enable("read_xde")
      shinyjs::enable("download_stat")
    } else {
      shinyjs::disable("read_xde")
      shinyjs::disable("download_stat")
    }
    
    if (is.null(heatmap_plot())) {
      shinyjs::disable("save_heatmap_plot")
    } else {
      shinyjs::enable("save_heatmap_plot")
    }
    
    if (is.null(heatmap_data())) {
      shinyjs::disable("save_heatmap_data")
    } else {
      shinyjs::enable("save_heatmap_data")
    }
  })
  
  observeEvent(input$read_xde, {
    showModal(modalDialog(
      title = "Please wait",
      tagList(
        p("Reading data..."),
        p("The waiting time is related to the size of data. Usually takes several minutes.")
      ),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      read_result <- visualization_read(input$result_folder)
      xde_result(read_result$xde_result)
      stat(read_result$stat)
      sce(read_result$sce)
      data_loaded(TRUE)
      removeModal()
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      error_msg <- paste0("Reading data error:\n", e$message)
      
      if (grepl("cannot open the connection", e$message, ignore.case = TRUE)) {
        error_msg <- paste0(error_msg, "\nThis is usually caused by typing the wrong path or your lamian program not being completed. Please ensure:\n",
                            "1. Your Lamian program has completed running.\n",
                            "2. The target folder contains 'xde_result.RDS' file.\n",
                            "3. The folder path is correct.")
      }
      
      showModal(modalDialog(
        title = "Error",
        HTML(gsub("\n", "<br/>", error_msg)),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
  })
  
  output$dynamic_content <- renderUI({
    if (!data_loaded()) {
      column(10,
             tags$div(
               style = "background: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #ddd;",
               h4("Please read data first!!!", style = "color: red;"),
               p("1. Enter the 'result folder' you selected previously."),
               p("2. click 'READ'."),
               hr(),
               h4("Download statistical result (Optional)"),
               p("1. Enter the 'result folder' you selected previously."),
               p("2. click 'Download statistical result'."),
               hr(),
               h4("Lamian_statistical_result Introduction."),
               p("The first column contains gene names."),
               p("Column name prefixes:"),
               tags$ul(
                 tags$li(strong("fdr"), " - False discovery rate adjusted p-values (Benjamini-Hochberg correction).", 
                         span("[Recommend]", style = "color:red; font-weight:bold;")),
                 tags$li(strong("pval"), " - Raw p-values indicating statistical significance of differences between conditions."),
                 tags$li(strong("z"), " - Z-scores quantifying standardized deviation from null distribution."),
                 tags$li(strong("log.pval"), " - Log-transformed p-values using kernel density estimation.")
               ),
               p("Column name suffixes:"),
               tags$ul(
                 tags$li(strong("overall"), " - Combined test considering both meanDiff and trendDiff effects."),
                 tags$li(strong("meanDiff"), " - Tests for absolute expression differences (vertical curve shifts)."),
                 tags$li(strong("trendDiff"), " - Tests for pattern/shape differences (slope changes).")
               ),
               p("Analysis workflow:"),
               tags$ol(
                 tags$li("First performs overall testing (fdr.overall)."),
                 tags$li("Only genes with fdr.overall < 0.05 proceed to meanDiff/trendDiff testing."),
                 tags$li("Three possible significant patterns: mean changes only, trend changes only, or both.")
               )
             )
      )
    } else {
      tagList(
        tags$h5("1. Use \"HETAMAP\" to draw heatmap and copy the gene of the gene pattern you are interested in.", style = "color:grey;"),
        tags$h5("2. Put the gene of interest into \"MULTI GENE ANALYSIS\" and save all the results", style = "color:grey;"),
        tags$h5("3. Use \"DOWNLOAD VISUALIZATION\" download all results to local.", style = "color:grey;"),
        tags$h5("4. For genes you want to further distribute use \"SINGLE GENE ANALYSIS\" and then use \"DOWNLOAD VISUALIZATION\" download all results to local.", style = "color:grey;"),
        navlistPanel(
          id = "navlist",
          widths = c(2, 10),
          
          # Information Tab,
          tabPanel("PARAMETER SETTING",
                   fluidRow(
                     column(12,
                            uiOutput("parameterDisplay")
                     )
                   )
          ),
          # HEATMAP Tab
          tabPanel("HEATMAP", 
                   tags$div(
                     id = "heatmap_desc", class = "description-panel",
                     style = "display: none;",
                     h3("HEATMAP Description"),
                     h4("Draw a heatmap for specified genes.", style = "color: grey;"),
                     hr(),
                     h4("1. Enter the genes list"),
                     p("The default gene list is for genes with an fdr.overall < 0.05. You can also enter your own gene list, separated by line breaks (you can directly copy and paste multiple rows of gene names from the Excel file downloaded using \"Download statistical results\")."),
                     p("a. Click \"",
                       tags$span(style = "color:#4678B2; font-weight:bold", "fdr.overall < 0.05"),
                       "\" to use the default list of genes with an fdr.overall value < 0.05."),
                     p("b. Click \"",
                       tags$span(style = "color:#C95C54; font-weight:bold", "Clear All"),
                       "\" to clear all genes."),
                     tags$hr(),
                     h4("2. Draw heatmap"),
                     p("a. ", tags$b("Clustering Method."), " Clustering Method includes Hierarchical clustering and K-means clustering. You can click to select."
                     ),
                     p("b. ", tags$b("Scale By."), " Scale By has two options: Both groups together and Each group separately. Selecting Both groups together will combine data from two conditions first and then perform z-scaling, allowing you to observe the absolute expression difference between the two conditions. Selecting Each group separately will perform z-scaling separately on each condition's data before merging, which helps observe trend differences between the two conditions. The final choice affects the expression values (from low to high) shown in the heatmap for the two conditions."
                     ),
                     p("c. ", tags$b("Cluster By."), " Cluster By has two options: Both groups together and Each group separately. It uses the scaled data obtained previously to perform classification, impacting the clustering distribution shown on the left side of the heatmap (and also changes the gene order)."
                     ),
                     p("d. ", tags$b("Select Color."), " Select color consists of three parts: low, mid, and high. You can click on the respective color box to reselect the color. The default range is [-2, 2], meaning scaled values below -2 will be colored as Low color, above 2 as High color, and 0 as Mid color."
                     ),
                     p("e. ", tags$b("Number of Clusters."), " Number of Clusters is used to select the number of clusters. When the Clustering Method is Hierarchical, Number of Clusters only changes the cluster distribution on the left side without changing gene order. When the Clustering Method is K-means, Number of Clusters changes both the cluster distribution and gene order (since K-means recalculates based on the number of clusters)."
                     ),
                     p("f. ", tags$b("Draw Plot."), " After setting parameters, click the \"",
                       tags$span(style = "color:#74B666; font-weight:bold", "Draw Plot"),
                       "\" button to draw the heatmap. Please wait patiently for rendering to complete (after the heatmap appears) before performing other operations."
                     ),
                     p("g. ", tags$b("Save Plot."), " You can adjust each parameter repeatedly and draw different heatmaps. When you get a satisfactory heatmap, click the \"",
                       tags$span(style = "color:#84A4CA; font-weight:bold", "Save Plot"),
                       "\" button. The image will be saved to the Result folder at visualization_YYYYmmddHHMM/HEATMAP, including two pdf files with and without gene names. The naming convention is whether it `includes gene names`_`Clustering Method`.`Number of Clusters`_`Scale By`.`Cluster By`."
                     ),
                     hr(),
                     h4("3. Download heatmap data"),
                     p("a. ", tags$b("Number of Pseudotime samples."), " Set the number of sampling points. The maximum value is 1000 and the minimum value is 3. Heatmap data will be generated based on the selected number of points."),
                     p("b. ", tags$b("Save Data."), "Click the \"",
                       tags$span(style = "color:#95D3E6; font-weight:bold", "Save Data"),
                       "\" button to save the final data. This includes gene names, their assigned clusters, scaled data, and raw data. The file will be saved to the Result folder at visualization_YYYYmmddHHMM/HEATMAP folder, named as Data_`Clustering Method`.`Number of Clusters`_`Scale By`.`Cluster by`.xlsx")
                   ),
                   uiOutput("heatmap_ui")
          ),
          
          
          # MULTI GENE ANALYSIS Tab
          tabPanel("MULTI GENE ANALYSIS", 
                   tags$div(
                     id = "multi_gene_analysis_desc", class = "description-panel",
                     style = "display: none;",
                     h3("Multi Gene Analysis Description"),
                     h4("Directly save the visualization analysis results of specified genes.", style = "color: grey;"),
                     hr(),
                     h4("1. Enter the genes list"),
                     p("The default gene list is for genes with an fdr.overall < 0.05. You can also enter your own gene list, separated by line breaks (you can directly copy and paste multiple rows of gene names from the Excel file downloaded using \"Download statistical results\")."),
                     p("a. Click \"",
                       tags$span(style = "color:#4678B2; font-weight:bold", "fdr.overall < 0.05"),
                       "\" to use the default list of genes with an fdr.overall value < 0.05."),
                     p("b. Click \"",
                       tags$span(style = "color:#C95C54; font-weight:bold", "Clear All"),
                       "\" to clear all genes."),
                     hr(),
                     h4("2. Option for visualization"),
                     p("Input \"Number of pseudotime's bin\" for bin pseudotime t-test plot."),
                     p("Select \"Statistical significance\" to choose whether to perform FDR correction for bin pseudotime t-test plot."),
                     p("Click \"Find bins with max non-significant\" then select the search range and start searching to find bins with max non-siginifcant counts."),
                     p("Click \"Show individual data points\" to display cell points in the pseudotime variation plot of different Samples."),
                     p("Click \"Show error bar\" to show the error bar of pseudotime variation of different Samples plot."),
                     p("Click \"Include zero-expression cells\" to add cells with zero expression in the distribution plot of different Samples in the dimensionality reduction space."),
                     hr(),
                     h4("3. Start Analysis"),
                     p("Click \"",
                       tags$span(style = "color:#4678B2; font-weight:bold", "Start Analysis"),
                       "\" to start analysis. All the gene file will be save to the Result folder at visualization_YYYYmmddHHMM/HEATMAP/COMBINED, included gene name, information, population plot, bin pseudotime t-test plot, pseudotime variation of different Samples plot, and distribution of different Samples in the dimensionality reduction space plot. There will also be a \"gene_direction_summary.xlsx\" file which includes the count of bin pseudotime (Minus, NS, Plus, NA)")
                   ),
                   uiOutput("multi_gene_ui")
          ),
          
          # SINGLE GENE ANALYSIS Tab
          tabPanel("SINGLE GENE ANALYSIS", 
                   tags$div(
                     id = "single_gene_analysis_desc", class = "description-panel",
                     style = "display: none;",
                       h3("Single Gene Analysis Description"),
                       h4("Analyze specific genes.", style = "color: grey;"),
                       hr(),
                       h4("1. Select appropriate dimensions"),
                       p("Since harmony and PCA may have mixing of different biological categories under specific dimensions during visualization, before selecting genes for analysis, you can use \"X Axis\" and \"Y Axis\" to select dimensions where the biological categories in the Cell Cluster plot are separated and pseudotime is continuous."),
                     hr(),  
                     h4("2. Select genes"),
                       p("a. Click ",
                         tags$span(style = "color:#5CB85C; font-weight:bold", "\"Select Gene\""),
                         " to enter the gene selection page."
                       ),
                       p("b. In the selection page, you can sort ascending or descending and use the search function. Scroll horizontally to view all column names."),
                       p("c. Click the row of the gene you want to analyze, then click ",
                         tags$span(style = "color:#337AB7; font-weight:bold", "\"Confirm\""),
                         " to confirm the gene."
                       ),
                       h4("View cell cluster and pseudotime distributions"),
                       p("Besides the Cell Cluster plot and Pseudotime plot shown automatically at the start, you can also click ",
                         tags$span(style = "color:#929292; font-weight:bold", "\"Reset\""),
                         " after selecting genes to re-display these two plots."
                       ),
                     hr(),
                       h4("3. Different plots"),
                       p("When no gene is selected, two plots are displayed: the Cell Cluster plot and the Pseudotime plot, which are used to observe classifications of different biological cell types and pseudotime distributions."),
                       br(),
                       p("After gene selection, four plots are shown: population plot, bin pseudotime t-test plot, pseudotime variation of different Samples plot, and distribution of different Samples in the dimensionality reduction space plot."),
                       p("Input \"Number of pseudotime's bin\" for bin pseudotime t-test plot. The program applies a mixed model to each bin based on the number of bins to calculate the final p-value, followed by FDR correction. If the FDR is greater than 0.05, it is considered not significant. Otherwise, if the mean of test.variable=1 is greater than test.variable=0, it is marked as \"+\"; otherwise, it is marked as \"-\"."),
                       p("Select \"Statistical significance\" to choose whether to perform FDR correction for bin pseudotime t-test plot."),
                       p("Click \"Find bins with max non-significant\" then select the search range and start searching to find bins with max non-siginifcant counts."),
                       p("Click \"Show individual data points (Plot3)\" to display cell points in the pseudotime variation plot of different Samples."),
                       p("Click \"Show error bar (Plot3)\" to show the error bar of pseudotime variation of different Samples plot."),
                       p("Click \"Include zero-expression cells (Plot4)\" to add cells with zero expression in the distribution plot of different Samples in the dimensionality reduction space."),
                       p("Uncheck \"Sample included\" to remove the corresponding sample from the pseudotime variation plot, facilitating observation of relationships among other samples."),
                     hr(),  
                     h4("4. Save images"),
                       p("At any time (whether genes are selected or not), you can click ",
                         tags$span(style = "color:#337AB7; font-weight:bold", "\"Save Plot\""),
                         " to save displayed plots to the Result folder at visualization_YYYYmmddHHMM/HEATMAP. If no gene is selected, the two plots will be saved in the reduction folder; otherwise, they will be saved in the folder corresponding to the gene name. Meanwhile, a COMBINED folder will be generated including gene information and all four plots."
                       )
                     ),
                   uiOutput("gene_ui")
          ),
          
          
          # DOWNLOAD VISUALIZATION Tab
          tabPanel("DOWNLOAD VISUALIZATION", 
                   tags$div(
                     id = "download_visualization_desc", class = "description-panel",
                     style = "display: none;",
                     h3("Download Visualization Description"),
                     h4("Download all visualization results", style = "color: grey;"),
                     hr(),
                     p("Download the results of visualization_YYYYmmddHHMM in the result folder to the computer's default download folder.")
                   ),
                   uiOutput("download_tab_content")
          )
        )
      )
    }
  })
  
  observeEvent(input$navlist, {
  shinyjs::hide(".description-panel")
  
  desc_id <- switch(input$navlist,
    "HEATMAP" = "heatmap_desc",
    "MULTI GENE ANALYSIS" = "multi_gene_analysis_desc",
    "SINGLE GENE ANALYSIS" = "single_gene_analysis_desc",
    "DOWNLOAD VISUALIZATION" = "download_visualization_desc"
  )
  shinyjs::show(desc_id)
}, ignoreInit = FALSE)
  
  output$download_stat <- downloadHandler(
    filename = function() {"Lamian_statistical_result.xlsx"},
    content = function(file) {
      req(input$result_folder)
      source_file <- file.path(input$result_folder, "stat.xlsx")
      if (file.exists(source_file)) {
        file.copy(source_file, file)
        showNotification(
          "Download successful! Please check 'Lamian_statistical_result.xlsx' in the default download folder",
          type = "message",
          duration = 10
        )
      } else {
        showNotification("File not found!", type = "error")
      }
    }
  )
  
  parameters_content <- reactive({
    req(input$result_folder)
    file_path <- file.path(input$result_folder, "parameters.txt")
    if(file.exists(file_path)) {
      readLines(file_path)
    } else {
      NULL
    }
  })
  
  output$parameterDisplay <- renderUI({
    content <- parameters_content()
    if(is.null(content)) return(p("No parameters file found."))
    
    tags$div(
      id = "parameter_desc", 
      class = "description-panel",
      style = "display: block;",
      h4("Click the tab on the left to start data visualization analysis", style = "color: red;"),
      h3("ANALYSIS PARAMETERS", style = "color: #2c3e50;"),
      h4("Detailed configuration of the Lamian analysis", style = "color: #7f8c8d;"),
      hr(style = "border-color: #bdc3c7;"),
      
      tags$pre(
        style = "
        background-color: #f8f9fa;
        border: 1px solid #e9ecef;
        border-radius: 5px;
        padding: 20px;
        font-family: 'Courier New', monospace;
        font-size: 14px;  
        line-height: 1.6; 
        overflow-x: auto;
        white-space: pre-wrap;
        word-wrap: break-word;
        color: #2c3e50;
        margin: 0;
      ",
        paste(content, collapse = "\n")
      ),
      
      hr(style = "border-color: #bdc3c7; margin: 20px 0;"),
      p(
        style = "font-size: 14px; color: #7f8c8d;",
        icon("info-circle"), 
        " This panel displays the complete parameter configuration used in the Lamian analysis."
      )
    )
  })
  
  output$heatmap_ui <- renderUI({
    sig_genes <- rownames(stat()[stat()$fdr.overall < 0.05, ])
    fluidRow(
      column(
        width = 5,
        wellPanel(
          style = "height: 90vh; overflow-y: auto;",
          h4("Gene List Input", class = "text-primary"),
          hr(),
          
          textAreaInput(
            inputId = "heatmap_gene_list",
            label = "Enter genes (one per line):",
            value = if(length(sig_genes) > 0) paste(sig_genes, collapse = "\n") else "",
            rows = 10,
            placeholder = "Paste genes here:\nGene_A\nGene_B\nGene_C"
          ),
          actionButton("load_default", "fdr.overall < 0.05", 
                                   class = "btn-primary", width = "100%"),
          actionButton("clear_genes", "Clear All", 
                                   class = "btn-danger", width = "100%"),
          
          hr(),
          h4("Heatmap Parameters", class = "text-info"),
          
          radioButtons(
            "cluster_method",
            "Clustering Method:",
            choices = c("Hierarchical" = "hierarchical", 
                        "K-means" = "kmeans"),
            selected = "hierarchical",
            inline = TRUE
          ),
          
          radioButtons(
            "scale_by",
            "Scale By:",
            choices = c("Both groups together" = "both", 
                        "Each group separately" = "separately"),
            selected = "both",
            inline = FALSE
          ),
          
          radioButtons(
            "cluster_by",
            "Cluster By:",
            choices = c("Both groups together" = "both", 
                        "Each group separately" = "separately"),
            selected = "both",
            inline = FALSE
          ),
          
          fluidRow(
            column(4, colourpicker::colourInput("low_color", "Low:", value = "blue", 
                                                showColour = "background", width = "100%")),
            column(4, colourpicker::colourInput("mid_color", "Mid:", value = "white", 
                                                showColour = "background", width = "100%")),
            column(4, colourpicker::colourInput("high_color", "High:", value = "red", 
                                                showColour = "background", width = "100%"))
          ),
          
          uiOutput("cluster_num_ui"),
          
          actionButton("draw_plot", "Draw Plot", 
                       class = "btn-success", width = "100%"),
          br(),
          uiOutput("gene_tools_ui"),

          actionButton("save_heatmap_plot", "Save Plot", 
                       class = "btn-primary", width = "100%", disabled = TRUE),
          
          hr(),
          
          numericInput(
            "pseudotime_samples",
            "Number of Pseudotime samples:",
            value = 11,
            min = 3,
            max = 100,
            step = 1
          ),
         
          actionButton("save_heatmap_data", "Save Data", 
                       class = "btn-info", width = "100%", disabled = TRUE)
        )
      ),
      
      column(
        width = 7,
        h4("Heatmap", class = "text-primary"),
        hr(),
        plotOutput("heatmap_display")  
      )
    )
  })
  
  output$gene_ui <- renderUI({
    
    fluidRow(
      column(
        width = 2,
        style = "border-right: 1px solid #eee; padding-right: 15px;",
        if(!is.na(gene_analysis_selected())) {
          actionButton("original_btn", "Reset", 
                       class = "btn-default", 
                       style = "width: 100%; margin-bottom: 15px;")
        },
        actionButton("select_gene", "Select Gene", 
                     class = "btn-success", 
                     style = "width: 100%; margin-bottom: 15px;"),
        actionButton("save_plot", "Save Plot", 
                     class = "btn-primary", 
                     style = "width: 100%; margin-bottom: 15px;"),
        uiOutput("selected_gene_info")
      ),
      column(
        width = 10,
        uiOutput("gene_plot_ui", height = "800px")
      )
    )
  })
  
  output$multi_gene_ui <- renderUI({
    sig_genes <- rownames(stat()[stat()$fdr.overall < 0.05, ])
    fluidRow(
      column(
        width = 3,
        wellPanel(
          style = "height: 90vh; overflow-y: auto;",
          h4("Gene List Input", class = "text-primary"),
          hr(),
          textAreaInput(
            inputId = "multi_gene_list",
            label = "Enter genes (one per line):",
            value = if(length(sig_genes) > 0) paste(sig_genes, collapse = "\n") else "",
            rows = 10,
            placeholder = "Paste genes here:\nGene_A\nGene_B\nGene_C"
          ),
          actionButton("load_default_multi", "fdr.overall < 0.05", 
                       class = "btn-primary", width = "100%"),
          actionButton("clear_genes_multi", "Clear All", 
                       class = "btn-danger", width = "100%"),
          hr(),
          numericInput("multi_bin_number", "Number of Bin Input for the Analysis:", min = 3, max = 30, value = 10),
          radioButtons(
            "multi_using_p_adj",
            "Statistical significance:",
            choices = c(
              "Raw p-value" = "pvalue",
              "FDR adjusted" = "p_adjust"
            ),
            selected = "p_adjust",
            inline = TRUE
          ),
          checkboxInput("show_points_multi", "Show individual data points", 
                        value = FALSE, width = "100%"),
          checkboxInput("show_error_bar_multi", "Show error bar", 
                        value = FALSE, width = "100%"),
          conditionalPanel(
            condition = "input.show_error_bar_multi == true",
            sliderInput("confidence_level_multi", "Confidence Level:",
                        min = 90, max = 99.9, value = 95, step = 0.1,
                        post = "%", width = "100%")
          ),
          checkboxInput("included_zero_multi", "Include zero-expression cells", 
                        value = FALSE, width = "100%"),
          hr(),
          actionButton("save_multi_gene", "Start Analysis", 
                       class = "btn btn-success", width = "100%")
        )
      ),
      column(
        width = 9,
        wellPanel(
          style = "height: 90vh; overflow-y: auto;",
          h4("Analysis Progress", class = "text-primary"),
          hr(),
          p("Click 'Start Analysis' to begin processing. A progress window will appear.")
        )
      )
    )
  })
  
  output$download_tab_content <- renderUI({
    fluidRow(
      column(12,
             div(
               style = "color: red; font-weight: bold; margin-bottom: 15px;",
               icon("exclamation-triangle"), 
               "If you have not generated any visualization files yet. Please use \"HEATMAP\", \"SINGLE GENE ANALYSIS\" or \"MULTI GENE ANALYSIS\" to generate them first."
             ),
             downloadButton(
               "download_visualization", 
               label = "Download All Visualization Files",
               style = "width: 100%;"
             )
      )
    )
  })
  
  #HEATMAP UI
  observeEvent(input$load_default, {
    sig_genes <- rownames(stat()[stat()$fdr.overall < 0.05, ])
    updateTextAreaInput(
      session, 
      "heatmap_gene_list",
      value = if(length(sig_genes) == 0){"NO SIGINIFICANT GENE"}else{paste(sig_genes, collapse = "\n")}
    )
  })
  
  observeEvent(input$clear_genes, {
    updateTextAreaInput(session, "heatmap_gene_list", value = "")
  })
  
  output$cluster_num_ui <- renderUI({
    req(input$cluster_method)
    
    default_val <- ifelse(input$cluster_method == "kmeans", 6, 1)
    
    numericInput(
      "cluster_number",
      "Number of Clusters:",
      value = default_val,
      min = 1,
      max = 20,
      step = 1
    )
  })
  
  observeEvent(input$cluster_method, {
    default_val <- ifelse(input$cluster_method == "kmeans", 6, 1)
    updateNumericInput(session, "cluster_number", value = default_val)
  })
  
  observeEvent(input$draw_plot, {
    req(xde_result(), input$heatmap_gene_list)
    
    showModal(modalDialog(
      title = "Processing",
      h4("Generating heatmap...", style = "text-align: center;"),
      footer = NULL,
      size = "s",
      easyClose = FALSE
    ))
    
    genes <- unlist(strsplit(input$heatmap_gene_list, "\n")) %>% 
      trimws() %>% 
      .[. != ""]
    
    if(length(genes) == 0) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        "Gene list cannot be empty",
        footer = modalButton("Close"),
        easyClose = FALSE
      ))
      return()
    }
    
    result <- tryCatch({
      draw_heatmap(
        xde_result = xde_result(),
        gene_list = genes,
        cluster_method = input$cluster_method,
        scale_by = input$scale_by,
        cluster_by = input$cluster_by,
        low_color = input$low_color,
        mid_color = input$mid_color,
        high_color = input$high_color,
        cluster_number = input$cluster_number
      )
    }, error = function(e) {
      list(
        plot = NULL,
        df = NULL,
        warning_message = paste("Error:", e$message)
      )
    })
    
    heatmap_plot(result$plot)
    heatmap_data(result$df)
    heatmap_name(result$names)
  
    if(!is.null(result$plot)) {
      output$heatmap_display <- renderPlot({
        result$plot
      })
      removeModal()
      if(!is.na(result$warning_message)) {
        showModal(modalDialog(
          title = "Warning",
          result$warning_message,
          footer = modalButton("Close"),
          easyClose = FALSE
        ))
      }
    } else {
      showModal(modalDialog(
        title = "Error",
        result$warning_message,
        footer = modalButton("Close"),
        easyClose = FALSE
      ))
    }
    
    output$gene_tools_ui <- renderUI({
      req(input$cluster_number)  
      tagList(
        fluidRow(
          column(
            width = 3,
            tags$span("Cluster:", style = "line-height: 35px; font-weight: bold;")
          ),
          column(
            width = 9,
            selectInput("cluster_select_copy", NULL,
                        choices = seq_len(input$cluster_number),
                        width = "100%",
                        selected = 1)
          )
        ),
        fluidRow(
          column(
            width = 12,
            actionButton("copy_gene", "Copy Gene to Clipboard", 
                         class = "btn-warning", width = "100%")
          )
        )
      )
    })
  })
  
  observeEvent(input$copy_gene, {
    df <- heatmap_data() 
    req(df, input$cluster_select_copy)
    
    selected_genes <- df$GENE[df$Cluster_Number == input$cluster_select_copy]
    
    if (length(selected_genes) > 0) {
      gene_text <- paste(selected_genes, collapse = "\n")
      runjs(sprintf("navigator.clipboard.writeText(`%s`);", gene_text))
      showNotification(paste("Copied", length(selected_genes), "genes to clipboard."),
                       type = "message")
    } else {
      showNotification("No genes found for this cluster.", type = "warning")
    }
  })
  
  observeEvent(input$save_heatmap_plot, {
    req(heatmap_plot(), heatmap_data(), input$result_folder)
    if (is.null(result_visulization_folder())) {
      time_str <- format(Sys.time(), "%Y%m%d%H%M")
      vis_dir <- file.path(input$result_folder, paste0("visualization_", time_str))
      if (!dir.exists(vis_dir)) dir.create(vis_dir, recursive = TRUE)
      result_visulization_folder(vis_dir)
    }
    heatmap_folder <- file.path(result_visulization_folder(),"HEATMAP")
    if (!dir.exists(heatmap_folder)) dir.create(heatmap_folder, recursive = TRUE)
    
    col_fun <- circlize::colorRamp2(
      c(-2, 0, 2), 
      c(input$low_color, input$mid_color, input$high_color)
    )
    
    params_title <- paste(
      "Cluster Method:", input$cluster_method,
      "| Scale By:", input$scale_by,
      "| Cluster By:", input$cluster_by
    )
    
    tryCatch({
      n_rows <- nrow(heatmap_data())
      plot_height <- max(8, n_rows * 0.03 + 6)  
      
      pdf(file.path(heatmap_folder, paste0("WithoutGeneName_", input$cluster_method, ".", input$cluster_number, "_", input$scale_by,".", input$cluster_by, ".pdf")), 
          width = 10, height = plot_height) 
      grid.text(params_title, 
                x = 0.5, y = 0.98, 
                just = "top",
                gp = gpar(fontsize = 12, fontface = "bold"))
      
      pushViewport(viewport(y = 0.95, height = 0.9, just = c("center","top")))
      ComplexHeatmap::draw(heatmap_plot(), 
                           newpage = FALSE,
                           padding = unit(c(2, 2, 2, 2), "mm"))
      popViewport()
      dev.off()
    }, error = function(e) {
      showNotification(paste("Failed to save plot without gene names:", e$message), type = "error")
    })
    
    tryCatch({
      df <- heatmap_data()
      n_rows <- nrow(df)
      
      scale0_mat <- as.matrix(df[, grep("^scale0_", colnames(df), value = TRUE)])
      scale1_mat <- as.matrix(df[, grep("^scale1_", colnames(df), value = TRUE)])
      rownames(scale0_mat) <- rownames(scale1_mat) <- df$GENE
      
      base_height <- 8  
      row_height <- 0.15  
      gene_name_width <- max(nchar(df$GENE)) * 0.15 + 1  
      
      plot_height <- max(base_height, n_rows * row_height + 6)
      plot_width <- 10 + gene_name_width * 0.5  
      
      cluster_colors <- setNames(
        color_selected(length(unique(df$Cluster_Number))),
        sort(unique(df$Cluster_Number))
      )
      
      ht0 <- Heatmap(
        scale0_mat,
        name = "expression",
        column_title = heatmap_name()[1],
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = col_fun,
        left_annotation = rowAnnotation(
          Cluster = as.factor(df$Cluster_Number),
          col = list(Cluster = cluster_colors),
          show_annotation_name = FALSE
        ),
        heatmap_legend_param = list(title = "Gene expression level")
      )
      
      ht1 <- Heatmap(
        scale1_mat,
        name = "expression", 
        column_title = heatmap_name()[2],
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 6),  
        row_names_max_width = unit(gene_name_width, "cm"),  
        show_column_names = FALSE,
        col = col_fun,
        heatmap_legend_param = list(title = NULL)
      )
      
      pdf(file.path(heatmap_folder, paste0("WithGeneName_", input$cluster_method,".",input$cluster_number, "_", input$scale_by, ".", input$cluster_by, ".pdf")), 
          width = plot_width, height = plot_height)
      grid.text(params_title, 
                x = 0.5, y = 0.98, 
                just = "top",
                gp = gpar(fontsize = 12, fontface = "bold"))
      
      pushViewport(viewport(y = 0.95, height = 0.9, just = c("center","top")))
      ComplexHeatmap::draw(ht0 + ht1,
                           heatmap_legend_side = "bottom",
                           merge_legend = TRUE,
                           newpage = FALSE,
                           padding = unit(c(2, 2, 2, 2), "mm"))
      popViewport()
      dev.off()
      
      showNotification("Heatmaps saved successfully!", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Failed to save plot with gene names:", e$message), type = "error")
    })
  })
  
  
  observeEvent(input$save_heatmap_data, {
    req(heatmap_data(), input$result_folder, input$pseudotime_samples)
    
    tryCatch({
      df <- heatmap_data()
      calculate_sample_points <- function(total_points, sample_num) {
        if(sample_num == 1) return(1) 
        interval <- total_points / (sample_num - 1)
        
        points <- round(seq(interval, total_points-1, by = interval))
        points <- unique(c(1, points, total_points)) 
        points <- sort(points)[1:sample_num]  
        
        return(points)
      }
      sample_points <- calculate_sample_points(total_points = 100, sample_num = input$pseudotime_samples)
      
      extract_samples <- function(col_prefix) {
        cols <- grep(paste0("^", col_prefix), colnames(df), value = TRUE)
        if(length(cols) > 0) {
          mat <- as.matrix(df[, cols, drop = FALSE])
          sampled <- mat[, sample_points, drop = FALSE]
          colnames(sampled) <- paste0(col_prefix, sample_points)
          return(sampled)
        }
        return(NULL)
      }
      
      sampled_data <- cbind(
        df[, c("GENE", "Cluster_Number"), drop = FALSE],
        extract_samples("scale0_"),
        extract_samples("scale1_"),
        extract_samples("raw0_"),
        extract_samples("raw1_") 

      )
      
      if (is.null(result_visulization_folder())) {
        time_str <- format(Sys.time(), "%Y%m%d%H%M")
        vis_dir <- file.path(input$result_folder, paste0("visualization_", time_str))
        if (!dir.exists(vis_dir)) dir.create(vis_dir, recursive = TRUE)
        result_visulization_folder(vis_dir)
      }
      heatmap_folder <- file.path(result_visulization_folder(),"HEATMAP")
      if (!dir.exists(heatmap_folder)) dir.create(heatmap_folder, recursive = TRUE)
      
      
      writexl::write_xlsx(
        list(HeatmapData = sampled_data),
        path = file.path(heatmap_folder, paste0("Data_", input$cluster_method,".",input$cluster_number, "_", input$scale_by, ".", input$cluster_by, ".xlsx"))
      )
      
      showNotification(
        paste("Sampled data saved with", length(sample_points), "time points"),
        type = "message", 
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(
        paste("Failed to save sampled data:", e$message),
        type = "error"
      )
    })
  })
  

  #GENE UI
  observeEvent(input$select_gene, {
    req(xde_result())
    showModal(modalDialog(
      title = "Select Gene",
      DTOutput("gene_selection_table"),
      footer = tagList(
        actionButton("confirm_gene_selection", "Confirm", class = "btn-primary"),
        modalButton("Cancel")
      ),
      size = "l"
    ))
  })
  
  output$gene_selection_table <- renderDT({
    req(xde_result())
    stats_df <- as.data.frame(xde_result()$statistics) 
    stats_df <- cbind(Gene = rownames(stats_df), stats_df)
    rownames(stats_df) <- NULL
    
    datatable(
      stats_df,
      selection = 'single',
      rownames = FALSE,
      filter = 'top',
      options = list(
        scrollX = TRUE,  
        autoWidth = TRUE, 
        columnDefs = list(
          list(width = '100px', targets = "_all")  
        ),
        pageLength = 10,
        searchHighlight = TRUE
      )
    ) %>%
      formatSignif(columns = sapply(stats_df, is.numeric), digits = 3) %>%
      formatStyle('Gene', color = 'red')
  })
  
  observeEvent(input$confirm_gene_selection, {
    selected_row <- input$gene_selection_table_rows_selected
    if(length(selected_row)) {
      stats_df <- as.data.frame(xde_result()$statistics)
      stats_df <- cbind(Gene = rownames(stats_df), stats_df)
      gene_analysis_selected(stats_df[selected_row, "Gene"])
      removeModal()
    } else {
      showNotification("Please select a gene first!", type = "warning")
    }
  })
  
  observeEvent(input$original_btn, {
    gene_analysis_selected(NA)
  })
  
  output$selected_gene_info <- renderUI({
    req(gene_analysis_selected())
    if(!is.na(gene_analysis_selected())) {
      gene_stats <- xde_result()$statistics[gene_analysis_selected(), , drop = FALSE]
      tagList(
        h4(paste0("Selected Gene: ", gene_analysis_selected()), style = "color: red;"),
        lapply(colnames(gene_stats), function(col_name) {
          tagList(
            p(strong(col_name, ":")),
            p(style = "color: red; margin-left: 20px;", 
              formatC(gene_stats[, col_name], format = "e", digits = 2)),
            hr(style = "margin: 5px 0;")
          )
        })
      )
    }
  })
  
  output$gene_plot_ui <- renderUI({
    req(sce(), xde_result())
    
    if (is.na(gene_analysis_selected())) {
      dim_names <- colnames(reducedDim(sce(), "Reduction"))
      
      tagList(
        fluidRow(
          column(6, plotOutput("gene_plot1_NA", height = "400px")),
          column(6, plotOutput("gene_plot2_NA", height = "400px"))
        )
      )
    } else {
      tagList(
        wellPanel(
          fluidRow(
            column(
              6, 
              numericInput("bin_number", "Number of Bin Input for the Analysis (Plot2):", min = 3, max = 50, value = 10)
            ),
            column(
              6,
              radioButtons(
                "using_p_adj",
                "Statistical significance (Plot2):",
                choices = c(
                  "Raw p-value" = "pvalue",
                  "FDR adjusted" = "p_adjust"
                ),
                selected = "p_adjust",
                inline = TRUE
              ))
          ),
          checkboxInput(
            "search_nonsig_bins",
            tags$span(icon("search"), "Find bins with max non-significant counts"),
            value = FALSE
          ),
          
          fluidRow(
            column(9,
                   conditionalPanel(
                     condition = "input.search_nonsig_bins == true",
                     sliderInput(
                       "search_bin_range",
                       "",
                       min = 5,
                       max = 50,
                       value = c(5, 30),
                       step = 1,
                       width = "100%",
                       ticks = TRUE
                     )
                   )
            ),
            column(3,
                   conditionalPanel(
                     condition = "input.search_nonsig_bins == true",
                     div(style = "display: flex; align-items: center; height: 100%; min-height: 60px;",
                         actionButton("search_bin", "Search", 
                                      icon = icon("search"),
                                      class = "btn-primary",
                                      width = "100%")
                     )
                   )
            )
          ),
          splitLayout(
            cellWidths = c("50%", "50%"),
            checkboxInput("show_points", "Show individual data points (Plot3)", 
                          value = FALSE, width = "95%"),
            checkboxInput("show_error_bar", "Show error bar (Plot3)", 
                          value = FALSE, width = "95%")
          ),
          
          conditionalPanel(
            condition = "input.show_error_bar == true",
            fluidRow(
              column(
                12,
                sliderInput("confidence_level", "Confidence Level:",
                            min = 90, max = 99.9, value = 95, step = 0.1,
                            post = "%", width = "100%")
              )
            )
          ),
          fluidRow(
            column(
              6,
              checkboxInput("included_zero", "Include zero-expression cells (Plot4)", 
                            value = FALSE, width = "100%")
            )
          )
        ),  
        
        h4("Sample included"),
        uiOutput("sample_selector_ui"),
        
        fluidRow(
          column(4, plotOutput("gene_plot1", height = "400px")),
          column(4, plotOutput("gene_plot2", height = "400px")),
          column(4, plotOutput("gene_plot3", height = "360px"))
        ),
        fluidRow(
          column(12, plotOutput("gene_plot4", height = "400px"))
        )
      )
    }
  })
  
  observeEvent(input$search_bin, {
    req(input$search_nonsig_bins, input$search_bin_range, gene_analysis_selected(), xde_result())
    
    showModal(modalDialog(
      title = "Calculating",
      tags$h4("Searching for bin number with maximum significant points..."),
      tags$br(),
      progressBar(
        id = "search_progress",
        value = 0,
        title = "Progress:",
        display_pct = TRUE
      ),
      footer = NULL,
      easyClose = FALSE,
      size = "m"
    ))
    
    bin_range <- input$search_bin_range
    min_bins <- bin_range[1]
    max_bins <- bin_range[2]
    total_bins <- max_bins - min_bins + 1
    gene_name <- gene_analysis_selected()
    
    search_results <- data.frame(
      bin_number = integer(),
      significant_count = integer(),
      stringsAsFactors = FALSE
    )
    
    gene_expr <- xde_result()$expr[gene_name, ]
    for (n_bins in min_bins:max_bins) {
      
      up_value = round(((n_bins+1-min_bins) / total_bins) * 100)
      
      pseudotime_bins <- cut(xde_result()$pseudotime, breaks = n_bins, labels = FALSE)
      analysis_data <- data.frame(
        expr = as.numeric(gene_expr),
        cell_id = names(gene_expr),
        pseudotime_bin = pseudotime_bins,
        pseudotime = xde_result()$pseudotime  
      )
      
      analysis_data <- merge(analysis_data, xde_result()$cellanno, by.x = "cell_id", by.y = "Cell")
      analysis_data <- merge(analysis_data, xde_result()$design, by.x = "Sample", by.y = "row.names")
      
      results <- data.frame(
        bin = 1:n_bins,
        mean_0 = numeric(n_bins),
        mean_1 = numeric(n_bins),
        p_value = numeric(n_bins),
        stringsAsFactors = FALSE
      )
      
      for (i in 1:n_bins) {
        bin_data <- analysis_data[analysis_data$pseudotime_bin == i, ]
        
        intercept_index <- which(colnames(bin_data) == "intercept")
        treatment_col <- colnames(bin_data)[intercept_index + 1]
        
        confounder_cols <- colnames(bin_data)[(intercept_index + 2):ncol(bin_data)]
        
        if (length(confounder_cols) > 0) {
          model_formula <- as.formula(paste("expr ~", treatment_col, "+", 
                                            paste(confounder_cols, collapse = " + "),
                                            "+ (1 | Sample)"))
        } else {
          model_formula <- as.formula(paste("expr ~", treatment_col, "+ (1 | Sample)"))
        }
        
        if (nrow(bin_data) > n_bins) {  
          tryCatch({
            suppressMessages(suppressWarnings({
              model <- lmer(model_formula, data = bin_data)
            }))
            
            model_summary <- summary(model)
            treatment_coef <- which(rownames(model_summary$coefficients) == treatment_col)
            
            if (length(treatment_coef) > 0) {
              results$p_value[i] <- model_summary$coefficients[treatment_coef, "Pr(>|t|)"]
            } else {
              results$p_value[i] <- NA
            }
            
          }, error = function(e) {
            results$p_value[i] <- NA
          })
        } else {
          results$p_value[i] <- NA
        }
        
        expr_0 <- bin_data$expr[bin_data[[treatment_col]] == 0]
        expr_1 <- bin_data$expr[bin_data[[treatment_col]] == 1]
        results$mean_0[i] <- mean(expr_0, na.rm = TRUE)
        results$mean_1[i] <- mean(expr_1, na.rm = TRUE)
      }
      
      if(input$using_p_adj == "pvalue"){
        results$significant <- ifelse(is.na(results$p_value), NA, 
                                      ifelse(results$p_value < 0.05, TRUE, FALSE))
      } else {
        valid_p_indices <- which(!is.na(results$p_value))
        if (length(valid_p_indices) > 0) {
          results$p_adjusted <- NA
          results$p_adjusted[valid_p_indices] <- p.adjust(results$p_value[valid_p_indices], method = "fdr")
        }
        results$significant <- ifelse(is.na(results$p_adjusted), NA,
                                      ifelse(results$p_adjusted < 0.05, TRUE, FALSE))
      }
      
      sig_count <- sum(results$significant == TRUE, na.rm = TRUE)
      
      search_results <- rbind(search_results, data.frame(
        bin_number = n_bins,
        significant_count = sig_count
      ))
      updateProgressBar(
        session = session,
        id = "search_progress",
        value = up_value,
        title = paste("Processing bin", n_bins, "of", max_bins)
      )
      Sys.sleep(0.01)
    }
    
    max_sig_count <- max(search_results$significant_count, na.rm = TRUE)
    max_sig_bins <- search_results$bin_number[search_results$significant_count == max_sig_count]
    
    optimal_bins <- min(max_sig_bins)
    
    removeModal()
    
    showModal(modalDialog(
      title = "Search Results",
      tags$h4(paste("Optimal bin number:", optimal_bins)),
      tags$h5(paste("Maximum significant bins:", max_sig_count)),
      if (length(max_sig_bins) > 1) {
        tags$p(paste("Other bin numbers with same significance:", 
                     paste(setdiff(max_sig_bins, optimal_bins), collapse = ", ")),
               style = "color: blue;")
      },
      tags$hr(),
      tags$h4("All Results (Grouped by significance count):"),
      
      renderTable({
        grouped_results <- split(search_results$bin_number, search_results$significant_count)
        
        result_table <- data.frame(
          `Total Significant Bins` = as.numeric(names(grouped_results)),
          `Number of Bin Input for the Analysis` = sapply(grouped_results, function(x) {
            paste(sort(x), collapse = ", ") 
          }),
          stringsAsFactors = FALSE,
          check.names = FALSE  
        )
        
        result_table <- result_table[order(-result_table$`Total Significant Bins`), ]
        
        result_table
      }, digits = 0, na = "", align = 'c', width = "100%"),
      footer = modalButton("Close"),
      size = "l",
      easyClose = TRUE
    ))
    
    updateNumericInput(session, "bin_number", value = optimal_bins)
  })
  
  observeEvent(input$show_points, {
    show_points(input$show_points)
  })
  
  observeEvent(input$included_zero, {
    included_zero(input$included_zero)
  })
  
  coloredCheckboxGroup <- function(inputId, label, choices, selected, colors) {
    choices <- as.list(choices)
    names(choices) <- NULL
    
    tags$div(
      id = inputId,
      class = "form-group shiny-input-checkboxgroup shiny-input-container",
      tags$label(class = "control-label", `for` = inputId, label),
      tags$div(
        class = "shiny-options-group",
        lapply(seq_along(choices), function(i) {
          choice <- choices[[i]]
          choice_value <- choice
          choice_name <- choice
          
          tags$div(
            class = "checkbox",
            tags$label(
              style = paste0("color:", colors[choice], "; font-weight: bold;"),
              tags$input(
                type = "checkbox",
                name = inputId,
                value = choice_value,
                checked = if (choice_value %in% selected) "checked"
              ),
              choice_name
            )
          )
        })
      )
    )
  }
  
  output$sample_selector_ui <- renderUI({
    req(xde_result())
    samples <- unique(xde_result()$cellanno$Sample)
    
    colors_all <- color_selected(length(samples))
    names(colors_all) <- samples
    
    half_len <- ceiling(length(samples) / 2)
    samples_left <- samples[1:half_len]
    samples_right <- samples[(half_len + 1):length(samples)]
    
    fluidRow(
      column(6,
             coloredCheckboxGroup(
               inputId = "selected_samples_left",
               label = NULL,
               choices = samples_left,
               selected = samples_left,
               colors = colors_all
             )
      ),
      column(6,
             coloredCheckboxGroup(
               inputId = "selected_samples_right",
               label = NULL,
               choices = samples_right,
               selected = samples_right,
               colors = colors_all
             )
      )
    )
  })
  
  output$gene_plot1_NA <- renderPlot({
    req(sce())
    rd <- reducedDim(sce(), "UMAP")
    
    p <- plotReducedDim(sce(), 
                        dimred = "UMAP", 
                        colour_by = "cluster", 
                        point_alpha = 1,
                        point_size = 0.5) +
      labs(
        x = colnames(rd)[1],  
        y = colnames(rd)[2],  
        title = "Cell Clusters"
      ) +
      theme_light() +
      theme(
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.5, "cm"),  
        axis.title = element_text(size = 16)
      )+
      guides(colour = guide_legend(
        override.aes = list(size = 5)   
      ))
    
    gene_plot1_NA(p)
    return(p)
  })
  
  output$gene_plot2_NA <- renderPlot({
    req(sce(), xde_result())
    
    rd <- reducedDim(sce(), "UMAP")
    pt <- xde_result()$pseudotime
    rd <- rd[rownames(rd) %in% names(pt), ]
    
    p <- ggplot(data.frame(rd, pseudotime = pt), 
                aes(x = !!sym(colnames(rd)[1]), 
                    y = !!sym(colnames(rd)[2]), 
                    color = pseudotime)) +
      geom_point(size = 0.5) +
      scale_color_viridis_c() +
      ggtitle("Pseudotime") +
      theme_light()+
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            axis.title = element_text(size = 16))
    
    gene_plot2_NA(p)
    return(p)
  })
  
  output$gene_plot1 <- renderPlot({
    req(gene_analysis_selected())
    p <- plotGenePopulation(xde_result(), gene_analysis_selected(), type = "variable") + theme(
      legend.title = element_text(size = 14),     
      legend.text = element_text(size = 10),      
      axis.title.x = element_text(size = 14),      
      axis.title.y = element_text(size = 14),     
      axis.text.x = element_text(size = 12),       
      axis.text.y = element_text(size = 12),
      legend.position = "bottom"
    )
    gene_plot1(p)
    return(p)
  })
  
  output$gene_plot2 <- renderPlot({
    req(gene_analysis_selected(), sce(), xde_result())
    gene_name <- gene_analysis_selected()
    bin_number_debounced <- reactive({
      input$bin_number
    }) %>% debounce(500)
    n_bins <- bin_number_debounced()
    gene_expr <- xde_result()$expr[gene_name, ]
    pseudotime_bins <- cut(xde_result()$pseudotime, breaks = n_bins, labels = FALSE)
    
    analysis_data <- data.frame(
      expr = as.numeric(gene_expr),
      cell_id = names(gene_expr),
      pseudotime_bin = pseudotime_bins,
      pseudotime = xde_result()$pseudotime  
    )
    
    analysis_data <- merge(analysis_data, xde_result()$cellanno, by.x = "cell_id", by.y = "Cell")
    analysis_data <- merge(analysis_data, xde_result()$design, by.x = "Sample", by.y = "row.names")
    results <- data.frame(
      bin = 1:n_bins,
      mean_0 = numeric(n_bins),
      mean_1 = numeric(n_bins),
      p_value = numeric(n_bins),
      significance = character(n_bins),
      stringsAsFactors = FALSE
    )
    
    
    for (i in 1:n_bins) {
      bin_data <- analysis_data[analysis_data$pseudotime_bin == i, ]
      
      intercept_index <- which(colnames(bin_data) == "intercept")
      treatment_col <- colnames(bin_data)[intercept_index + 1]
      
      confounder_cols <- colnames(bin_data)[(intercept_index + 2):ncol(bin_data)]
      
      if (length(confounder_cols) > 0) {
        model_formula <- as.formula(paste("expr ~", treatment_col, "+", 
                                          paste(confounder_cols, collapse = " + "),
                                          "+ (1 | Sample)"))
      } else {
        model_formula <- as.formula(paste("expr ~", treatment_col, "+ (1 | Sample)"))
      }
      
      if (nrow(bin_data) > n_bins) {  
        tryCatch({
          suppressMessages(suppressWarnings({
            model <- lmer(model_formula, data = bin_data)
          }))
          
          model_summary <- summary(model)
          treatment_coef <- which(rownames(model_summary$coefficients) == treatment_col)
          
          if (length(treatment_coef) > 0) {
            results$p_value[i] <- model_summary$coefficients[treatment_coef, "Pr(>|t|)"]
          } else {
            results$p_value[i] <- NA
          }
          
        }, error = function(e) {
          results$p_value[i] <- NA
        })
      } else {
        results$p_value[i] <- NA
      }
      
      expr_0 <- bin_data$expr[bin_data[[treatment_col]] == 0]
      expr_1 <- bin_data$expr[bin_data[[treatment_col]] == 1]
      results$mean_0[i] <- mean(expr_0, na.rm = TRUE)
      results$mean_1[i] <- mean(expr_1, na.rm = TRUE)
    }
    
    valid_p_indices <- which(!is.na(results$p_value))
    
    if(input$using_p_adj == "pvalue"){
      for (i in 1:n_bins) {
        if (is.na(results$p_value[i])) {
          results$significance[i] <- "NA"
        } else {
          if (results$p_value[i] < 0.001) {
            results$significance[i] <- "***"
          } else if (results$p_value[i] < 0.01) {
            results$significance[i] <- "**"
          } else if (results$p_value[i] < 0.05) {
            results$significance[i] <- "*"
          } else {
            results$significance[i] <- "ns"
          }
        }
      }
    }else{if (length(valid_p_indices) > 0) {
      results$p_adjusted <- NA
      results$p_adjusted[valid_p_indices] <- p.adjust(results$p_value[valid_p_indices], method = "fdr")
    } else {
      results$p_adjusted <- NA
    }
      
      for (i in 1:n_bins) {
        if (is.na(results$p_adjusted[i])) {
          results$significance[i] <- "NA"
        } else {
          if (results$p_adjusted[i] < 0.001) {
            results$significance[i] <- "***"
          } else if (results$p_adjusted[i] < 0.01) {
            results$significance[i] <- "**"
          } else if (results$p_adjusted[i] < 0.05) {
            results$significance[i] <- "*"
          } else {
            results$significance[i] <- "ns"
          }
        }
      }
    }
    
    results$direction <- ifelse(results$significance == "ns", "ns", 
                                ifelse(results$significance == "NA", "NA",
                                       ifelse(results$mean_1 < results$mean_0, "-", "+")))
    
    results$direction[is.na(results$direction)] <- "NA"
    results$mean_0[is.na(results$mean_0)] <- 0
    results$mean_1[is.na(results$mean_1)] <- 0
    
    results$direction <- factor(results$direction, levels = c("-", "ns", "+", "NA"))
    
    bin_median_pseudotime <- tapply(analysis_data$pseudotime, analysis_data$pseudotime_bin, median, na.rm = TRUE)
    
    results$median_pseudotime <- bin_median_pseudotime[as.character(results$bin)]
    
    p <- ggplot(results, aes(x = median_pseudotime, y = mean_1 - mean_0)) +
      geom_point(data = subset(results, direction != "NA"),
                 aes(shape = direction, color = direction), size = 6) +
      geom_text(data = subset(results, direction == "NA"),
                aes(label = "NA"), color = "black", size = 6, fontface = "bold") +
      scale_shape_manual(values = c("-" = 45, "ns" = 1, "+" = 43, "NA" = 4),
                         labels = c("lower (-)", "No significant", "higher (+)", "NA"),
                         drop = FALSE) +
      scale_color_manual(values = c("-" = "blue", "ns" = "black", "+" = "red", "NA" = "black"),
                         labels = c("lower (-)", "No significant", "higher (+)", "NA"),
                         drop = FALSE) +
      ylim(-max(abs(results$mean_1 - results$mean_0), na.rm = TRUE) - 0.2, 
           max(abs(results$mean_1 - results$mean_0), na.rm = TRUE) + 0.2) +
      labs(x = "Pseudotime",
           y = "Expression Difference") +
      theme_light() +
      guides(
        shape = guide_legend(override.aes = list(
          shape = c(45, 1, 43, 4),
          color = c("blue", "black", "red", "black"),
          size = rep(4, 4)
        ), title = "Significance"),
        color = "none"  
      )+ theme(
        legend.title = element_text(size = 0),     
        legend.text = element_text(size = 10),      
        axis.title.x = element_text(size = 14),      
        axis.title.y = element_text(size = 14),     
        axis.text.x = element_text(size = 12),       
        axis.text.y = element_text(size = 12),
        legend.position = "bottom"
      )
    gene_plot2(p)
    return(p)
  })
  
  raw_combined_samples <- reactive({
    unique(c(
      if(is.null(input$selected_samples_left)) character(0) else input$selected_samples_left,
      if(is.null(input$selected_samples_right)) character(0) else input$selected_samples_right
    ))
  })
  
  combined_samples <- debounce(raw_combined_samples, 500)
  
  output$gene_plot3 <- renderPlot({
    # Require essential inputs
    req(gene_analysis_selected(), xde_result())
    
    # Get the debounced combined samples
    selected_samples <- combined_samples()
    
    # Handle empty selection case
    if (length(selected_samples) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "Please select at least one sample", 
                        size = 6) + 
               theme_void())
    }
    
    # Prepare expression data
    gene_name <- gene_analysis_selected()
    long_expr <- data.frame(
      Cell = colnames(xde_result()$expr),
      expr = as.numeric(xde_result()$expr[gene_name, colnames(xde_result()$expr)]),
      pseudotime = xde_result()$pseudotime[colnames(xde_result()$expr)]
    )
    
    # Merge with cell annotations
    cellanno_sub <- xde_result()$cellanno[
      xde_result()$cellanno$Cell %in% colnames(xde_result()$expr), ]
    plot_df <- merge(cellanno_sub, long_expr, by = "Cell")
    
    # Filter for selected samples
    plot_df <- plot_df[plot_df$Sample %in% selected_samples, ]
    
    # Prepare color scheme
    samples_all <- unique(xde_result()$cellanno$Sample)
    colors_all <- color_selected(length(samples_all))
    names(colors_all) <- samples_all
    
    # Create base plot
    p <- ggplot(plot_df, aes(x = pseudotime, y = expr)) +
      theme_light() +
      labs(x = "Pseudotime", y = "Expression Level") + 
      theme(
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "right"
      )
    
    # Add smoothed lines with consistent legend
    if (length(selected_samples) == 1) {
      p <- p + 
        geom_smooth(
          aes(color = "Selected Sample", fill = "Selected Sample"),  
          se = input$show_error_bar, 
          method = "gam",
          formula = y ~ s(x, bs = "cs"),
          show.legend = TRUE,
          alpha = 0.2,
          level = ifelse(input$show_error_bar, input$confidence_level / 100, 0.95)#AD
        ) +
        scale_color_manual(
          name = "Samples",
          values = setNames(colors_all[selected_samples], "Selected Sample"),
          labels = setNames(selected_samples, "Selected Sample")
        ) +
        scale_fill_manual(  
          name = "Samples",
          values = setNames(alpha(colors_all[selected_samples], 0.2), "Selected Sample"),
          labels = setNames(selected_samples, "Selected Sample")
        )
    } else {
      fill_colors <- alpha(colors_all[selected_samples], 0.2)
      names(fill_colors) <- selected_samples
      
      p <- p + 
        geom_smooth(
          aes(color = Sample, fill = Sample),  
          se = input$show_error_bar,
          method = "gam",
          formula = y ~ s(x, bs = "cs"),
          show.legend = TRUE,
          alpha = 0.2,
          level = ifelse(input$show_error_bar, input$confidence_level / 100, 0.95)#AD
        ) +
        scale_color_manual(values = colors_all[selected_samples]) +
        scale_fill_manual(values = fill_colors) 
    }
    
    # Conditionally add points
    if (show_points()) {
      if (length(selected_samples) == 1) {
        p <- p + geom_point(
          aes(color = "Selected Sample"),
          alpha = 0.3, 
          size = 1,
          show.legend = FALSE
        )
      } else {
        p <- p + geom_point(
          aes(color = Sample),
          alpha = 0.3, 
          size = 1,
          show.legend = FALSE
        )
      }
    }
    gene_plot3(p)
    return(p+theme(legend.position = "none"))
  })
  
  output$gene_plot4 <- renderPlot({
    req(gene_analysis_selected(), sce(), xde_result())
    
    gene_name <- gene_analysis_selected()
    reduced_dim <- reducedDim(sce(), "UMAP")
    expr_values <- as.numeric(xde_result()$expr[gene_name, colnames(sce())])
    
    plot_df <- data.frame(
      Dim1 = reduced_dim[,1],
      Dim2 = reduced_dim[,2],
      Expression = expr_values,
      Sample = factor(sce()$sample) 
    )
    
    all_samples <- levels(plot_df$Sample)
    
    if (!included_zero()) {
      plot_df <- plot_df[plot_df$Expression > 0, ]
      plot_df$Sample <- factor(plot_df$Sample, levels = all_samples)
    }
    
    p <- ggplot(plot_df, aes(x = Dim1, y = Dim2)) +
      {
        if (nrow(plot_df) > 0) {
          geom_point(aes(color = Expression), size = 0.5, alpha = 0.8)
        } else {
          geom_blank()
        }
      } +
      scale_color_viridis_c() +
      labs(x = colnames(reduced_dim)[1], 
           y = colnames(reduced_dim)[2]) +
      theme_light() +
      facet_wrap(~ Sample, nrow = 1, drop = FALSE) +  
      {if (included_zero() && length(setdiff(all_samples, unique(plot_df$Sample)))) {
        geom_text(
          data = data.frame(Sample = setdiff(all_samples, unique(plot_df$Sample))),
          x = mean(range(reduced_dim[,1])),
          y = mean(range(reduced_dim[,2])),
          label = "No expression",
          color = "black",
          inherit.aes = FALSE,
          label.size = 0.5
        )
      }
      } +
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        theme(legend.position = "right", legend.width = unit(2, "inches"))
      )
    
    gene_plot4(p)
    return(p)
  })
  
  observeEvent(input$save_plot, {
    req(input$result_folder)  
    if (is.null(result_visulization_folder())) {
      time_str <- format(Sys.time(), "%Y%m%d%H%M")
      vis_dir <- file.path(input$result_folder, paste0("visualization_", time_str))
      if (!dir.exists(vis_dir)) dir.create(vis_dir, recursive = TRUE)
      result_visulization_folder(vis_dir)
    }
    vis_dir <- result_visulization_folder()
    
    tryCatch({
      if (is.na(gene_analysis_selected())) {
        reduction_dir <- file.path(vis_dir, "reduction")
        if (!dir.exists(reduction_dir)) dir.create(reduction_dir)
        
        pdf(file.path(reduction_dir, "Cluster.pdf"), width = 8, height = 6)
        print(gene_plot1_NA())
        dev.off()
        
        pdf(file.path(reduction_dir, "Pseudotime.pdf"), width = 8, height = 6)
        print(gene_plot2_NA())
        dev.off()
        
        sample_number <- length(unique(xde_result()$cellanno$Sample))
        pdf(file.path(reduction_dir, "Combined_Plots.pdf"), width = 2 * (13.5 / sample_number) + 4, height = 6)
        plot1_with_theme <- gene_plot1_NA() + theme(legend.position = "right", legend.width = unit(2, "inches"))
        plot2_with_theme <- gene_plot2_NA() + theme(legend.position = "right", legend.width = unit(2, "inches"))
        
        grid.arrange(plot1_with_theme, plot2_with_theme, nrow = 1)
        dev.off()
        
        showNotification("Cluster and Pseudotime plots saved successfully!", type = "message")
      } else {
        gene_dir <- file.path(vis_dir, gene_analysis_selected())
        if (!dir.exists(gene_dir)) dir.create(gene_dir)
        combined_dir <- file.path(vis_dir, "COMBINED")
        if (!dir.exists(combined_dir)) dir.create(combined_dir)
        
        pdf(file.path(gene_dir, "PopulationPlot.pdf"), width = 8, height = 6)
        print(gene_plot1())
        dev.off()
        
        pdf(file.path(gene_dir, "Bin.pdf"), width = 8, height = 6)
        print(gene_plot2())
        dev.off()
        
        pdf(file.path(gene_dir, "SamplePlot.pdf"), width = 8, height = 6)
        print(gene_plot3())
        dev.off()
        
        pdf(file.path(gene_dir, "ReductionSplitPlot.pdf"), width = 16, height = 6)
        print(gene_plot4())
        dev.off()
        
        write.xlsx(gene_plot2()$data, file.path(gene_dir, "Bin_result.xlsx"))
        
        pdf(file.path(combined_dir, paste0(gene_analysis_selected(), ".pdf")), width = 20, height = 12)
        
        grid.newpage()

        pushViewport(viewport(x = 0.05, y = 0.98, width = 0.15, height = 0.95, just = c("left","top")))
        
        gene_stats <- xde_result()$statistics[gene_analysis_selected(), , drop = FALSE]
        info_text <- paste0(
          "Gene:\n\t\t", gene_analysis_selected(), "\n\n",
          paste(
            sapply(colnames(gene_stats), function(col) {
              paste0(col, ":\n\t\t", formatC(gene_stats[, col], format = "e", digits = 2))
            }),
            collapse = "\n\n"
          )
        )
        
        grid.text(
          info_text,
          x = 0,                 
          y = 1,                 
          gp = gpar(col = "darkred", fontface = "bold", cex = 1.2),
          just = c("left", "top") 
        )
        
        popViewport() 

        
        # p1
        pushViewport(viewport(x = 0.18, y = 0.98, width = 0.24, height = 0.45, just = c("left","top")))
        print(gene_plot1(), newpage = FALSE)
        popViewport()
        
        # p2
        pushViewport(viewport(x = 0.42, y = 0.98, width = 0.24, height = 0.45, just = c("left","top")))
        print(gene_plot2(), newpage = FALSE)
        popViewport()
        
        # p3
        pushViewport(viewport(x = 0.66, y = 0.98, width = 0.24, height = 0.405, just = c("left","top")))
        p3 <- gene_plot3()
        p3_no_legend <- p3 + theme(legend.position = "none")
        print(p3_no_legend, newpage = FALSE)
        popViewport()
        
        pushViewport(viewport(x = 0.9, y = 0.98, width = 0.1, height = 0.405, just = c("left","top")))
        legend <- get_legend(p3)
        grid.draw(legend)
        popViewport()
        
        # p4
        pushViewport(viewport(x = 0.18, y = 0.5, width = 0.8, height = 0.48, just = c("left","top")))
        print(gene_plot4(), newpage = FALSE)
        popViewport()
        
        dev.off()
        
        
        
        showNotification(paste("Gene", gene_analysis_selected(), "plots saved successfully!"), 
                         type = "message")
      }
    }, error = function(e) {
      showNotification(paste("Save failed:", e$message), type = "error")
    })
  })
  
  #MULTI GENE UI
  observeEvent(input$load_default_multi, {
    sig_genes <- rownames(stat()[stat()$fdr.overall < 0.05, ])
    updateTextAreaInput(
      session, 
      "multi_gene_list",
      value = if(length(sig_genes) == 0){"NO SIGINIFICANT GENE"}else{paste(sig_genes, collapse = "\n")}
    )
  })
  
  observeEvent(input$clear_genes_multi, {
    updateTextAreaInput(session, "multi_gene_list", value = "")
  })
  
  observeEvent(input$save_multi_gene, {
    req(input$multi_gene_list)
    if (is.null(result_visulization_folder())) {
      time_str <- format(Sys.time(), "%Y%m%d%H%M")
      vis_dir <- file.path(input$result_folder, paste0("visualization_", time_str))
      if (!dir.exists(vis_dir)) dir.create(vis_dir, recursive = TRUE)
      result_visulization_folder(vis_dir)
    }
    
    combined_dir <- file.path(result_visulization_folder(), "COMBINED")
    if (!dir.exists(combined_dir)) dir.create(combined_dir, recursive = TRUE)
    
    gene_list <- unlist(strsplit(input$multi_gene_list, "\n"))
    gene_list <- trimws(gene_list)  
    gene_list <- gene_list[gene_list != ""]  
    
    if (length(gene_list) == 0) {
      showNotification("No genes to analyze!", type = "warning")
      return()
    }
    
    progress_log <- paste0(Sys.time(), ": Starting multi-gene analysis...\n")
    
    showModal(modalDialog(
      title = "Processing Genes",
      tagList(
        div(style = "height: 400px; overflow-y: auto; background: #f8f9fa; padding: 10px; border-radius: 5px;",
            pre(id = "progress-log", style = "white-space: pre-wrap; word-break: break-all; font-family: monospace;",
                progress_log)
        )
      ),
      footer = NULL,
      size = "l",
      easyClose = FALSE
    ))
    
    show_points_value <- input$show_points_multi
    included_zero_value <- input$included_zero_multi
    tryCatch({
      result_summary <- data.frame(
        Gene = character(),
        Number_of_Minus_Bins = integer(),
        Number_of_NonSig_Bins = integer(),
        Number_of_Plus_Bins = integer(),
        Number_of_NA_Bins = integer(),
        stringsAsFactors = FALSE
      )
      for (i in seq_along(gene_list)) {
        
        gene <- gene_list[i]
        
        progress_log <- paste0(progress_log, Sys.time(), ": Processing gene '", gene, "' (", i, "/", length(gene_list), ")...\n")
        
        shinyjs::html(id = "progress-log", html = progress_log, add = FALSE)
        
        shinyjs::runjs("var elem = document.getElementById('progress-log'); elem.scrollTop = elem.scrollHeight;")
        
        Sys.sleep(0.05)
        shinyjs::runjs("1+1")  
        
        tryCatch({
          if (!gene %in% rownames(xde_result()$statistics)) {
            progress_log <- paste0(progress_log, Sys.time(), ": Gene '", gene, "' not found in results!\n")
            next
          }
          
          # Create PDF file for this gene
          pdf_file <- file.path(combined_dir, paste0(gene, ".pdf"))
          pdf(pdf_file, width = 20, height = 12)
          
          # Ensure PDF device closes even on error
          on.exit({
            if (dev.cur() > 1) dev.off()
          }, add = TRUE)
          
          # Gene info panel
          pushViewport(viewport(x = 0.05, y = 0.98, width = 0.15, height = 0.95, just = c("left","top")))
          
          gene_stats <- xde_result()$statistics[gene, , drop = FALSE]
          info_text <- paste0(
            "Gene:\n", gene, "\n\n",
            paste(
              sapply(colnames(gene_stats), function(col) {
                paste0(col, ":\n\t\t", formatC(gene_stats[, col], format = "e", digits = 2))
              }),
              collapse = "\n\n"
            )
          )
          
          grid.text(
            info_text,
            x = 0,                 
            y = 1,                 
            gp = gpar(col = "darkred", fontface = "bold", cex = 1.2),
            just = c("left", "top") 
          )
          
          popViewport() 
          
          # Population plot
          pushViewport(viewport(x = 0.18, y = 0.98, width = 0.24, height = 0.45, just = c("left","top")))
          p1 <- plotGenePopulation(xde_result(), gene, type = "variable") + theme(
            legend.title = element_text(size = 14),     
            legend.text = element_text(size = 10),      
            axis.title.x = element_text(size = 14),      
            axis.title.y = element_text(size = 14),     
            axis.text.x = element_text(size = 12),       
            axis.text.y = element_text(size = 12),
            legend.position = "bottom"
          )
          print(p1, newpage = FALSE)
          popViewport()
          
          # Bin Plot
          pushViewport(viewport(x = 0.42, y = 0.98, width = 0.24, height = 0.45, just = c("left","top")))
          
          n_bins <- input$multi_bin_number
          gene_expr <- xde_result()$expr[gene, ]
          pseudotime_bins <- cut(xde_result()$pseudotime, breaks = n_bins, labels = FALSE)
          
          analysis_data <- data.frame(
            expr = as.numeric(gene_expr),
            cell_id = names(gene_expr),
            pseudotime_bin = pseudotime_bins,
            pseudotime = xde_result()$pseudotime  
          )
          
          analysis_data <- merge(analysis_data, xde_result()$cellanno, by.x = "cell_id", by.y = "Cell")
          analysis_data <- merge(analysis_data, xde_result()$design, by.x = "Sample", by.y = "row.names")
          results <- data.frame(
            bin = 1:n_bins,
            mean_0 = numeric(n_bins),
            mean_1 = numeric(n_bins),
            p_value = numeric(n_bins),
            significance = character(n_bins),
            stringsAsFactors = FALSE
          )
          
          
          for (i in 1:n_bins) {
            bin_data <- analysis_data[analysis_data$pseudotime_bin == i, ]
            
            intercept_index <- which(colnames(bin_data) == "intercept")
            treatment_col <- colnames(bin_data)[intercept_index + 1]
            
            confounder_cols <- colnames(bin_data)[(intercept_index + 2):ncol(bin_data)]
            
            if (length(confounder_cols) > 0) {
              model_formula <- as.formula(paste("expr ~", treatment_col, "+", 
                                                paste(confounder_cols, collapse = " + "),
                                                "+ (1 | Sample)"))
            } else {
              model_formula <- as.formula(paste("expr ~", treatment_col, "+ (1 | Sample)"))
            }
            
            if (nrow(bin_data) > n_bins) {  
              tryCatch({
                suppressMessages(suppressWarnings({
                  model <- lmer(model_formula, data = bin_data)
                }))
                
                model_summary <- summary(model)
                treatment_coef <- which(rownames(model_summary$coefficients) == treatment_col)
                
                if (length(treatment_coef) > 0) {
                  results$p_value[i] <- model_summary$coefficients[treatment_coef, "Pr(>|t|)"]
                } else {
                  results$p_value[i] <- NA
                }
                
              }, error = function(e) {
                results$p_value[i] <- NA
              })
            } else {
              results$p_value[i] <- NA
            }
            
            expr_0 <- bin_data$expr[bin_data[[treatment_col]] == 0]
            expr_1 <- bin_data$expr[bin_data[[treatment_col]] == 1]
            results$mean_0[i] <- mean(expr_0, na.rm = TRUE)
            results$mean_1[i] <- mean(expr_1, na.rm = TRUE)
          }
          valid_p_indices <- which(!is.na(results$p_value))
          
          if(input$multi_using_p_adj == "pvalue"){
            for (i in 1:n_bins) {
              if (is.na(results$p_value[i])) {
                results$significance[i] <- "NA"
              } else {
                if (results$p_value[i] < 0.001) {
                  results$significance[i] <- "***"
                } else if (results$p_value[i] < 0.01) {
                  results$significance[i] <- "**"
                } else if (results$p_value[i] < 0.05) {
                  results$significance[i] <- "*"
                } else {
                  results$significance[i] <- "ns"
                }
              }
            }
          }else{if (length(valid_p_indices) > 0) {
            results$p_adjusted <- NA
            results$p_adjusted[valid_p_indices] <- p.adjust(results$p_value[valid_p_indices], method = "fdr")
          } else {
            results$p_adjusted <- NA
          }
            
            for (i in 1:n_bins) {
              if (is.na(results$p_adjusted[i])) {
                results$significance[i] <- "NA"
              } else {
                if (results$p_adjusted[i] < 0.001) {
                  results$significance[i] <- "***"
                } else if (results$p_adjusted[i] < 0.01) {
                  results$significance[i] <- "**"
                } else if (results$p_adjusted[i] < 0.05) {
                  results$significance[i] <- "*"
                } else {
                  results$significance[i] <- "ns"
                }
              }
            }
          }
          
          results$direction <- ifelse(results$significance == "ns", "ns", 
                                      ifelse(results$significance == "NA", "NA",
                                             ifelse(results$mean_1 < results$mean_0, "-", "+")))
          
          results$direction[is.na(results$direction)] <- "NA"
          results$mean_0[is.na(results$mean_0)] <- 0
          results$mean_1[is.na(results$mean_1)] <- 0
          
          results$direction <- factor(results$direction, levels = c("-", "ns", "+", "NA"))
          
          bin_median_pseudotime <- tapply(analysis_data$pseudotime, analysis_data$pseudotime_bin, median, na.rm = TRUE)
          
          results$median_pseudotime <- bin_median_pseudotime[as.character(results$bin)]
          direction_counts <- table(factor(results$direction, levels = c("-", "ns", "+", "NA")))
          result_summary <- rbind(result_summary, data.frame(
            Gene = gene,
            Number_of_Minus_Bins = as.integer(direction_counts["-"]),
            Number_of_NonSig_Bins = as.integer(direction_counts["ns"]),
            Number_of_Plus_Bins = as.integer(direction_counts["+"]),
            Number_of_NA_Bins = as.integer(direction_counts["NA"]),
            stringsAsFactors = FALSE
          ))
          p2 <- ggplot(results, aes(x = median_pseudotime, y = mean_1 - mean_0)) +
            geom_point(data = subset(results, direction != "NA"),
                       aes(shape = direction, color = direction), size = 6) +
            geom_text(data = subset(results, direction == "NA"),
                      aes(label = "NA"), color = "black", size = 6, fontface = "bold") +
            scale_shape_manual(values = c("-" = 45, "ns" = 1, "+" = 43, "NA" = 4),
                               labels = c("lower (-)", "No significant", "higher (+)", "NA"),
                               drop = FALSE) +
            scale_color_manual(values = c("-" = "blue", "ns" = "black", "+" = "red", "NA" = "black"),
                               labels = c("lower (-)", "No significant", "higher (+)", "NA"),
                               drop = FALSE) +
            ylim(-max(abs(results$mean_1 - results$mean_0), na.rm = TRUE) - 0.2, 
                 max(abs(results$mean_1 - results$mean_0), na.rm = TRUE) + 0.2) +
            labs(x = "Pseudotime",
                 y = "Expression Difference") +
            theme_light() +
            guides(
              shape = guide_legend(override.aes = list(
                shape = c(45, 1, 43, 4),
                color = c("blue", "black", "red", "black"),
                size = rep(4, 4)
              ), title = "Significance"),
              color = "none"  
            )+ theme(
              legend.title = element_text(size = 0),     
              legend.text = element_text(size = 10),      
              axis.title.x = element_text(size = 14),      
              axis.title.y = element_text(size = 14),     
              axis.text.x = element_text(size = 12),       
              axis.text.y = element_text(size = 12),
              legend.position = "bottom"
            )
          print(p2, newpage = FALSE)  
          popViewport()
          
          # Pseudotime plot
          pushViewport(viewport(x = 0.66, y = 0.98, width = 0.24, height = 0.405, just = c("left","top")))
          long_expr <- data.frame(
            Cell = colnames(xde_result()$expr),
            expr = as.numeric(xde_result()$expr[gene, colnames(xde_result()$expr)]),
            pseudotime = xde_result()$pseudotime[colnames(xde_result()$expr)]
          )
          
          cellanno_sub <- xde_result()$cellanno[xde_result()$cellanno$Cell %in% colnames(xde_result()$expr), ]
          plot_df <- merge(cellanno_sub, long_expr, by = "Cell")
          
          samples_all <- unique(xde_result()$cellanno$Sample)
          colors_all <- color_selected(length(samples_all))
          names(colors_all) <- samples_all
          
          p3 <- ggplot(plot_df, aes(x = pseudotime, y = expr, color = Sample)) +
            scale_color_manual(values = colors_all) + 
            theme_light() +
            labs(x = "Pseudotime", y = "Expression Level") + 
            theme(
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              plot.title = element_text(size = 16, hjust = 0.5),
              legend.position = "right"
            )
          p3 <- p3 + 
            geom_smooth(
              aes(color = Sample, fill = Sample),  
              se = input$show_error_bar_multi,
              method = "gam",
              formula = y ~ s(x, bs = "cs"),
              show.legend = TRUE,
              alpha = 0.2,
              level = ifelse(input$show_error_bar_multi, input$confidence_level_multi / 100, 0.95) #AD
            ) 
          if (show_points_value) {
            p3 <- p3 + geom_point(alpha = 0.3, size = 1)
          }
          print(p3 + theme(legend.position = "none"), newpage = FALSE)  
          popViewport()
          
          #legend
          pushViewport(viewport(x = 0.9, y = 0.98, width = 0.1, height = 0.405, just = c("left","top")))
          legend <- get_legend(p3)
          grid.draw(legend)
          popViewport()
          
          
          # Reduction plot
          pushViewport(viewport(x = 0.18, y = 0.5, width = 0.8, height = 0.48, just = c("left","top")))
          reduced_dim <- reducedDim(sce(), "UMAP")
          expr_values <- as.numeric(xde_result()$expr[gene, colnames(sce())])
          
          plot_df <- data.frame(
            Dim1 = reduced_dim[,1],
            Dim2 = reduced_dim[,2],
            Expression = expr_values,
            Sample = factor(sce()$sample) 
          )
          
          all_samples <- levels(plot_df$Sample)
          
          if (!included_zero_value) {
            plot_df <- plot_df[plot_df$Expression > 0, ]
            plot_df$Sample <- factor(plot_df$Sample, levels = all_samples)
          }
          
          p4 <- ggplot(plot_df, aes(x = Dim1, y = Dim2)) +
            {
              if (nrow(plot_df) > 0) {
                geom_point(aes(color = Expression), size = 0.5, alpha = 0.8)
              } else {
                geom_blank()
              }
            } +
            scale_color_viridis_c() +
            labs(x = colnames(reduced_dim)[1], 
                 y = colnames(reduced_dim)[2]) +
            theme_light() +
            facet_wrap(~ Sample, nrow = 1, drop = FALSE) +  
            {if (included_zero_value && length(setdiff(all_samples, unique(plot_df$Sample)))) {
              geom_text(
                data = data.frame(Sample = setdiff(all_samples, unique(plot_df$Sample))),
                x = mean(range(reduced_dim[,1])),
                y = mean(range(reduced_dim[,2])),
                label = "No expression",
                color = "black",
                inherit.aes = FALSE,
                label.size = 0.5
              )
            }} +
            theme(
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              strip.text.x = element_text(size = 10)
            )
          
          print(p4, newpage = FALSE) 
          popViewport()
          
          dev.off()
          
          progress_log <- paste0(progress_log, Sys.time(), ": Gene '", gene, "' finished.\n")
          
        }, error = function(e) {
          if (dev.cur() > 1) dev.off()
          progress_log <<- paste0(progress_log, Sys.time(), ": Error processing gene '", gene, "' - ", e$message, "\n")
        })
        
        shinyjs::html(id = "progress-log", html = progress_log, add = FALSE)
        shinyjs::runjs("var elem = document.getElementById('progress-log'); elem.scrollTop = elem.scrollHeight;")
        
        Sys.sleep(0.05)
      }
      
      stats_data <- xde_result()$statistics
      additional_stats <- as.data.frame(stats_data[, -1, drop = FALSE])  
      gene_indices <- match(result_summary$Gene, rownames(additional_stats))
      
      stat_colnames <- colnames(additional_stats)
      for (col in stat_colnames) {
        result_summary[[col]] <- NA 
      }
    
      valid_matches <- !is.na(gene_indices)
      if (any(valid_matches)) {
        for (col in stat_colnames) {
          result_summary[valid_matches, col] <- additional_stats[gene_indices[valid_matches], col]
        }
      }
      
      description_df <- data.frame(
        Colname = c(
          "Number_of_Minus_Bins",
          "Number_of_NonSig_Bins", 
          "Number_of_Plus_Bins",
          "Number_of_NA_Bins",
          "pval.overall",
          "z.overall",
          "log.pval.overall",
          "fdr.meanDiff", 
          "pval.meanDiff",
          "z.meanDiff",
          "log.pval.meanDiff",
          "fdr.trendDiff",
          "pval.trendDiff", 
          "z.trendDiff",
          "log.pval.trendDiff"
        ),
        Description = c(
          "Number of bins where test.variable=1 mean is significantly less than test.variable=0 mean (p < 0.05)",
          "Number of bins with no significant difference between conditions (p ≥ 0.05)",
          "Number of bins where test.variable=1 mean is significantly greater than test.variable=0 mean (p < 0.05)", 
          "Number of bins where statistical test could not be performed (missing data or model convergence issues)",
          "Raw p-value from combined test considering both meanDiff and trendDiff effects",
          "Z-score from combined test quantifying standardized deviation from null distribution",
          "Log-transformed p-value from combined test using kernel density estimation",
          "False discovery rate adjusted p-values for mean difference test (Benjamini-Hochberg correction) [Recommend]",
          "Raw p-value indicating statistical significance of absolute expression differences between conditions (vertical curve shifts)",
          "Z-score for mean difference test quantifying standardized deviation from null distribution",
          "Log-transformed p-value for mean difference test using kernel density estimation",
          "False discovery rate adjusted p-values for trend difference test (Benjamini-Hochberg correction) [Recommend]",
          "Raw p-value indicating statistical significance of pattern/shape differences between conditions (slope changes)",
          "Z-score for trend difference test quantifying standardized deviation from null distribution", 
          "Log-transformed p-value for trend difference test using kernel density estimation"
        ),
        stringsAsFactors = FALSE
      )
      
      wb <- createWorkbook()
      
      addWorksheet(wb, "Data")
      writeData(wb, sheet = "Data", result_summary, rowNames = FALSE)
      
      addWorksheet(wb, "Description_of_Column")
      writeData(wb, sheet = "Description_of_Column", description_df, rowNames = FALSE)
      
      setColWidths(wb, sheet = "Description_of_Column", cols = 1, widths = 20)  
      setColWidths(wb, sheet = "Description_of_Column", cols = 2, widths = 80)  
      saveWorkbook(wb, file.path(combined_dir, "gene_direction_summary.xlsx"), overwrite = TRUE)
      progress_log <- paste0(progress_log, Sys.time(), ": DONE! Processed ", length(gene_list), " genes.\n The window will automatically close after 3 seconds.")

      
      shinyjs::html(id = "progress-log", html = progress_log, add = FALSE)
      shinyjs::runjs("var elem = document.getElementById('progress-log'); elem.scrollTop = elem.scrollHeight;")
      
      delay(3000, {
        removeModal()
        
      })
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error in processing:", e$message), type = "error")
    })
  })
  
  #DOWNLOAD UI
  output$download_visualization <- downloadHandler(
    filename = function() {
      paste0("visualization_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      if (is.null(result_visulization_folder())) {
        time_str <- format(Sys.time(), "%Y%m%d%H%M")
        vis_dir <- file.path(input$result_folder, paste0("visualization_", time_str))
        if (!dir.exists(vis_dir)) dir.create(vis_dir, recursive = TRUE)
        result_visulization_folder(vis_dir)
      }
      
      vis_dir <- result_visulization_folder()
      req(dir.exists(vis_dir))
      showModal(modalDialog(
        title = "Processing",
        div(style = "text-align: center;",
            h4("Preparing visualization files for download..."),
            tags$div(style = "margin-top: 20px;")  
        ),
        footer = NULL,
        easyClose = FALSE,
        size = "s"
      ))
      
      temp_zip <- tempfile(fileext = ".zip")
      tryCatch({
        zip::zipr(
          zipfile = temp_zip,
          files = vis_dir,
          include_directories = FALSE
        )
        
        file.copy(temp_zip, file)
        
        removeModal()
        
        showModal(modalDialog(
          title = "Download Complete",
          div(style = "font-size: 16px;",
              p("Download successful!"),
              p("Please check your default download folder."),
              p("You can also find the files in:"),
              tags$pre(style = "background: #f5f5f5; padding: 10px; border-radius: 5px;",
                       vis_dir)
          ),
          footer = modalButton("OK"),
          easyClose = TRUE
        ))
      }, error = function(e) {
        removeModal()
        showModal(modalDialog(
          title = "Error",
          div(style = "color: red;",
              p("Failed to create download package:"),
              tags$pre(e$message)
          ),
          footer = modalButton("Close"),
          easyClose = TRUE
        ))
        return(NULL)
      })
    },
    contentType = "application/zip"
  )
}