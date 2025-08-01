dataVisualizationUI <- function() {
  fluidPage(
    shinyjs::useShinyjs(),
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
  show_points <- reactiveVal(FALSE)
  remove_zero <- reactiveVal(FALSE)
  current_dims <- reactiveVal(c(1, 2))
  gene_list <- reactiveVal()
  heatmap_plot <- reactiveVal(NULL)
  heatmap_data <- reactiveVal(NULL)
  
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
      navlistPanel(
        widths = c(2, 10),
        tabPanel("HEATMAP", uiOutput("heatmap_ui")),
        tabPanel("GENE ANALYSIS", uiOutput("gene_ui"))
      )
    }
  })
  
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
  
  output$heatmap_ui <- renderUI({
    sig_genes <- rownames(stat()[stat()$fdr.overall < 0.05, ])
    fluidRow(
      column(
        width = 3,
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
          fluidRow(
            column(6, actionButton("load_default", "fdr.overall < 0.05", 
                                   class = "btn-primary", width = "100%")),
            column(6, actionButton("clear_genes", "Clear All", 
                                   class = "btn-danger", width = "100%"))
          ),
          
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
          
          selectInput(
            "scale_method", 
            "Scale:",
            choices = c("Together", "Group Split (TREND)"),
            selected = "Together"
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
          
          actionButton("save_heatmap_plot", "Save Plot", 
                       class = "btn-primary", width = "100%", disabled = TRUE),
          
          hr(),
          
          numericInput(
            "pseudotime_samples",
            "Number of Pseudotime samples:",
            value = 10,
            min = 3,
            max = 100,
            step = 1
          ),
         
          actionButton("save_heatmap_data", "Save Data", 
                       class = "btn-info", width = "100%", disabled = TRUE)
        )
      ),
      
      column(
        width = 9,
        h4("Heatmap", class = "text-primary"),
        hr(),
        plotOutput("heatmap_display")  
      )
    )
  })
  
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
  
  observeEvent(input$submit_genes, {
    req(input$heatmap_gene_list)
    
    genes <- strsplit(input$heatmap_gene_list, "\n")[[1]] %>% 
      trimws() %>% 
      .[. != ""]  
    
    if (length(genes) > 0) {
      gene_list(genes)
      showNotification(paste("Submitted", length(genes), "genes"), type = "message")
    } else {
      showNotification("Gene list is empty!", type = "warning")
    }
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
        scale_method = input$scale_method,
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
    default_val <- ifelse(input$cluster_method == "kmeans", 6, 0)
    updateNumericInput(session, "cluster_number", value = default_val)
  })
  
  output$gene_ui <- renderUI({
    fluidRow(
      column(
        width = 3,
        style = "border-right: 1px solid #eee; padding-right: 15px;",
        if(!is.na(gene_analysis_selected())) {
          actionButton("original_btn", "Reset", 
                       class = "btn-default", 
                       style = "width: 100%; margin-bottom: 15px;")
        },
        actionButton("select_gene", "Select Gene", 
                     class = "btn-success", 
                     style = "width: 100%; margin-bottom: 15px;"),
        actionButton("save_plot", "Save Plot to result folder", 
                     class = "btn-primary", 
                     style = "width: 100%; margin-bottom: 15px;"),
        uiOutput("selected_gene_info")
      ),
      column(
        width = 9,
        uiOutput("gene_plot_ui", height = "800px")
      )
    )
  })
  
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
      current <- current_dims()  
      
      tagList(
        fluidRow(
          column(6, selectInput("dim_x", "X Axis:", 
                                choices = dim_names, 
                                selected = dim_names[current[1]], width = "100%")),  
          column(6, selectInput("dim_y", "Y Axis:", 
                                choices = dim_names, 
                                selected = dim_names[current[2]], width = "100%"))
        ),
        fluidRow(
          column(6, plotOutput("gene_plot1_NA", height = "400px")),
          column(6, plotOutput("gene_plot2_NA", height = "400px"))
        )
      )
    } else {
      tagList(
        fluidRow(
          column(6, 
                 checkboxInput("show_points", "Show individual data points", 
                               value = FALSE, width = "100%")),
          column(6, 
                 checkboxInput("remove_zero", "Remove zero-expression cells", 
                               value = FALSE, width = "100%"))
        ),
        fluidRow(
          column(6, plotOutput("gene_plot1", height = "400px")),
          column(6, plotOutput("gene_plot2", height = "400px"))
        ),
        fluidRow(
          column(12, plotOutput("gene_plot3", height = "400px"))
        )
      )
    }
  })
  
  observeEvent(c(input$dim_x, input$dim_y), {
    req(input$dim_x, input$dim_y)
    dim_indices <- match(c(input$dim_x, input$dim_y), colnames(reducedDim(sce(), "Reduction")))
    current_dims(dim_indices)
  })
  
  observeEvent(input$show_points, {
    show_points(input$show_points)
  })
  
  observeEvent(input$remove_zero, {
    remove_zero(input$remove_zero)
  })
  
  observeEvent(input$save_plot, {
    req(input$result_folder)  
    vis_dir <- file.path(input$result_folder, "visualization")
    if (!dir.exists(vis_dir)) dir.create(vis_dir)
    
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
        
        showNotification("Cluster and Pseudotime plots saved successfully!", type = "message")
      } else {
        gene_dir <- file.path(vis_dir, gene_analysis_selected())
        if (!dir.exists(gene_dir)) dir.create(gene_dir)
        
        pdf(file.path(gene_dir, "PopulationPlot.pdf"), width = 8, height = 6)
        print(gene_plot1())
        dev.off()
        
        pdf(file.path(gene_dir, "SamplePlot.pdf"), width = 8, height = 6)
        print(gene_plot2())
        dev.off()
        
        pdf(file.path(gene_dir, "ReductionSplitPlot.pdf"), width = 16, height = 6)
        print(gene_plot3())
        dev.off()
        
        showNotification(paste("Gene", gene_analysis_selected(), "plots saved successfully!"), 
                         type = "message")
      }
    }, error = function(e) {
      showNotification(paste("Save failed:", e$message), type = "error")
    })
  })
  
  output$gene_plot1_NA <- renderPlot({
    req(sce(), current_dims())
    dims <- current_dims()
    rd <- reducedDim(sce(), "Reduction")[, dims]  
    
    p <- plotReducedDim(sce(), 
                        dimred = "Reduction", 
                        colour_by = "cluster",
                        ncomponents = dims,  
                        point_alpha = 1) +
      labs(
        x = colnames(rd)[1],  
        y = colnames(rd)[2],  
        title = "Cell Clusters"
      ) +
      theme_minimal() +
      theme(
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 16)
      )
    
    gene_plot1_NA(p)
    return(p)
  })
  
  output$gene_plot2_NA <- renderPlot({
    req(sce(), xde_result(), current_dims())
    dims <- current_dims()
    
    rd <- reducedDim(sce(), "Reduction")[, dims]  
    pt <- xde_result()$pseudotime
    rd <- rd[rownames(rd) %in% names(pt), ]
    
    p <- ggplot(data.frame(rd, pseudotime = pt), 
                aes(x = !!sym(colnames(rd)[1]), 
                    y = !!sym(colnames(rd)[2]), 
                    color = pseudotime)) +
      geom_point(size = 0.5) +
      scale_color_viridis_c() +
      ggtitle("Pseudotime") +
      theme_minimal()+
      theme(legend.title = element_text(size = 14),
            axis.title = element_text(size = 16))
    
    gene_plot2_NA(p)
    return(p)
  })
  
  output$gene_plot1 <- renderPlot({
    req(gene_analysis_selected())
    p <- plotGenePopulation(xde_result(), gene_analysis_selected(), type = "variable") + theme(
      legend.title = element_text(size = 14),     
      legend.text = element_text(size = 12),      
      axis.title.x = element_text(size = 16),      
      axis.title.y = element_text(size = 16),     
      axis.text.x = element_text(size = 14),       
      axis.text.y = element_text(size = 14)        
    )
    gene_plot1(p)
    return(p)
  })
  
  output$gene_plot2 <- renderPlot({
    req(gene_analysis_selected(), xde_result())
    
    gene_name <- gene_analysis_selected()
    long_expr <- data.frame(
      Cell = colnames(xde_result()$expr),
      expr = as.numeric(xde_result()$expr[gene_name, colnames(xde_result()$expr)]),
      pseudotime = xde_result()$pseudotime[colnames(xde_result()$expr)]
    )
    
    cellanno_sub <- xde_result()$cellanno[xde_result()$cellanno$Cell %in% colnames(xde_result()$expr), ]
    plot_df <- merge(cellanno_sub, long_expr, by = "Cell")
    
    p <- ggplot(plot_df, aes(x = pseudotime, y = expr, color = Sample)) +
      geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
      scale_color_brewer(palette = "Set2") + 
      theme_minimal() +
      labs(x = "Pseudotime", y = "Expression Level") + 
      theme(
        legend.title = element_text(size = 14),     
        legend.text = element_text(size = 12),      
        axis.title.x = element_text(size = 16),      
        axis.title.y = element_text(size = 16)      
      )
    
    if (show_points()) {
      p <- p + geom_point(alpha = 0.3, size = 1)
    }
    gene_plot2(p)
    return(p)
  })
  
  output$gene_plot3 <- renderPlot({
  req(gene_analysis_selected(), sce(), xde_result(), current_dims())
  dims <- current_dims()
  
  gene_name <- gene_analysis_selected()
  reduced_dim <- reducedDim(sce(), "Reduction")[, dims]
  expr_values <- as.numeric(xde_result()$expr[gene_name, colnames(sce())])
  
  plot_df <- data.frame(
    Dim1 = reduced_dim[,1],
    Dim2 = reduced_dim[,2],
    Expression = expr_values,
    Sample = factor(sce()$sample) 
  )
  
  all_samples <- levels(plot_df$Sample)
  
  if (remove_zero()) {
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
         y = colnames(reduced_dim)[2],
         title = gene_name) +
    theme_minimal() +
    facet_wrap(~ Sample, nrow = 1, drop = FALSE) +  
    {if (remove_zero() && length(setdiff(all_samples, unique(plot_df$Sample)))) {
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
    }
  
  gene_plot3(p)
  return(p)
})
  
  observeEvent(input$original_btn, {
    gene_analysis_selected(NA)  
  })
}