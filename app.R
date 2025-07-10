library(shiny)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(slingshot)
library(harmony)
library(dplyr)
library(tidyr)
library(rhandsontable)
source("step1.R")
source("step3.R")
source("step4.R")
source("step5.R")

ui <- fluidPage(
  titlePanel("Lamian - Waxman's lab"),
  tags$h4("produced by Bingtian Ye(btye@bu.edu)"),
  tabsetPanel(
    tabPanel("Data preprocessing",
             fluidRow(
               column(8,
                      uiOutput("step_ui"),
                      uiOutput("step_content"),
                      fluidRow(uiOutput("nav_buttons"))
               )
             )
    ),
    tabPanel("Data visualization",
             h3("Coming soon...")
    )
  )
)

server <- function(input, output, session) {
  current_step <- reactiveVal(1)
  total_steps <- 5
  
  seuratObj <- reactiveVal(NULL)
  processedObj <- reactiveVal(NULL)
  dataLoaded <- reactiveVal(FALSE)
  runDoneStep1 <- reactiveVal(FALSE)
  
  data_temp <- reactiveVal(NULL)
  subsetDone <- reactiveVal(FALSE)
  
  step3_result <- reactiveVal(NULL)
  step3_run_done <- reactiveVal(FALSE)
  
  step4_result <- reactiveVal(NULL)
  filter_run_done <- reactiveVal(FALSE)
  
  obj_filtered_temp <- reactiveVal(NULL) # Temporary filtered object for Step4
  
  output$step_ui <- renderUI({
    h4(paste0("Step ", current_step(), " of ", total_steps))
  })
  
  output$step_content <- renderUI({
    step <- current_step()
    if (step == 1) {
      fluidPage(
        fluidRow(
          column(8, textInput("data_path", "Enter path to RDS file", value = "")),
          column(4, br(), actionButton("read_btn_ui", "READ"))
        ),
        checkboxInput("do_harmony", "Perform Harmony correction", TRUE),
        conditionalPanel(
          condition = "input.do_harmony == true",
          uiOutput("harmony_ui")
        )
      )
    } else if (step == 2) {
      req(processedObj())
      obj <- processedObj()
      is_harmony <- isTRUE(input$do_harmony)
      group_by_vars <- if (is_harmony) input$group_by_vars else NULL
      mdcols <- colnames(obj@meta.data)
      fluidRow(
        column(6,
               h4("Select Sample"),
               selectInput("sample_select", "Sample variable", choices = mdcols,
                           selected = if (is_harmony && !is.null(group_by_vars)) group_by_vars else NULL),
               plotOutput("sample_dimplot"),
               uiOutput("sample_subset_ui")
        ),
        column(6,
               h4("Select Cluster"),
               selectInput("cluster_select", "Cluster variable", choices = mdcols,
                           selected = mdcols[1]),
               plotOutput("cluster_dimplot"),
               uiOutput("cluster_subset_ui")
        )
      )
    } else if (step == 3) {
      req(processedObj())
      obj <- processedObj()
      mdcols <- colnames(obj@meta.data)
      reductions <- names(obj@reductions)
      is_harmony <- isTRUE(input$do_harmony)
      default_reduction <- if (is_harmony && "harmony" %in% reductions) "harmony" else "pca"
      fluidRow(
        column(4,
               selectInput("step3_reduction", "Select reduction:", choices = reductions,
                           selected = default_reduction),
               uiOutput("start_cluster_ui"),
               actionButton("run_sling_btn", "RUN SLINGSHOT"),
               br(),
               uiOutput("step3_next_ui")
        ),
        column(8,
               uiOutput("step3_plots_ui")
        )
      )
    } else if (step == 4) {
      req(processedObj(), step3_result())
      fluidRow(
        column(4,
               numericInput("non_zero_num", "Minimum Non-zero Cells:", value = 100, min = 1),
               numericInput("lower_quantile", "Lower Quantile:", min = 0, max = 1, value = 0.001, step = 0.001),
               numericInput("upper_quantile", "Upper Quantile:", min = 0, max = 1, value = 0.999, step = 0.001),
               hr(),
               tags$div(strong("BEFORE FILTERING"), style = "margin-bottom:5px;"),
               tags$div("GENE NUMBER: ", span(textOutput("gene_num_before"), style = "color:red;font-weight:bold;")),
               tags$div("CELL NUMBER: ", span(textOutput("cell_num_before"), style = "color:red;font-weight:bold;")),
               br(),
               actionButton("run_filtering_btn", "RUN FILTERING"),
               hr(),
               tags$div(strong("AFTER FILTERING"), style = "margin-bottom:5px;"),
               tags$div("GENE NUMBER: ", span(textOutput("gene_num_after"), style = "color:red;font-weight:bold;")),
               tags$div("CELL NUMBER: ", span(textOutput("cell_num_after"), style = "color:red;font-weight:bold;")),
               br(),
               actionButton("next_btn", "NEXT STEP", class = "btn btn-primary", disabled = !filter_run_done())
        ),
        column(8,
               plotOutput("pseudotime_boxplot")
        )
      )
    } else if (step == 5) {
      fluidRow(
        column(5,
               h3("Lamian Parameters"),
               textInput("output_path", "Output file path:", value = ""),
               hr(),
               h4("Design Matrix Settings"),
               numericInput("design_ncol", "Number of columns:", 
                            min = 1, max = 30, value = 2),
               actionButton("generate_design", "Generate Design Matrix",  class = "btn btn-primary"),
               hr(),
               numericInput("max_knot", "Maximum knot for spline:", min = 4, value = 10),
               numericInput("perm_iter", "Permutation iterations:", min = 11, value = 100),
               hr(),
               textInput("project_name", "Project name:", value = ""),
               textInput("email", "Email address:", value = ""),
               actionButton("run_lamian", "Run Lamian Analysis", 
                            class = "btn btn-primary")
        ),
        column(7,
               rHandsontableOutput("design_matrix"))
      )
    }
  })
  
  output$harmony_ui <- renderUI({
    obj <- seuratObj()
    if (is.null(obj)) return(tags$em("Please read a valid Seurat object first."))
    md_cols <- tryCatch(colnames(obj@meta.data), error = function(e) NULL)
    if (is.null(md_cols)) return(tags$em("Cannot extract meta.data columns."))
    tagList(
      selectInput("group_by_vars", "Select group.by.vars", choices = md_cols, multiple = FALSE),
      numericInput("dims_use", "Set dims.use", value = 10, min = 1, max = 100)
    )
  })
  
  observeEvent(input$read_btn_ui, {
    path <- input$data_path
    if (!nzchar(path) || !file.exists(path)) {
      showModal(modalDialog(title = "Error", "Please enter a valid existing path!", easyClose = TRUE))
      dataLoaded(FALSE)
      runDoneStep1(FALSE)
      return()
    }
    showModal(modalDialog(title = "Please wait", "Reading...", footer = NULL, easyClose = FALSE))
    tryCatch({
      obj <- readRDS(path)
      seuratObj(obj)
      removeModal()
      showNotification("RDS file loaded.", type = "message")
      dataLoaded(TRUE)
      runDoneStep1(FALSE)
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
      seuratObj(NULL)
      dataLoaded(FALSE)
      runDoneStep1(FALSE)
    })
  })
  
  observeEvent(input$run_btn_ui, {
    if (!dataLoaded()) {
      showModal(modalDialog(title = "Warning", "Please read data first.", easyClose = TRUE))
      return()
    }
    obj <- seuratObj()
    is_harmony <- isTRUE(input$do_harmony)
    group_by_vars <- if (is_harmony) input$group_by_vars else NULL
    dims_use <- if (is_harmony) input$dims_use else NULL
    showModal(modalDialog(title = "Please wait", "Running step1...", footer = NULL, easyClose = FALSE))
    tryCatch({
      obj_processed <- step1(data = obj, is_harmony, group_by_vars, dims_use)
      processedObj(obj_processed)
      data_temp(obj_processed)
      removeModal()
      showNotification("Step 1 completed.", type = "message")
      runDoneStep1(TRUE)
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", paste0("Step1 error:\n", e$message), easyClose = TRUE))
      runDoneStep1(FALSE)
    })
  })
  
  output$sample_subset_ui <- renderUI({
    dt <- processedObj()
    req(dt, input$sample_select)
    vals <- sort(unique(dt@meta.data[[input$sample_select]]))
    selected_vals <- isolate(input$sample_subset_vals)
    if (is.null(selected_vals) || !all(selected_vals %in% vals)) selected_vals <- vals
    checkboxGroupInput("sample_subset_vals", "Subset samples", choices = vals, selected = selected_vals)
  })
  
  output$cluster_subset_ui <- renderUI({
    dt <- processedObj()
    req(dt, input$cluster_select)
    vals <- sort(unique(dt@meta.data[[input$cluster_select]]))
    selected_vals <- isolate(input$cluster_subset_vals)
    if (is.null(selected_vals) || !all(selected_vals %in% vals)) selected_vals <- vals
    checkboxGroupInput("cluster_subset_vals", "Subset clusters", choices = vals, selected = selected_vals)
  })
  
  output$sample_dimplot <- renderPlot({
    dt <- data_temp()
    req(dt, input$sample_select)
    DimPlot(dt, group.by = input$sample_select) + ggtitle(paste("DimPlot by", input$sample_select))
  })
  
  output$cluster_dimplot <- renderPlot({
    dt <- data_temp()
    req(dt, input$cluster_select)
    DimPlot(dt, group.by = input$cluster_select) + ggtitle(paste("DimPlot by", input$cluster_select))
  })
  
  observeEvent(input$run_subset_btn, {
    obj <- processedObj()
    req(obj)
    sample_var <- input$sample_select
    cluster_var <- input$cluster_select
    sample_vals <- input$sample_subset_vals
    cluster_vals <- input$cluster_subset_vals
    showModal(modalDialog(title = "Please wait", "Subsetting data...", footer = NULL, easyClose = FALSE))
    tryCatch({
      keep_sample <- colnames(obj)[obj@meta.data[[sample_var]] %in% sample_vals]
      keep_cluster <- colnames(obj)[obj@meta.data[[cluster_var]] %in% cluster_vals]
      keep_cells <- intersect(keep_sample, keep_cluster)
      obj_sub <- subset(obj, cells = keep_cells)
      data_temp(obj_sub)
      removeModal()
      showNotification("Subset completed.", type = "message")
      subsetDone(TRUE)
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  
  output$start_cluster_ui <- renderUI({
    obj <- processedObj()
    req(obj)
    cls_col <- input$cluster_select
    req(cls_col)
    vals <- unique(obj@meta.data[[cls_col]])
    selectInput("step3_start_cluster", "Select start cluster:", choices = vals, selected = vals[1])
  })
  
  observeEvent(input$run_sling_btn, {
    req(processedObj())
    obj <- processedObj()
    sample <- input$sample_select
    clusters <- input$cluster_select
    reduction <- input$step3_reduction
    start_cluster <- input$step3_start_cluster
    showModal(modalDialog(title = "Please wait", "Running Slingshot...", footer = NULL, easyClose = FALSE))
    tryCatch({
      res <- step3(obj, sample, clusters, reduction, start_cluster)
      step3_result(res)
      removeModal()
      showNotification("Slingshot completed.", type = "message")
      step3_run_done(TRUE)
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
      step3_run_done(FALSE)
    })
  })
  
  output$step3_plots_ui <- renderUI({
    req(step3_result())
    tagList(
      plotOutput("step3_plot_cluster"),
      plotOutput("step3_plot_pseudotime")
    )
  })
  
  output$step3_plot_cluster <- renderPlot({
    req(step3_result())
    step3_result()$plot1
  })
  
  output$step3_plot_pseudotime <- renderPlot({
    req(step3_result())
    step3_result()$plot2
  })
  
  output$step3_next_ui <- renderUI({
    disabled <- !step3_run_done()
    actionButton("next_btn", "Next step",
                 class = if (disabled) "btn btn-secondary" else "btn btn-primary",
                 disabled = disabled)
  })
  
  observeEvent(input$run_filtering_btn, {
    req(processedObj(), step3_result())
    obj <- processedObj()
    sce <- step3_result()$sce
    non_zero_num <- input$non_zero_num
    lower_quantile <- input$lower_quantile
    upper_quantile <- input$upper_quantile
    if (upper_quantile <= lower_quantile) {
      showModal(modalDialog(title = "Error", "Upper quantile must be greater than lower quantile.", easyClose = TRUE))
      return()
    }
    showModal(modalDialog(title = "Please wait", "Filtering...", footer = NULL, easyClose = FALSE))
    tryCatch({
      res <- step4(obj, sce, non_zero_num, lower_quantile, upper_quantile)
      step4_result(res)
      filter_run_done(TRUE)
      obj_filtered_temp(res$obj)
      
      updateNumericInput(session, "lower_quantile", value = res$lower_quantile)
      updateNumericInput(session, "upper_quantile", value = res$upper_quantile)
      
      removeModal()
      showNotification("Filtering completed.", type = "message")
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", paste0("Filtering error:\n", e$message), easyClose = TRUE))
      filter_run_done(FALSE)
    })
  })
  
  # design matrix
  design_matrix <- reactiveVal(NULL)
  observeEvent(input$generate_design, {
    n_data_cols <- input$design_ncol
    if (is.null(n_data_cols) || n_data_cols < 1) {
      showNotification("Number of data columns must be ≥1", type = "error")
      return()
    }
    
    # Get sample names from the processed data
    sample_names <- unique(processedObj()@meta.data[[input$sample_select]])
    
    nr <- length(sample_names)
    # Create data frame with Sample column
    df <- data.frame(
      Sample = sample_names,
      stringsAsFactors = FALSE
    )
    
    # Add data columns with specific naming convention
    if (n_data_cols >= 1) {
      df[["test.variable"]] <- rep(FALSE, nr)  # First column after Sample
      
      if (n_data_cols > 1) {
        # Add confounder columns
        for (i in 1:(n_data_cols-1)) {
          df[[paste0("confounder", i)]] <- rep(FALSE, nr)
        }
      }
    }
    design_matrix(df)
  })
  
  output$removed_cells_text <- renderText({
    req(filter_run_done())
    res <- step4_result()
    paste("Number of removed cells:", res$remove_cells_num)
  })
  
  output$gene_num_before <- renderText({
    obj <- processedObj()
    req(obj)
    expr <- GetAssayData(obj, layer = "data")
    nrow(expr)
  })
  
  output$cell_num_before <- renderText({
    obj <- processedObj()
    req(obj)
    expr <- GetAssayData(obj, layer = "data")
    ncol(expr)
  })
  
  output$gene_num_after <- renderText({
    obj <- obj_filtered_temp()
    if (is.null(obj)) return("0")
    expr <- GetAssayData(obj, layer = "data")
    nrow(expr)
  })
  
  output$cell_num_after <- renderText({
    obj <- obj_filtered_temp()
    if (is.null(obj)) return("0")
    expr <- GetAssayData(obj, layer = "data")
    ncol(expr)
  })
  
  
  output$pseudotime_boxplot <- renderPlot({
    sce <- step3_result()$sce
    pseudotime_matrix <- assay(sce@colData$slingshot, "pseudotime")
    pseudotime_vec <- pseudotime_matrix[,1]
    
    sample_id_vec <- sce@colData$sample
    
    df <- data.frame(
      pseudotime = pseudotime_vec,
      sample_id = sample_id_vec,
      stringsAsFactors = FALSE
    )
    
    df_all <- df %>%
      mutate(sample_id = "All_samples")
    df_plot <- bind_rows(df, df_all)
    Q1_all <- quantile(df$pseudotime, probs = input$lower_quantile, na.rm = TRUE)
    Q3_all <- quantile(df$pseudotime, probs = input$upper_quantile, na.rm = TRUE)
    
    x_min <- 0
    x_max <- max(df$pseudotime, na.rm = TRUE)
    
    ggplot(df_plot, aes(x = pseudotime, y = sample_id)) +
      geom_boxplot(outlier.shape = 1, fill = "lightblue") + 
      geom_vline(xintercept = Q1_all, color = "red", linetype = "dashed") +
      geom_vline(xintercept = Q3_all, color = "red", linetype = "dashed") +
      scale_x_continuous(limits = c(x_min, x_max)) +
      xlab("Pseudotime") +
      ylab("Sample ID") +
      theme_bw(base_size = 14) +  # increase base font size here
      theme(
        axis.text.x = element_text(size = 15),  # customize x-axis text size
        axis.text.y = element_text(size = 15)   # y-axis text size
      )
  })
  
  output$design_matrix <- renderRHandsontable({
    df <- design_matrix()
    req(df)
    
    rh <- rhandsontable(df, useTypes = TRUE, stretchH = "all") %>%
      hot_col("Sample", readOnly = TRUE)  # Make Sample column uneditable
    
    if (ncol(df) > 1) {
      # Set all columns except Sample as checkboxes
      for (colname in colnames(df)[-1]) {
        rh <- hot_col(rh, colname, type = "checkbox")
      }
    }
    rh
  })
  
  observeEvent(input$design_matrix, {
    df <- hot_to_r(input$design_matrix)
    design_matrix(df)
  })
  
  observeEvent(input$run_lamian, {
    req(processedObj(), step3_result(), design_matrix())
    
    # Get all parameters
    output_file_path <- input$output_path
    design <- hot_to_r(input$design_matrix)
    maximumknotallowed <- input$max_knot
    permuiter <- input$perm_iter
    project_name <- input$project_name
    email_address <- input$email
    
    # Validate inputs
    if (nchar(project_name) == 0) {
      showNotification("Please enter a project name", type = "error")
      return()
    }
    
    if (nchar(email_address) == 0 || !grepl("@", email_address)) {
      showNotification("Please enter a valid email address", type = "error")
      return()
    }
    
    showModal(modalDialog(
      title = "Running Lamian Analysis",
      "This may take several minutes...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      design_conv <- design
      if (ncol(design_conv) > 1) {
        for (i in 2:ncol(design_conv)) {
          if (is.logical(design_conv[[i]])) {
            design_conv[[i]] <- as.integer(design_conv[[i]])
          }
        }
      }
      
      step5_result <- step5(
        data = processedObj(),
        sce = step3_result()$sce,
        output_file_path = output_file_path,
        design = design_conv,
        maximumknotallowed = maximumknotallowed,
        permuiter = permuiter,
        project_name = project_name,
        email_address = email_address
      )
      removeModal() 
      
      showNotification("Lamian analysis completed.", type = "message")
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", paste0("Lamian error:\n", e$message), easyClose = TRUE))
    })
  })
  
  output$nav_buttons <- renderUI({
    step <- current_step()
    if (step == 1) {
      fluidRow(
        column(4,
               actionButton("run_btn_ui", "RUN",
                            class = if (dataLoaded()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !dataLoaded())
        ),
        column(4,
               actionButton("next_btn", "Next step",
                            class = if (runDoneStep1()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !runDoneStep1())
        )
      )
    } else if (step == 2) {
      fluidRow(
        column(3,
               actionButton("run_subset_btn", "RUN Subset",
                            class = if (runDoneStep1()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !runDoneStep1())
        ),
        column(9,
               actionButton("next_btn", "Next step",
                            class = "btn btn-primary",
                            disabled = FALSE)
        )
      )
    } else if (step == 3) {
      NULL
    } else if (step == 4) {
      NULL
    } else if (step == 5) {
      NULL
    }
  })
  
  output$step3_next_ui <- renderUI({
    disabled <- !step3_run_done()
    actionButton("next_btn", "Next step",
                 class = if (disabled) "btn btn-secondary" else "btn btn-primary",
                 disabled = disabled)
  })
  
  observeEvent(input$next_btn, {
    step <- current_step()
    if (step == 1) {
      req(runDoneStep1())
      current_step(step + 1)
    } else if (step == 2) {
      if (subsetDone()) processedObj(data_temp())
      current_step(step + 1)
    } else if (step == 3) {
      req(step3_run_done())
      current_step(step + 1)
    } else if (step == 4) {
      req(filter_run_done())
      res <- step4_result()
      if (!is.null(res)) processedObj(res$obj)
      current_step(step + 1)
    }
  })
  
  observeEvent(input$start_lamian_btn, {
    showModal(modalDialog(
      title = "Notice",
      "Lamian starting... (please implement your logic here)",
      easyClose = TRUE
    ))
  })
}

shinyApp(ui, server)