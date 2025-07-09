library(shiny)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(slingshot)
library(harmony)

source("step1.R")  # Your step1() function
source("step3.R")  # Your step3() function
source("step4.R")  # Your step4() function

ui <- fluidPage(
  titlePanel("Multi-step Preprocessing & Slingshot Analysis"),
  tabsetPanel(
    tabPanel("Data preprocessing",
             fluidRow(
               column(width = 8,
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
  total_steps <- 6
  
  # Variables for steps
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
  
  # UI step header
  output$step_ui <- renderUI({
    h4(paste0("Step ", current_step(), " of ", total_steps))
  })
  
  output$step_content <- renderUI({
    step <- current_step()
    if(step == 1) {
      fluidPage(
        fluidRow(
          column(8, textInput("data_path", "Enter path to RDS file", value = "/projectnb/wax-es/ForVallari/G183_G193.rds")),
          column(4, br(), actionButton("read_btn_ui", "READ"))
        ),
        checkboxInput("do_harmony", "Perform Harmony correction", FALSE),
        conditionalPanel(
          condition="input.do_harmony == true",
          uiOutput("harmony_ui")
        )
      )
    } else if(step == 2) {
      req(processedObj())
      obj <- processedObj()
      is_harmony <- isTRUE(input$do_harmony)
      group_by_vars <- if(is_harmony) input$group_by_vars else NULL
      mdcols <- colnames(obj@meta.data)
      fluidRow(
        column(6,
               h4("Select Sample"),
               selectInput("sample_select", "Sample variable", choices = mdcols,
                           selected = if(is_harmony && !is.null(group_by_vars)) group_by_vars else NULL),
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
    } else if(step == 3) {
      req(processedObj())
      obj <- processedObj()
      mdcols <- colnames(obj@meta.data)
      reductions <- names(obj@reductions)
      is_harmony <- isTRUE(input$do_harmony)
      default_reduction <- if(is_harmony && "harmony" %in% reductions) "harmony" else "pca"
      fluidRow(
        column(4,
               selectInput("step3_clusters", "Select clusters:", choices = mdcols,
                           selected = if(!is.null(input$cluster_select)) input$cluster_select else mdcols[1]),
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
    } else if(step == 4) {
      req(processedObj(), step3_result())
      fluidRow(
        column(4,
               numericInput("non_zero_num", "Minimum Non-zero Cells:", value = 100, min = 1),
               sliderInput("lower_quantile", "Lower Quantile:", 0, 1, 0.1),
               sliderInput("upper_quantile", "Upper Quantile:", 0, 1, 0.9),
               numericInput("IQR_coefficient", "IQR Coefficient:", value = 1.5, min = 0),
               actionButton("run_filtering_btn", "RUN FILTERING"),
               br(),
               verbatimTextOutput("removed_cells_text")
        ),
        column(8,
               plotOutput("gene_counts_barplot"),
               plotOutput("pseudotime_boxplot")
        )
      )
    } else {
      h3(paste("Step", step, "- To be implemented"))
    }
  })
  
  # Step1 Harmony UI
  output$harmony_ui <- renderUI({
    obj <- seuratObj()
    if(is.null(obj)) return(tags$em("Please read a valid Seurat object first."))
    md_cols <- tryCatch(colnames(obj@meta.data), error = function(e) NULL)
    if(is.null(md_cols)) return(tags$em("Cannot extract meta.data columns."))
    tagList(
      selectInput("group_by_vars", "Select group.by.vars", choices = md_cols, multiple = FALSE),
      numericInput("dims_use", "Set dims.use", value = 10, min = 1, max = 100)
    )
  })
  
  # Step1 READ button
  observeEvent(input$read_btn_ui, {
    path <- input$data_path
    if(!nzchar(path) || !file.exists(path)) {
      showModal(modalDialog(title="Error", "Please enter a valid existing path!", easyClose=TRUE))
      dataLoaded(FALSE)
      runDoneStep1(FALSE)
      return()
    }
    showModal(modalDialog(title="Please wait", "Reading...", footer=NULL, easyClose=FALSE))
    tryCatch({
      obj <- readRDS(path)
      seuratObj(obj)
      removeModal()
      showNotification("RDS file loaded.", type="message")
      dataLoaded(TRUE)
      runDoneStep1(FALSE)
    }, error = function(e){
      removeModal()
      showModal(modalDialog(title="Error", e$message, easyClose=TRUE))
      seuratObj(NULL); dataLoaded(FALSE); runDoneStep1(FALSE)
    })
  })
  
  # Step1 RUN button
  observeEvent(input$run_btn_ui, {
    if(!dataLoaded()) {
      showModal(modalDialog(title="Warning","Please READ data first.",easyClose=TRUE))
      return()
    }
    obj <- seuratObj()
    is_harmony <- isTRUE(input$do_harmony)
    group_by_vars <- if(is_harmony) input$group_by_vars else NULL
    dims_use <- if(is_harmony) input$dims_use else NULL
    showModal(modalDialog(title="Please wait","Running step1...",footer=NULL,easyClose=FALSE))
    tryCatch({
      obj_processed <- step1(data=obj, is_harmony, group_by_vars, dims_use)
      processedObj(obj_processed)
      data_temp(obj_processed)
      removeModal()
      showNotification("Step 1 completed.", type="message")
      runDoneStep1(TRUE)
    }, error=function(e){
      removeModal()
      showModal(modalDialog(title="Error", e$message, easyClose=TRUE))
      runDoneStep1(FALSE)
    })
  })
  
  # Step2 subset UI and plots
  output$sample_subset_ui <- renderUI({
    dt <- data_temp()
    req(dt, input$sample_select)
    vals <- sort(unique(dt@meta.data[[input$sample_select]]))
    selected_vals <- isolate(input$sample_subset_vals)
    if(is.null(selected_vals) || !all(selected_vals %in% vals)) selected_vals <- vals
    checkboxGroupInput("sample_subset_vals", "Subset samples", choices=vals, selected=selected_vals)
  })
  
  output$cluster_subset_ui <- renderUI({
    dt <- data_temp()
    req(dt, input$cluster_select)
    vals <- sort(unique(dt@meta.data[[input$cluster_select]]))
    selected_vals <- isolate(input$cluster_subset_vals)
    if(is.null(selected_vals) || !all(selected_vals %in% vals)) selected_vals <- vals
    checkboxGroupInput("cluster_subset_vals", "Subset clusters", choices=vals, selected=selected_vals)
  })
  
  output$sample_dimplot <- renderPlot({
    dt <- data_temp()
    req(dt, input$sample_select)
    DimPlot(dt, group.by=input$sample_select) + ggtitle(paste("DimPlot by", input$sample_select))
  })
  
  output$cluster_dimplot <- renderPlot({
    dt <- data_temp()
    req(dt, input$cluster_select)
    DimPlot(dt, group.by=input$cluster_select) + ggtitle(paste("DimPlot by", input$cluster_select))
  })
  
  observeEvent(input$run_subset_btn, {
    obj <- processedObj()
    req(obj)
    sample_var <- input$sample_select
    cluster_var <- input$cluster_select
    sample_vals <- input$sample_subset_vals
    cluster_vals <- input$cluster_subset_vals
    showModal(modalDialog(title="Please wait", "Subsetting data...", footer=NULL, easyClose=FALSE))
    tryCatch({
      keep_sample <- colnames(obj)[obj@meta.data[[sample_var]] %in% sample_vals]
      keep_cluster <- colnames(obj)[obj@meta.data[[cluster_var]] %in% cluster_vals]
      keep_cells <- intersect(keep_sample, keep_cluster)
      obj_sub <- subset(obj, cells=keep_cells)
      data_temp(obj_sub)
      removeModal()
      showNotification("Subset completed.", type="message")
      subsetDone(TRUE)
    }, error=function(e){
      removeModal()
      showModal(modalDialog(title="Error", e$message, easyClose=TRUE))
    })
  })
  
  # Step 3 UI components
  output$start_cluster_ui <- renderUI({
    obj <- processedObj()
    req(obj)
    cls_col <- input$step3_clusters
    req(cls_col)
    vals <- unique(obj@meta.data[[cls_col]])
    selectInput("step3_start_cluster", "Select start cluster:", choices=vals, selected=vals[1])
  })
  
  observeEvent(input$run_sling_btn, {
    req(processedObj())
    obj <- processedObj()
    clusters <- input$step3_clusters
    reduction <- input$step3_reduction
    start_cluster <- input$step3_start_cluster
    showModal(modalDialog(title="Please wait","Running Slingshot...",footer=NULL,easyClose=FALSE))
    tryCatch({
      res <- step3(obj, clusters, reduction, start_cluster)
      step3_result(res)
      removeModal()
      showNotification("Slingshot completed.", type="message")
      step3_run_done(TRUE)
    }, error=function(e){
      removeModal()
      showModal(modalDialog(title="Error", e$message, easyClose=TRUE))
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
                 class = if(disabled) "btn btn-secondary" else "btn btn-primary",
                 disabled = disabled)
  })
  
  # Step 4 observe for filtering run
  observeEvent(input$run_filtering_btn, {
    req(processedObj(), step3_result())
    obj <- processedObj()
    sce <- step3_result()$sce
    non_zero_num <- input$non_zero_num
    lower_quantile <- input$lower_quantile
    upper_quantile <- input$upper_quantile
    IQR_coefficient <- input$IQR_coefficient
    if(upper_quantile <= lower_quantile) {
      showModal(modalDialog(title="Error", "Upper quantile must be greater than lower quantile.", easyClose=TRUE))
      return()
    }
    showModal(modalDialog(title="Please wait", "Filtering...", footer=NULL, easyClose=FALSE))
    tryCatch({
      res <- step4(obj, sce, non_zero_num, lower_quantile, upper_quantile, IQR_coefficient)
      step4_result(res)
      filter_run_done(TRUE)
      removeModal()
      showNotification("Filtering completed.", type="message")
    }, error=function(e){
      removeModal()
      showModal(modalDialog(title="Error", e$message, easyClose=TRUE))
      filter_run_done(FALSE)
    })
  })
  
  output$removed_cells_text <- renderText({
    req(filter_run_done())
    res <- step4_result()
    paste("Number of removed cells:", res$remove_cells_num)
  })
  
  output$gene_counts_barplot <- renderPlot({
    req(step4_result())
    res <- step4_result()
    obj_filtered <- res$obj
    expr <- GetAssayData(obj_filtered, layer="data")
    expr_count <- rowSums(expr > 0)
    cnt <- table(expr_count)
    df <- as.data.frame(cnt)
    colnames(df) <- c("NonZeroCells", "GeneCount")
    df$NonZeroCells <- as.numeric(as.character(df$NonZeroCells))
    ggplot(df, aes(x=NonZeroCells,y=GeneCount)) +
      geom_bar(stat="identity", fill="steelblue") +
      labs(title="Gene counts by non-zero cell number",
           x="Non-zero cells", y="Number of genes") +
      theme_minimal()
  })
  
  output$pseudotime_boxplot <- renderPlot({
    req(step4_result())
    res <- step4_result()
    obj_filtered <- res$obj
    sce <- step3_result()$sce
    pseudotime <- slingPseudotime(sce, na=TRUE)[,1]
    names(pseudotime) <- colnames(sce)
    pt_filtered <- pseudotime[colnames(obj_filtered)]
    sample_colname <- "samples" # adjust if needed
    samples <- obj_filtered@meta.data[,sample_colname]
    df <- data.frame(Sample=samples, Pseudotime=pt_filtered)
    Q1 <- res$Q1
    Q3 <- res$Q3
    IQR <- res$IQR
    lower_cutoff <- Q1 - input$IQR_coefficient * IQR
    upper_cutoff <- Q3 + input$IQR_coefficient * IQR
    p <- ggplot(df, aes(x=Sample,y=Pseudotime,fill=Sample)) +
      geom_boxplot() +
      geom_hline(yintercept=lower_cutoff, color="red", linetype="dashed", size=1) +
      geom_hline(yintercept=upper_cutoff, color="red", linetype="dashed", size=1) +
      coord_flip() +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title="Pseudotime distribution by sample")
    p_all <- ggplot(df, aes(x=factor(1), y=Pseudotime)) +
      geom_boxplot(fill="grey70") +
      geom_hline(yintercept=lower_cutoff, color="red", linetype="dashed", size=1) +
      geom_hline(yintercept=upper_cutoff, color="red", linetype="dashed", size=1) +
      coord_flip() +
      theme_minimal() +
      labs(title="Overall pseudotime distribution", x="", y="Pseudotime") +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
    if(requireNamespace("patchwork", quietly = TRUE)) {
      patchwork::wrap_plots(p, p_all, ncol = 1, heights = c(3,1))
    } else if(requireNamespace("cowplot", quietly = TRUE)) {
      cowplot::plot_grid(p, p_all, ncol = 1, rel_heights = c(3,1))
    } else {
      gridExtra::grid.arrange(p, p_all, ncol = 1, heights = c(3,1))
    }
  })
  
  # Navigation buttons
  output$nav_buttons <- renderUI({
    step <- current_step()
    if(step == 1) {
      fluidRow(
        column(4,
               actionButton("run_btn_ui", "RUN",
                            class = if(dataLoaded()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !dataLoaded())
        ),
        column(4,
               actionButton("next_btn", "Next step",
                            class = if(runDoneStep1()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !runDoneStep1())
        )
      )
    } else if(step == 2) {
      fluidRow(
        column(3,
               actionButton("run_subset_btn", "RUN Subset",
                            class = if(runDoneStep1()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !runDoneStep1())
        ),
        column(9,
               actionButton("next_btn", "Next step",
                            class = "btn btn-primary",
                            disabled = FALSE)
        )
      )
    } else if(step == 3) {
      NULL # Step 3 next handled separately in UI
    } else if(step == 4) {
      fluidRow(
        column(4,
               actionButton("run_filtering_btn", "RUN FILTERING",
                            class = if(runDoneStep1()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !runDoneStep1())
        ),
        column(8,
               actionButton("next_btn", "Next step",
                            class = if(filter_run_done()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !filter_run_done())
        )
      )
    } else if(step == total_steps) {
      fluidRow(
        column(12,
               actionButton("start_lamian_btn", "Start Lamian")
        )
      )
    } else {
      fluidRow(
        column(12,
               actionButton("next_btn", "Next step")
        )
      )
    }
  })
  
  # Step 3 Next step button under RUN button
  output$step3_next_ui <- renderUI({
    disabled <- !step3_run_done()
    actionButton("next_btn", "Next step",
                 class = if(disabled) "btn btn-secondary" else "btn btn-primary",
                 disabled = disabled)
  })
  
  # Next step logic
  observeEvent(input$next_btn, {
    step <- current_step()
    if(step == 1) {
      req(runDoneStep1())
      current_step(step + 1)
    } else if(step == 2) {
      if(subsetDone()) processedObj(data_temp())
      current_step(step + 1)
    } else if(step == 3) {
      req(step3_run_done())
      current_step(step + 1)
    } else if(step == 4) {
      req(filter_run_done())
      res <- step4_result()
      if(!is.null(res)) processedObj(res$obj)
      current_step(step + 1)
    } else if(step < total_steps) {
      current_step(step + 1)
    }
  })
  
  # Start Lamian button
  observeEvent(input$start_lamian_btn, {
    showModal(modalDialog(
      title = "Notice",
      "Lamian starting... (please implement your logic here)",
      easyClose = TRUE
    ))
  })
  
}

shinyApp(ui, server)