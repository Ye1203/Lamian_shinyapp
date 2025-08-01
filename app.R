library(shiny)
library(shinyjs)
library(DT)
library(Seurat)
library(ggplot2)
library(scater)
library(viridis)
library(SingleCellExperiment)
library(slingshot)
library(harmony)
library(dplyr)
library(tibble)
library(tidyr)
library(rhandsontable)
library(plotly)
devtools::load_all("/projectnb/wax-es/Bingtian/Lamian/renv/library/linux-almalinux-8.10/R-4.4/x86_64-pc-linux-gnu/Lamian")
source("/projectnb/wax-es/00_shinyapp/Lamian/lamian/step1.R")
source("/projectnb/wax-es/00_shinyapp/Lamian/lamian/step3.R")
source("/projectnb/wax-es/00_shinyapp/Lamian/lamian/step4.R")
source("/projectnb/wax-es/00_shinyapp/Lamian/lamian/step5.R")
source("/projectnb/wax-es/00_shinyapp/Lamian/lamian/visualization.R")
source("/projectnb/wax-es/00_shinyapp/Lamian/lamian/visualization_read.R")
source("/projectnb/wax-es/00_shinyapp/Lamian/lamian/draw_heatmap.R")
set.seed(123)
ui <- fluidPage(
  titlePanel("Lamian - Waxman's lab"),
  tags$h4("Made by Bingtian Ye(btye@bu.edu)"),
  tabsetPanel(
    tabPanel("Data preprocessing",
             fluidRow(
               column(12,
                      uiOutput("step_ui"),
                      uiOutput("step_content")
               )
             )
    ),
    tabPanel("Data visualization",
             dataVisualizationUI()
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
    step <- current_step()
    if (step == 1) {
      h4("Step 1 of 5: Reading data, Normalizing, Harmonizing")
    }else if (step == 2){
      h4("Step 2 of 5: Selecting sample and cell cluster index, Subseting")
    }else if (step == 3){
      h4("Step 3 of 5: Trajectory analysis - Slingshot")
    }else if (step == 4){
      h4("Step 4 of 5: Filtering")
    }else if (step == 5){
      h4("Step 5 of 5: Lamian parameter setting")
    }
  })
  
  output$step_content <- renderUI({
    step <- current_step()
    if (step == 1) {
      fluidPage(
        fluidRow(
          column(4,
                 fluidRow(
                   column(8, textInput("data_path", "Enter path to RDS (Seurat object) file", value = "")),
                   column(4, br(), actionButton("read_btn_ui", "READ"))
                 ),
                 checkboxInput("do_harmony", "Perform Harmony correction", TRUE),
                 conditionalPanel(
                   condition = "input.do_harmony == true",
                   uiOutput("harmony_ui")
                 ),
                 uiOutput("nav_buttons")
          ),
          column(4,
                 plotOutput("step1_plot")
          ),
          column(4,
                 uiOutput("step1_description")
          )
        )
      )
    } else if (step == 2) {
      req(processedObj())
      obj <- processedObj()
      is_harmony <- isTRUE(input$do_harmony)
      group_by_vars <- if (is_harmony) input$group_by_vars else NULL
      mdcols <- names(sapply(obj@meta.data, function(x) !is.numeric(x))[sapply(obj@meta.data, function(x) !is.numeric(x))])
      if (!is.null(mdcols) && "CB" %in% mdcols) {
        mdcols <- setdiff(mdcols, "CB")
      }
      tagList(
        fluidRow(
          column(4,
                 h4("Select Sample index"),
                 selectInput("sample_select", "", choices = mdcols,
                             selected = if (is_harmony && !is.null(group_by_vars)) group_by_vars else NULL),
                 plotOutput("sample_dimplot"),
                 uiOutput("sample_subset_ui")
          ),
          column(4,
                 h4("Select Cell Cluster index"),
                 selectInput("cluster_select", "", choices = mdcols,
                             selected = mdcols[1]),
                 plotOutput("cluster_dimplot"),
                 uiOutput("cluster_subset_ui")
          ),
          column(4,
                 uiOutput("step2_description")
          )
        ),
        fluidRow(
          column(12,
                 uiOutput("nav_buttons")
          )
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
               selectInput("step3_reduction", "Select dimensionality reduction method:", choices = reductions,
                           selected = default_reduction),
               uiOutput("start_cluster_ui"),
               actionButton("run_sling_btn", "INITIATE SLINGSHOT"),
               br(),
               uiOutput("lineage_select_ui"),
               br(),
               uiOutput("step3_next_ui")
        ),
        column(4,
               uiOutput("step3_plots_ui")
        ),
        column(4,
               uiOutput("step3_description")
        )
      )
    } else if (step == 4) {
      req(processedObj(), step3_result())
      fluidRow(
        column(4,
               numericInput("non_zero_num", "Gene Filter: Minimum of Non-zero Cells", value = 100, min = 1),
               numericInput("lower_quantile", "Cell Filter: Lower Pseudotime Quantile (%)", min = 0, max = 100, value = 0.01, step = 0.01),
               numericInput("upper_quantile", "Cell Filter: Upper Pseudotime Quantile (%)", min = 0, max = 100, value = 99.99, step = 0.01),
               hr(),
               tags$div(strong("BEFORE FILTERING"), style = "margin-bottom:5px;"),
               tags$div("GENE NUMBER: ", span(textOutput("gene_num_before"), style = "color:red;font-weight:bold;")),
               tags$div("CELL NUMBER: ", span(textOutput("cell_num_before"), style = "color:red;font-weight:bold;")),
               br(),
               actionButton("run_filtering_btn", "FILTERING"),
               hr(),
               tags$div(strong("AFTER FILTERING"), style = "margin-bottom:5px;"),
               tags$div("GENE NUMBER: ", span(textOutput("gene_num_after"), style = "color:red;font-weight:bold;")),
               tags$div("CELL NUMBER: ", span(textOutput("cell_num_after"), style = "color:red;font-weight:bold;")),
               br(),
               actionButton("next_btn", "NEXT STEP", class = "btn btn-primary", disabled = !filter_run_done())
        ),
        column(4,
               plotOutput("pseudotime_boxplot")
        ),
        column(4,
               uiOutput("step4_description")
        )
      )
    } else if (step == 5) {
      fluidRow(
        column(3,
               h3("Lamian Parameters"),
               textInput("output_path", "Output file path:", value = ""),
               hr(),
               h4("Design Matrix Settings"),
               numericInput("design_ncol", "Number of columns:", 
                            min = 1, max = 30, value = 2),
               actionButton("generate_design", "Generate Design Matrix",  class = "btn btn-primary"),
               hr(),
               numericInput("max_knot", "Maximum knot for spline:", min = 4, value = 10),
               numericInput("perm_iter", "Permutation iterations:", min = 11, value = 110),
               hr(),
               textInput("project_name", "Project name:", value = ""),
               textInput("email", "Email address:", value = ""),
               actionButton("run_lamian", "Run Lamian Analysis", 
                            class = "btn btn-primary")
        ),
        column(5,
               rHandsontableOutput("design_matrix")),
        column(4,
               uiOutput("step5_description"))
      )
    }
  })
  
  output$harmony_ui <- renderUI({
    obj <- seuratObj()
    if (is.null(obj)) return(tags$em("Please read a valid Seurat object first."))
    md_cols <- tryCatch(names(sapply(obj@meta.data, function(x) !is.numeric(x))[sapply(obj@meta.data, function(x) !is.numeric(x))]), error = function(e) NULL)
    if (!is.null(md_cols) && "CB" %in% md_cols) {
      md_cols <- setdiff(md_cols, "CB")
    }
    if (is.null(md_cols)) return(tags$em("Cannot extract meta.data columns."))
    tagList(
      selectInput("group_by_vars", "Select sample index", choices = md_cols, multiple = FALSE),
      numericInput("dims_use", "Set the number of pca to use", value = length(obj@reductions$pca), min = 1, max = length(obj@reductions$pca))
    )
  })
  
  output$step1_description <- renderUI({
    tags$div(
      style = "background: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #ddd;",
      h4("1.Read data"),
      p("a. Input your file path to your RDS file. (e.g. /projectnb/wax-es/ForVallari/G183_G193.rds)"),
      p("b. Click \"READ\" to read your RDS file."),
      tags$hr(),
      h4("2.Normalize and Harmonize"),
      p("a. After reading the data, if you need to perform harmony (highly recommended), make sure \"Perform Harmony correction\" is checked, otherwise uncheck it."),
      p("b. If you \"Perform Harmony correction\" is checked, select the sample index of the RDS file and select the number of pca dimension to use (default is recommended). "),
      p("c. Click \"INITIATE Normalization\" to start normalization and harmonization."),
      tags$hr(),
      h4("3. Next Step"),
      p("Click \"NEXT STEP\" to go to the step 2."),
      em("")
    )
  })
  output$step2_description <- renderUI({
    tags$div(
      style = "background: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #ddd;",
      h4("1. Select the Sample index and Cell cluster index"),
      p("Select the appropriate index in Select Sample index and Select Cell Cluster. Sample index is generally the sample number and information used to distinguish different samples; Cell Cluster refers to different cell types in biology. After selection, the UMAP plot in the umap dimension will be automatically updated."),
      tags$hr(),
      h4("2. Subset"),
      p("a. Uncheck the samples or cell you don't want to remove these cells from subsequent analysis. Removal can significantly increase the speed of subsequent runs and reduce the interference of unwanted cells on the results."),
      p(style = "color: red;", "The nunber of sample must greater than or equal to 4 and the number of cell cluster must be greater than or equal to 2!!!!"),
      p("b. Click \"INITIATE SUBSET\" to view the results of subset. DI"),
      tags$hr(),
      h4("3. Next Step"),
      p("After confirmation the susbet result, click \"NEXT STEP\" to go to the step 3."),
      em("")
    )
  })
  
  output$step3_description <- renderUI({
    tags$div(
      style = "background: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #ddd;",
      h4("1. Select dimensionality reduction method"),
      p("Select dimensionality reduction method, It is strongly recommended to use harmony."),
      tags$hr(),
      h4("2. Select the start cell cluster and start slingshot analysis"),
      p("a. Select the starting cell type, Slingshot will automatically calculate the starting cell, the pseudotime of this cell will be 0."),
      p("b. Click \"INITIATE SLINGSHOT\" to start slingshot analysis."),
      p("After Slingshot finishes, two plots appear: the top shows cell type distribution, and the bottom shows pseudotime distribution in the chosen dimensionality reduction space."),
      tags$hr(),
      h4("3. Select lineage"),
      p("If there are multiple cell types, multiple trajectory(Lineage) may appear. The currently selected lineage will be marked in red in the figure, and the unselected ones will be marked in black. If the pseudotime is NA (gray), the cell will be automatically removed."),
      tags$hr(),
      h4("4. Next Step"),
      p("After confirming the selected lineage, click \"NEXT STEP\" to go to the step 4."),
      em("")
    )
  })
  output$step4_description <- renderUI({
    tags$div(
      style = "background: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #ddd;",
      h4("1. Set Gene filtering criteria of munimum non-zero cells"),
      p("Calculate the number of cells with non-zero expression for each gene in all selected samples, and remove genes with low-level set values."),
      tags$hr(),
      h4("2. Set Cell filtering criteria of pseudotime"),
      p("Set the cells to be removed based on the upper and lower quantiles of the distribution. The middle boxplot will automatically update the red dotted line based on the selected quantile."),
      tags$hr(),
      h4("3. Filtering"),
      p("Click \"FILTERING\" to filtering cell and gene based on criteria setting."),
      tags$hr(),
      h4("4. Next Step"),
      p("After confirming, click \"NEXT STEP\" to go to the step 5."),
      em("")
    )
  })
  
  output$step5_description <- renderUI({
    tags$div(
      style = "background: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #ddd;",
      h4("1. Input the output file path"),
      p("Enter an absolute path address on scc, and all generated files will be under this file."),
      tags$hr(),
      h4("2. Set design matrix"),
      p("a. Set the number of columns required in the design matrix."),
      p("b. Click \"Generate Design Matrix\" to get the design matrix."),
      p("c. Check the design matrix. Test.variable is the covariate to test; confounders are other factors. For example, if the sample includes gender and injection, focus on gender (test.variable) and treat injection as confounder. Check for yes (1), uncheck for no (0)."),
      tags$hr(),
      h4("3. Other value setting"),
      p("a. \"Maximum knot of spline\" controls spline flexibility; default is recommended."),
      p("b. \"Permutation iterations\" set simulation runs; higher gives more accuracy but takes longer. Default is recommended. Since qsub limits to 24 hours, too large values may cause timeout. If so, contact Bingtian."),
      p("c. \"Project name\" is the project name you can use on scc, such as wax-es, wax-dk, etc."),
      p("d. \"Email address\" is empty by default. Since Lamian runs for hours, provide your email to get results. If no email after 48 hours, it may have timed out. Try reducing cells, genes, or permutation times."),
      tags$hr(),
      h4("4. Start lamian analysis"),
      p("click \"Run Lamian Analysis\" to start lamian data preparing. Please close Shinyapp after the Lamian Job Submitted pop-up window appears."),
      em("")
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
    showModal(modalDialog(title = "Please wait", 
                          tagList(
                            p("Reading..."),
                            p("The waiting time is related to the size of the RDS file and usually takes one minute.")
                          ), footer = NULL, easyClose = FALSE))
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
    if (is_harmony) {showModal(modalDialog(title = "Please wait",
                                          tagList(
                                            p("Running Normalization and Harmonization..."),
                                            p("The waiting time is related to the size of data and the number of dimention. Usually takes several minutes.")
                                          ), footer = NULL, easyClose = FALSE))
      }else if(!is_harmony){showModal(modalDialog(title = "Please wait",
                                                  tagList(
                                                    p("Running Normalization..."),
                                                    p("Usually takes less than 1 minute.")
                                                  ), footer = NULL, easyClose = FALSE))}
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
    checkboxGroupInput("sample_subset_vals", "Subset Samples", choices = vals, selected = selected_vals)
  })
  
  output$cluster_subset_ui <- renderUI({
    dt <- processedObj()
    req(dt, input$cluster_select)
    vals <- sort(unique(dt@meta.data[[input$cluster_select]]))
    selected_vals <- isolate(input$cluster_subset_vals)
    if (is.null(selected_vals) || !all(selected_vals %in% vals)) selected_vals <- vals
    checkboxGroupInput("cluster_subset_vals", "Subset Cell clusters", choices = vals, selected = selected_vals)
  })
  
  output$step1_plot <- renderPlot({
    req(dataLoaded())
    req(input$do_harmony == TRUE)
    data <- seuratObj()
    DimPlot(data, group.by = input$group_by_vars) + ggtitle(paste("UMAP Plot by", input$group_by_vars))
  })
  
  output$sample_dimplot <- renderPlot({
    dt <- data_temp()
    req(dt, input$sample_select)
    DimPlot(dt, group.by = input$sample_select) + ggtitle(paste("UMAP Plot by", input$sample_select))
  })
  
  output$cluster_dimplot <- renderPlot({
    dt <- data_temp()
    req(dt, input$cluster_select)
    DimPlot(dt, group.by = input$cluster_select) + ggtitle(paste("UMAP Plot by", input$cluster_select))
  })
  
  observeEvent(input$run_subset_btn, {
    obj <- processedObj()
    req(obj)
    sample_var <- input$sample_select
    cluster_var <- input$cluster_select
    sample_vals <- input$sample_subset_vals
    cluster_vals <- input$cluster_subset_vals
    showModal(modalDialog(title = "Please wait",
                          tagList(
                            p("Subsetting..."),
                            p("This step doesn't take long time.")
                          ), footer = NULL, easyClose = FALSE))
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
    selectInput("step3_start_cluster", "Select root Cell cluster:", choices = vals, selected = vals[1])
  })
  
  observeEvent(input$run_sling_btn, {
    req(processedObj())
    obj <- processedObj()
    sample <- input$sample_select
    clusters <- input$cluster_select
    reduction <- input$step3_reduction
    start_cluster <- input$step3_start_cluster
    showModal(modalDialog(title = "Please wait",
                          tagList(
                            p("Running Slingshot..."),
                            p("The waiting time is related to the the number of cell. It takes several minutes to half hour (if the input data is large and no subset is performed).")
                          ), footer = NULL, easyClose = FALSE))
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
  
  output$lineage_select_ui <- renderUI({
    req(step3_run_done())    
    req(step3_result())
    sce <- step3_result()
    
    lin_names <- names(slingCurves(sce))
    if(length(lin_names) == 0) return(NULL)
    
    selectInput(
      inputId = "selected_lineage",
      label = "Select Lineage:",
      choices = lin_names,
      selected = lin_names[1]
    )
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
    req(input$selected_lineage) 
    
    sce <- step3_result()
    reduction <- input$step3_reduction
    xvar <- paste0(reduction, "_1")
    yvar <- paste0(reduction, "_2")
    
    rd_df <- as.data.frame(reducedDims(sce)$Reduction[, 1:2])
    colnames(rd_df) <- c(xvar, yvar)
    rd_df$cluster <- sce$cluster
    
    lineages <- slingCurves(sce)
    curve_data <- data.frame()
    for (i in seq_along(lineages)) {
      curve_i <- as.data.frame(lineages[[i]]$s[,1:2])
      colnames(curve_i) <- c(xvar, yvar)
      curve_i$lineage <- names(lineages)[i]
      curve_data <- rbind(curve_data, curve_i)
    }
    
    curve_data_red <- subset(curve_data, lineage == input$selected_lineage)
    curve_data_black <- subset(curve_data, lineage != input$selected_lineage)
    
    umap_plot <- ggplot(rd_df, aes(x = !!rlang::sym(xvar), y = !!rlang::sym(yvar), color = cluster)) +
      geom_point(size = 0.5) +
      labs(x = xvar, y = yvar, color = "Cluster Labels") +  
      theme_minimal() +
      theme(
        legend.title = element_text(size = 14),     
        legend.text = element_text(size = 12)        
      )+
      guides(color = guide_legend(
        override.aes = list(size = 5) 
      ))
    
    umap_plot + 
      geom_path(data = curve_data_black,
                aes(x = !!rlang::sym(xvar), y = !!rlang::sym(yvar), group=lineage),
                color = "black",
                linewidth = 1) +
      geom_path(data = curve_data_red,
                aes(x = !!rlang::sym(xvar), y = !!rlang::sym(yvar), group=lineage),
                color = "red",
                linewidth = 1.5)
  })
  
  
  output$step3_plot_pseudotime <- renderPlot({
    req(step3_result())
    req(input$selected_lineage)
    
    sce <- step3_result()
    reduction <- input$step3_reduction
    xvar <- paste0(reduction, "_1")
    yvar <- paste0(reduction, "_2")
    
    lineages <- slingCurves(sce)
    sel_lineage_name <- input$selected_lineage
    sel_index <- which(names(lineages) == sel_lineage_name)
    if(length(sel_index) == 0) sel_index <- 1  
    
    pseudotime_col <- paste0("slingPseudotime_", sel_index)
    
    df <- data.frame(
      Dim_1 = reducedDims(sce)$Reduction[,1],
      Dim_2 = reducedDims(sce)$Reduction[,2],
      cluster = colData(sce)$cluster
    )
    
    pt <- colData(sce)[[pseudotime_col]]
    if(is.null(pt)){
      pt <- rep(NA, nrow(df))
    }
    df$pseudotime <- pt
    
    ggplot(df, aes(Dim_1, Dim_2, color = pseudotime)) +
      geom_point(size = 1.5) +
      scale_color_viridis_c(option = "D", na.value = "grey50") +
      labs(title = paste0("Slingshot colored by pseudotime: ", sel_lineage_name), color = "pseudotime") +
      theme_classic()+
      theme(
        legend.title = element_text(size = 14),     
        legend.text = element_text(size = 12),      
        legend.key.size = unit(1, "cm")              
      )
  })
  
  output$step3_next_ui <- renderUI({
    disabled <- !step3_run_done()
    actionButton("next_btn", "NEXT STEP",
                 class = if (disabled) "btn btn-secondary" else "btn btn-primary",
                 disabled = disabled)
  })
  
  observeEvent(input$run_filtering_btn, {
    req(processedObj(), step3_result())
    obj <- processedObj()
    sce <- step3_result()
    non_zero_num <- input$non_zero_num
    lower_quantile <- input$lower_quantile
    upper_quantile <- input$upper_quantile
    if (upper_quantile <= lower_quantile) {
      showModal(modalDialog(title = "Error", "Upper quantile must be greater than lower quantile.", easyClose = TRUE))
      return()
    }
    showModal(modalDialog(title = "Please wait", "Filtering...", footer = NULL, easyClose = FALSE))
    tryCatch({
      res <- step4(obj, sce, non_zero_num, lower_quantile, upper_quantile, input$selected_lineage)
      step4_result(res)
      filter_run_done(TRUE)
      obj_filtered_temp(res$obj)
      
      updateNumericInput(session, "non_zero_num", value = res$non_zero_num)
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
    sce <- step3_result()
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
    Q1_all <- quantile(df$pseudotime, probs = input$lower_quantile/100, na.rm = TRUE)
    Q3_all <- quantile(df$pseudotime, probs = input$upper_quantile/100, na.rm = TRUE)
    
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
      title = "Lamian Data Preprocessing",
      tagList(
        p("Your data is being prepared. Please wait..."),
        p("This step usually takes a few minutes.")
      ),
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
        sce = step3_result(),
        selected_lineage = input$selected_lineage,
        output_file_path = output_file_path,
        design = design_conv,
        maximumknotallowed = maximumknotallowed,
        permuiter = permuiter,
        project_name = project_name,
        email_address = email_address
      )
      
      removeModal()
      
      showModal(
        modalDialog(
          title = "Lamian Job Submitted",
          tagList(
            p("The Lamian analysis is now running in the background."),
            p("Your job ID is:"),
            tags$pre(style = "color: red; user-select: text;", step5_result),
            p("You can copy the following command in the SCC terminal to check the job status:"),
            tags$pre(style = "color: red; user-select: text;", paste0("qstat -j ", step5_result)),
            if (nchar(email_address) > 0) {
              p("The results will be sent to your email address: ", email_address)
            } else {
              NULL
            }
          ),
          footer = modalButton("Close"),
          easyClose = TRUE
        )
      )
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", paste0("Lamian error:\n", e$message), easyClose = TRUE))
    })
  })
  
  output$nav_buttons <- renderUI({
    step <- current_step()
    if (step == 1) {
      btn_label <- if (isTRUE(input$do_harmony)) {
        "INITIATE Normalization with Harmony"
      } else {
        "INITIATE Normalization"
      }
      
      fluidRow(
        column(8,
               div(
                 actionButton("run_btn_ui", btn_label,
                              class = if (dataLoaded()) "btn btn-primary" else "btn btn-secondary",
                              disabled = !dataLoaded()),
                 style = "margin-bottom: 15px;"
               ),
               actionButton("next_btn", "NEXT STEP",
                            class = if (runDoneStep1()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !runDoneStep1())
        )
      )
    } else if (step == 2) {
      fluidRow(
        column(8,
               actionButton("run_subset_btn", "INITIATE SUBSET",
                            class = if (runDoneStep1()) "btn btn-primary" else "btn btn-secondary",
                            disabled = !runDoneStep1(),
                            style = "width: 100%; text-align: center;")
        ),
        column(4,
               actionButton("next_btn", "NEXT STEP",
                            class = "btn btn-primary",
                            disabled = FALSE,
                            style = "width: 100%; text-align: center;")
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
    actionButton("next_btn", "NEXT STEP",
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
      req(step3_result())
      req(input$selected_lineage)
      req(processedObj())
      sce <- step3_result()
      obj <- processedObj()
      pseudotime <- slingPseudotime(sce, na = TRUE)[,input$selected_lineage]
      keep_cells <- names(pseudotime)[!is.na(pseudotime)]
      sce_sub <- sce[, colnames(sce) %in% keep_cells]
      obj_sub <- subset(obj, cells = keep_cells)
      processedObj(obj_sub)
      step3_result(sce_sub)
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
  dataVisualization_server(input, output, session)
}

shinyApp(ui, server)