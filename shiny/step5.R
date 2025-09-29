step5 <- function(data, sce, selected_lineage, output_file_path, design, maximumknotallowed, permuiter, task_duration, project_name, email_address = NA, path_to_lamian, auto_run){
  if (!startsWith(output_file_path, "/")) {
    output_file_path <- paste0("/", output_file_path)
  }
  
  if (!dir.exists(output_file_path)) {
    dir.create(output_file_path, recursive = TRUE)
  }
  
  data@assays$SCT <- NULL
  pseudotime <- slingPseudotime(sce, na = TRUE)[,selected_lineage]
  names(pseudotime) <- colnames(sce)
  expr <- as.matrix(GetAssayData(data, layer = "data"))
  cellanno <- data.frame(
    Cell = colnames(expr),
    Sample  = as.character(data@meta.data$sample_id),
    stringsAsFactors = FALSE
  )
  pseudotime <- pseudotime[names(pseudotime) %in% cellanno$Cell]
  
  if ("Sample" %in% colnames(design)) {
    rownames(design) <- design$Sample
    design$Sample <- NULL
  }
  if(!("intercept" %in% colnames(design))) {
    design <- cbind(intercept = 1, design)
  }
  
  saveRDS(expr, file.path(output_file_path,"expr.RDS"))
  saveRDS(pseudotime, file.path(output_file_path, "pseudotime.RDS"))
  saveRDS(design, file.path(output_file_path, "design.RDS"))
  saveRDS(cellanno, file.path(output_file_path, "cellanno.RDS"))
  saveRDS(sce, file.path(output_file_path, "sce.RDS"))
  r_script_path <- file.path(output_file_path, "lamian.R")
  
  lamian_r_code <- sprintf(
    'options(warn=-1) # suppress warning messages
library(devtools)
# Change Here to Lamian Folder
devtools::load_all("%s")

library(circlize)
library(RColorBrewer)
library(dplyr)
library(openxlsx)
#source("renv/activate.R")

send_email <- function(to, subject, body) {
  f <- tempfile(fileext = ".txt")
  mail_content <- c(
    paste0("To: ", to),
    paste0("Subject: ", subject),
    "",
    body
  )
  writeLines(mail_content, f)
  cmd <- sprintf("/usr/sbin/sendmail -t < %%s", shQuote(f))  
  ret <- system(cmd)
  unlink(f)
  if (ret != 0) {
    warning("sendmail command failed with exit code ", ret)
  }
}

output_file_path <- "%s"
email_address <- "%s"
maximumknotallowed <- %d
permuiter <- %d


cellanno <- readRDS(file.path(output_file_path,"cellanno.RDS"))
expr <- readRDS(file.path(output_file_path,"expr.RDS"))
design <- readRDS(file.path(output_file_path,"design.RDS"))
pseudotime <- readRDS(file.path(output_file_path,"pseudotime.RDS"))

start_time <- Sys.time()

res <- tryCatch({
  xde_result <- lamian_test(
    expr = expr,
    cellanno = cellanno,
    pseudotime = pseudotime,
    design = design,
    test.type = "variable",
    maxknotallowed = maximumknotallowed,
    testvar = 2,
    permuiter = permuiter,
    test.method = "permutation",
    ncores = 28,
    verbose.output = TRUE,
    sd.adjust = 1e-5
  )
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  saveRDS(xde_result, file.path(output_file_path,"xde_result.RDS"))
  stat <- xde_result$statistics
  write.xlsx(stat, file.path(output_file_path,"stat.xlsx"), rowNames = TRUE)
  
  if (!is.na(email_address)) {
    body <- c(
      "Hi,",
      "",
      sprintf("The procedure has been completed successfully, the result is under folder %%s. Spending %%s %%s", output_file_path, round(as.numeric(elapsed_time),2), attr(elapsed_time, "units")),
      "",
      "Best,",
      "Bingtian"
    )
    send_email(to = email_address, subject = "Program Completed", body = body)
  }
  
  list(success = TRUE, elapsed = elapsed_time)
}, error = function(e) {
  err_msg <- conditionMessage(e)
  
  if (!is.na(email_address)) {
    body <- c(
      "Hi,",
      "",
      sprintf("There is an error (%%s). You can copy this email and send to btye@bu.edu.", err_msg),
      "",
      sprintf("output_file_path = %%s", output_file_path),
      "",
      "Best,",
      "Bingtian"
    )
    send_email(to = email_address, subject = "Program Completed", body = body)
  }
  line  <- sprintf("Spending %%s %%s", as.numeric(elapsed_time), attr(elapsed_time, "units"))
  writeLines(
  line,
  con = file.path(output_file_path, "time_lamian.txt")
)
  list(success = FALSE, error = err_msg)
})
',
path_to_lamian,
output_file_path,
email_address,
maximumknotallowed,
permuiter
  )
  
  writeLines(lamian_r_code, r_script_path)
  
  # generate qsub
  log_path <- file.path(output_file_path, "output.log")
  qsub_script_path <- file.path(output_file_path, "lamian.qsub")
  
  qsub_content <- paste0(
    "#!/bin/bash\n",
    "#$ -N lamian\n",
    "#$ -cwd\n",
    "#$ -j y\n",
    "#$ -o ", log_path, "\n",
    "#$ -pe omp 28\n",
    "#$ -l h_rt=",task_duration,":00:00\n",
    "#$ -l mem_per_core=18G\n",
    "#$ -V\n",
    "#$ -P ", project_name, "\n",
    "\n",
    "module load R/4.4.3\n",
    "source /projectnb/wax-es/Bingtian/Lamian/miniforge_env/etc/profile.d/conda.sh\n",
    "conda activate /projectnb/wax-es/Bingtian/Lamian/miniforge_env\n",
    "\n",
    "Rscript ", r_script_path, "\n"
  )
  
  writeLines(qsub_content, qsub_script_path)
  if(auto_run){
    qsub_output <- system(paste("qsub", qsub_script_path), intern = TRUE)
    job_id <- sub(".*?([0-9]+).*", "\\1", qsub_output)
    return(job_id)
  }
  
}