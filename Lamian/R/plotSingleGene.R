#' Plot gene(s) by showing both of the original cellcular gene expression and the sample-level, population-level fitting values.
#'
#' This function is used to plot gene(s) by showing both of the original cellcular gene expression and the sample-level, population-level fitting values.
#'
#' @import ggplot2 RColorBrewer splines gridExtra viridis
#' @importFrom grDevices colorRampPalette
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param testobj object returned from lamian_test().
#' @param gene a character vector of gene names. It can be of length 1 or > 1.
#' @param type One of c('Time', 'Variable').
#' @param variable character, the variable (covariate) to color the samples, should be null or one of the column names of design matrix. Default is NULL, meaning each sample is colored differently. Otherwise, samples are colored by the variable (covariate) values.
#' @param variable.text a character vector. The text for the legend of the plot, corresponding to each variable values.
#' @param facet.sample logical. If TRUE (default), facet_wrap the samples.
#' @param plot.point point size
#' @param line.alpha alpha value of the curves
#' @param continuous if TRUE, samples are colored using viridis continuous colors. If FALSE, RColorBrewer "Dark2" discrete palette.
#' @param cellProp logical. If FALSE (default), plot gene expression. If TRUE, it is cell proportion.
#' @param x.lab a string to indicates x-axis label
#' @param y.lab a string to indicates y-axis label
#' @param point.size the size value of points.
#' @param point.alpha the alpha value of points.
#' @param ylim y-axis limits to be passed to ggplot.
#' @param xlim x-axis limits to be passed to ggplot.
#' @param sep a string in the gene names that needs to replaced with blank.
#' @param free.scale logical. If TRUE, the y-axis is on free scale.
#' @param palette a RColorBrewer palette name.
#' @param ncol number of colums for organizing all genes' plot.
#' @param line.size the size of the curves.
#' @param axis.text.blank logical. If TRUE, leave axis text as blank.
#' @examples
#' data(mantestobj)
#' plotGene(testobj = mantestobj, gene = rownames(mantestobj$populationFit[[1]])[1], variable = 'gender')
plotSingleGene <-
  function(testobj,
           gene,
           variable = NULL,
           variable.text = NULL,
           free.scale = TRUE,
           facet.sample = FALSE,
           plot.point = FALSE,
           line.alpha = 1,
           line.size = 1,
           point.alpha = 1,
           point.size = 0.5,
           continuous = TRUE,
           sep = NA,
           palette = 'Dark2',
           ncol = NULL,
           axis.text.blank = FALSE,
           x.lab = 'Pseudotime',
           y.lab = 'Expression') {
    pseudotime <- testobj[['pseudotime']]
    cellanno <- testobj[['cellanno']]
    colnames(cellanno) <- c('Cell', 'Sample')
    expression <- testobj[['expr']][gene, , drop = FALSE]
    
    predict.values <-
      predict_fitting(testobj, gene = gene, test.type = testobj$test.type)

    
  }
