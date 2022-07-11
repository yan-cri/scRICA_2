##--------------------------------------------------------------------------------------##
#' getChiResults() Function
#' @details
#' This function is used to perform chi-square test on the cell number decomposition table to generate the corresponding contribution plot.
#'
#' @param perTable cell number decomposion table, where row represent *** and column represent.
#' @param expCondCheck specify which columns in a regular expression pattern to be test and plot.
#' @param chisqMonteCarlo a logical indicating whether to compute p-values by Monte Carlo simulation, by default = Fasle.
#' @param plotLegendPer plot specification on the legend bar, specify the percentage of legend bar to occupy the entire plot, by default = 0.4.
#' @param plotCex plot specification on the cex size of the plot, by default = 2.
#' @param legendCex plot specification on the cex size of the legend, by default = 2.
#' @param legendDistance plot specification on the distance between legend text and legend bar, by default = 0.3.
#' @param axisLableDistance plot specification on the distance between axis label text and plot, by default = 1.
#' @param topTextDegree plot specification on the rotation degree of top label text, by default = 90 (vertical align).
#' @param plotFnamePrefix plot specification on the plot file name prefix, by default a plot named 'test.pdf' will be generated.
#' @param plotVertical plot specification on the plot display whether to rotate column and row displays on the plot.
#' @param plotWidth plot specification on the plot width size, by default = 7.
#' @param plotHeight plot specification on the plot height size, by default = 7.
#' @param resReturn specify whether chi-square test results to be returned in a list format, by default = False.
#' @param debug by default False, for development purpose to check the analysis progress.
#'
#' @importFrom dplyr %>%
#' @importFrom corrplot corrplot
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom utils write.table
#' @importFrom ggplot2 element_text
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#'
#' @keywords getChiResults
#' @examples getChiResults(perTable)
#' @export
#' @return
#' if specified to return, a list of 3 items can be returned, they are 'chiRes', 'residual' and 'contrib'.
## ------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------ ##
getChiResults <- function(perTable, expCondCheck = NULL, chisqMonteCarlo = F, plotLegendPer = 0.4, plotCex = 2, legendCex = 2, legendDistance = 0.3, axisLableDistance = 1, topTextDegree = 90, plotFnamePrefix = 'test', plotVertical = F, plotWidth = 7, plotHeight = 7, resReturn = F, debug = F) {
  # library(dplyr)
  # library(corrplot)
  # library(ggplot2)
  if (is.null(expCondCheck)) {
    res4test <- perTable
  } else {
    res4test <- perTable[, grep(pattern = sprintf('%s|clusters', expCondCheck), colnames(perTable))]
  }
  rownames(res4test) <- res4test$clusters
  res4test2 <- res4test %>% dplyr::select(-clusters)
  rownames(res4test2) <- res4test$clusters
  colnames(res4test2) <- gsub(pattern = 'cellNo_|_Per', '', colnames(res4test2))
  if (plotVertical) {
    res4test2 <- t(res4test2)
  }
  if (debug) {
    print(head(res4test2))
  }
  ## ---
  if (chisqMonteCarlo) {
    chires    <- chisq.test(x = res4test2, simulate.p.value = T, rescale.p = T)
  } else {
    chires    <- chisq.test(x = res4test2)
  }
  print(chires)
  ## ---
  print(sprintf("Chi-square test p-value is %s with method: %s ", chires$p.value, chires$method))
  chi.obs <- chires$observed
  chi.exp <- chires$expected
  ## calculate number of dependence between the row and column with pearson residuals (r)
  chi.residual <- round(chires$residuals, 3)
  pdf(file = sprintf('%s_associationPlot.pdf', plotFnamePrefix), width = plotWidth, height = plotHeight)
  corrplot(corr = chi.residual, is.cor = F, diag = F, tl.col = 1, tl.cex = plotCex, cl.cex = legendCex, cl.align.text = 'l', cl.ratio = plotLegendPer, cl.offset = legendDistance, tl.offset = axisLableDistance, tl.srt = topTextDegree)
  dev.off()
  ## calculate contribution of a give cell to total chi-square score
  contrib  <- round(100*chi.residual^2/chires$statistic, 3)
  pdf(file = sprintf('%s_contributionPlot.pdf', plotFnamePrefix), width = plotWidth, height = plotHeight)
  corrplot(corr = contrib, is.cor = F, diag = F, tl.col = 1, tl.cex = plotCex, cl.cex = legendCex, cl.align.text = 'l', cl.ratio = plotLegendPer, cl.offset = legendDistance, tl.offset = axisLableDistance, tl.srt = topTextDegree)
  dev.off()
  ## -
  return(list('chiRes'=chires, 'residual'=chi.residual, 'contrib' = contrib))
}
## ------------------------------------------------------------------------------------ ##
