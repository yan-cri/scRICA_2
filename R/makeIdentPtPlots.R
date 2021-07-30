## makeIdentPlots(): used to generate pseudotime times PCA plots of 3 analysis methods, calling plotPseudotime() from plotPseudotime.R  ##
## Developed by Yan Li, June, 2021
##--------------------------------------------------------------------------------------##
#' makeIdentPlots() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param pseudoRes required, results from fn getClusterPseudo()
#' @param plotname pseudo time GMM clustering plot name prefix, by default 'TEST'
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme
#' @importFrom stats cor
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom ggthemes scale_color_tableau
#' @importFrom cowplot save_plot
#' @importFrom cowplot plot_grid
#' @importFrom dplyr rename_at
#' @importFrom dplyr %>%
#' @importFrom corrplot corrplot.mixed
#'
#' @keywords makeIdentPlots
#' @examples makeIdentPlots()
#' @export
#'
#' @return
#' 2 plots (correlation plots and 3 combined pseudo time plots) with this specified cluster name in input 'pseudoRes$clusterName'
#'
#' ## ------------------------------------------------------------------------------------ ##
makeIdentPlots <- function(pseudoRes, plotname = 'TEST') {
  ## ---
  if (missing(pseudoRes)) stop("Please provide pseudotime analysis results in 'pseudoRes' for makeIdentPlots() to make pseudotime PC/DC plots")
  ## ---
  ## Add a warning notice, to be done
  plotTheme       <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                              axis.title = element_text(size = 20),
                              axis.text = element_text(size = 25),
                              legend.position="right",
                              legend.text = element_text(size = 20) )
  ## ---
  R1  <- plotPseudotime(pseudoRes = pseudoRes, xCol = 'slingPseudotime_1', xLabel = 'slingshot pseudotime', yCol = 'PCApc1', yLabel = 'PC1', colCol = 'seurat_clusters')
  # R1  <- R1 + ggtitle(paste(sprintf('PCA: %s', clusterName))) + theme(plot.title = element_text(hjust=0.5))
  R1  <- R1 + plotTheme + NoLegend()
  ## -
  R2  <- plotPseudotime(pseudoRes = pseudoRes, xCol = 'slingPseudotime_1', xLabel = 'slingshot pseudotime', yCol = 'dmPc1', yLabel = 'DC1', colCol = 'seurat_clusters')
  R2  <- R2 + plotTheme
  ## -
  # df_pseudotime           <- as.data.frame(as.data.frame(colData(pseudoRes$sceObj))[, c("ptPC1", "dmapDptRank", "slingPseudotime_1")]) ##consistent with tutorial exploration with rank result for PCA and DM
  # colnames(df_pseudotime) <- c("PCA", "Diffusion", "Slingshot")
  # pdf(file = sprintf('%s_pseudotime_correlation.pdf', plotname), width = 6, height = 6.5)
  # corrplot.mixed(cor(df_pseudotime, use = "na.or.complete"),
  #                order = "hclust", tl.col = "black",
  #                # main = "Correlation matrix for pseudotime results",
  #                mar = c(0, 0, 2, 3))
  # dev.off()
  ## -
  p <- plot_grid(R1, R2, align = 'h', ncol = 2, nrow = 1, rel_widths = c(1,1.2))
  save_plot(sprintf('ptSeuratIndent_%s_plots.pdf', plotFnamePrefix), p, base_height = 5, base_width = 10)
  ## ---
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
