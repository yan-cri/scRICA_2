## plotPseudotime(): used to perform functional pseudotime analysis                      ##
## Developed by Yan Li, May, 2021
##--------------------------------------------------------------------------------------##
#' plotPseudotime() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param getClusterPseudoRes required, results from fn getClusterPseudo()
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
#' @keywords plotPseudotime
#' @examples plotPseudotime()
#' @export
#'
#' @return
#' 2 plots (correlation plots and 3 combined pseudotime plots) with this specified cluster name in input 'getClusterPseudoRes$clusterName'
#'
#' ## ------------------------------------------------------------------------------------ ##
plotObjPseudotime  <- function(getClusterPseudoRes) {
  ## ---
  pseudoRes        <- getClusterPseudoRes$pseudoRes
  resDir           <- getClusterPseudoRes$resDir
  clusterName      <- getClusterPseudoRes$clusterName
  plotFnamePrefix  <- paste(resDir, clusterName, sep = '/')
  ##' Note: PCA, diffusion map used ranked results as pseudotime,
  ##' whereas slingshot used slingshot() returned results - column 'slingPseudotime_1'.
  ##' row1: 3 methods PC1 vs PC2 pseudotime plot
  ##' row2: 3 methods yCol/colCol on rank(PC1) or slingPseudotime_1
  ##' correlation plot of 3 methods
  ## ---
  R1  <- plotPseudotime(pseudoRes = pseudoRes, xCol = 'PCApc1', xLabel = 'PC1', yCol = 'seurat_clusters', yLabel = 'clusters', colCol = 'seurat_clusters')
  R1  <- R1 + ggtitle(paste(sprintf('PCA: %s', clusterName))) + theme(plot.title = element_text(hjust=0.5))
  R1  <- R1 + theme(legend.position = "none")
  ## -
  R2  <- plotPseudotime(pseudoRes = pseudoRes, xCol = 'dmPc1', xLabel = 'DC1', yCol = 'seurat_clusters', yLabel = 'clusters', colCol = 'seurat_clusters')
  R2  <- R2 + ggtitle(paste(sprintf('Diffusion Map: %s', clusterName))) + theme(plot.title = element_text(hjust=0.5))
  R2  <- R2 + theme(legend.position = "none")
  ## -
  R3  <- plotPseudotime(pseudoRes = pseudoRes, xCol = 'slingPseudotime_1', xLabel = 'slingshot pseudotime', yCol = 'seurat_clusters', yLabel = 'clusters', colCol = 'seurat_clusters')
  R3  <- R3 + ggtitle(paste(sprintf('Slingshot: %s', clusterName))) + theme(plot.title = element_text(hjust=0.5))
  R3  <- R3 + theme(legend.position = "none")
  ## ---
  p1  <- plotPseudotime(pseudoRes = pseudoRes, xCol = 'PCApc1', xLabel = 'PC1', yCol = 'PCApc2', yLabel = 'PC2', colCol = 'seurat_clusters')
  p1  <- p1 + theme(legend.position = "none")
  ## -
  p2  <- plotPseudotime(pseudoRes = pseudoRes, xCol = 'dmPc1', xLabel = 'DC1', yCol = 'PCApc2', yLabel = 'DC2', colCol = 'seurat_clusters')
  ## -
  df_pseudotime           <- as.data.frame(as.data.frame(colData(pseudoRes))[, c("ptPC1", "dmapDptRank", "slingPseudotime_1")]) ##consistent with tutorial exploration with rank result for PCA and DM
  colnames(df_pseudotime) <- c("PCA", "Diffusion", "Slingshot")
  pdf(file = sprintf('%s_pseudotime_correlation.pdf', plotFnamePrefix), width = 6, height = 6.5)
  corrplot.mixed(cor(df_pseudotime, use = "na.or.complete"),
                 order = "hclust", tl.col = "black",
                 # main = "Correlation matrix for pseudotime results",
                 mar = c(0, 0, 2, 3))
  dev.off()
  ## -
  p <- plot_grid(R1, R2, R3, p1, p2, align = 'h', ncol = 3, nrow = 2)
  save_plot(sprintf('%s_pseudotime_plots.pdf', plotFnamePrefix), p, base_height = 8, base_width = 12)
  ## ---
  return(p)
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
