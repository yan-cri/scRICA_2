## plotPseudotime(): used to perform functional pseudotime analysis                      ##
## Developed by Yan Li, May, 2021
##--------------------------------------------------------------------------------------##
#' plotPseudotime() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param sceObj2plot required, returned results from above fn calcFnPseudo()
#' @param yCol required, which column in the provided object needs to be visualized
#' @param yLabel optional, y axis lable, by default = 'yCol'
#' @param colCol required, which column in the provided object needs to be visualized via different colour
#' @param colLabel optional, legend color lable, by default = 'colCol'
#' @param plotFnamePrefix optional, returned 2 plots (correlation plots and 3 combined pseudotime plots) with this defined 'plotFnamePrefix'
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
#' 2 plots (correlation plots and 3 combined pseudotime plots) with this defined 'plotFnamePrefix'
#'
#' ## ------------------------------------------------------------------------------------ ##
plotPseudotime     <- function(sceObj2plot, yCol, yLabel, colCol, colLabel, plotFnamePrefix) {
  ## ---
  ##' Note: PCA, diffusion map used ranked results as pseudotime,
  ##' whereas slingshot used slingshot() returned results - column 'slingPseudotime_1'.
  if(missing(yLabel)) yLabel <- as.character(yCol)
  if(missing(colLabel)) colLabel <- as.character(colCol)
  if(missing(plotFnamePrefix)) plotFnamePrefix <- 'TEST'
  ## -
  ## rename defined y (yCluster) and col column into 'y' and 'col
  df2plotOrg       <- as.data.frame(colData(sceObj2plot))
  if (yCol == colCol) {
    df2plot        <- df2plotOrg %>% rename_at(match(yCol, colnames(df2plotOrg)), ~'y')
  } else {
    df2plot        <- df2plotOrg %>% rename_at(match(yCol, colnames(df2plotOrg)), ~'y')
    df2plot        <- df2plot %>% rename_at(match(colCol, colnames(df2plot)), ~'col')
  }

  # df2plot$slingshotRank <- rank(df2plot$slingPseudotime_1)
  ## -
  if (yCol == colCol) {
    p1 <- ggplot(df2plot, aes(x = ptPC1, y = y, colour = y))
  } else {
    p1 <- ggplot(df2plot, aes(x = ptPC1, y = y, colour = col))
  }
  p1 <- p1 + geom_quasirandom(groupOnX = FALSE)
  p1 <- p1 + scale_color_tableau() + theme_classic()
  p1 <- p1 + labs(color = as.character(colLabel))
  p1 <- p1 + xlab("PCA PC1 pseudotime") + ylab(as.character(yLabel))
  p1 <- p1 + ggtitle("PCA") + theme(plot.title = element_text(hjust=0.5))
  p1 <- p1 + theme(legend.position = "none")
  print(p1)
  ## -
  if (yCol == colCol) {
    p2 <- ggplot(df2plot, aes(x = dmapDptRank, y = y, colour = y))
  } else {
    p2 <- ggplot(df2plot, aes(x = dmapDptRank, y = y, colour = col))
  }
  p2 <- p2 + geom_quasirandom(groupOnX = FALSE)
  p2 <- p2 + scale_color_tableau() + theme_classic()
  p2 <- p2 + labs(color = as.character(colLabel))
  p2 <- p2 + xlab("Diffusion map pseudotime (dpt)") + ylab(yLabel)
  p2 <- p2 + ggtitle("Diffusion map") + theme(plot.title = element_text(hjust=0.5))
  p2 <- p2 + theme(legend.position = "none")
  print(p2)
  ## -
  if (yCol == colCol) {
    p3 <- ggplot(df2plot, aes(x = slingPseudotime_1, y = y, colour = y))
  } else {
    p3 <- ggplot(df2plot, aes(x = slingPseudotime_1, y = y, colour = col))
  }
  p3 <- p3 + geom_quasirandom(groupOnX = FALSE)
  p3 <- p3 + scale_color_tableau() + theme_classic()
  p3 <- p3 + labs(color = as.character(colLabel))
  p3 <- p3 + xlab("Slingshot pseudotime") + ylab(yLabel)
  p3 <- p3 + ggtitle("Slingshot") + theme(plot.title = element_text(hjust=0.5))
  print(p3)
  ## -
  df_pseudotime <- as.data.frame(df2plot[, c("dmPc1", "dmapDpt1", "slingPseudotime_1")])
  colnames(df_pseudotime) <- c("PC1", "diffusion", "slingshot")
  pdf(file = sprintf('%s_pseudotime_correlation.pdf', plotFnamePrefix), width = 6, height = 6.5)
  corrplot.mixed(cor(df_pseudotime, use = "na.or.complete"),
                 order = "hclust", tl.col = "black",
                 # main = "Correlation matrix for pseudotime results",
                 mar = c(0, 0, 2, 3))
  dev.off()
  ## -
  p <- plot_grid(p1, p2, p3, align = 'hv', ncol = 3, nrow = 1)
  save_plot(sprintf('%s_pseudotime_plots.pdf', plotFnamePrefix), p, base_height = 4, base_width = 12)
  ## ---
  return(p)
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
