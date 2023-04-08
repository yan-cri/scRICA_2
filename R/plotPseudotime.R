## plotPseudotime(): used to generate functional pseudo time plots, called by makeIdentPlots.R ##
## Developed by Yan Li, May, 2021
##--------------------------------------------------------------------------------------##
#' plotPseudotime() Function
#' @details
#' This function is used to perform functional pseudo time analysis via PCA, Diffusion Map, and slingshot
#' @param pseudoRes required, returned results from function calcFnPseudo()$pseudoRes
#' @param xCol required, which column in the provided object needs to be visualized as x-axis
#' @param xLabel optional, x axis label, by default = 'xCol'
#' @param yCol required, which column in the provided object needs to be visualized as y-axis
#' @param yLabel optional, y axis label, by default = 'yCol'
#' @param colCol required, which column in the provided object needs to be visualized via different color
#' @param colLabel optional, legend color label, by default = 'colCol'
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
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom ggthemes scale_color_tableau
#' @importFrom cowplot save_plot
#' @importFrom cowplot plot_grid
#' @importFrom dplyr rename_at
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#'
#' @keywords plotPseudotime
#' @examples plotPseudotime()
#' @export
#'
#' @return
#' a plot with corresponding specified x, y-axis used for plot
#'
#' ## ------------------------------------------------------------------------------------ ##
plotPseudotime     <- function(pseudoRes, xCol, xLabel=NULL, yCol, yLabel=NULL, colCol, colLabel=NULL) {
  ## ---
  # library(ggplot2)
  # library(cowplot)
  # library(ggthemes)
  # library(ggbeeswarm)
  # library(dplyr)
  # library(SingleCellExperiment)
  ## -
  if(is.null(xLabel)) xLabel <- as.character(xCol)
  if(is.null(yLabel)) yLabel <- as.character(yCol)
  if(is.null(colLabel)) colLabel <- as.character(colCol)
  ## -
  ## rename defined x, y (yCluster) and col column into 'x', 'y' and 'col
  df2plotOrg       <- as.data.frame(colData(pseudoRes$sceObj))
  # df2plotOrg       <- as.data.frame(res2@colData@listData[1:37]  )
  if (yCol == colCol) {
    df2plot        <- df2plotOrg %>% rename_at(match(xCol, colnames(df2plotOrg)), ~'x')
    df2plot        <- df2plot %>% rename_at(match(yCol, colnames(df2plotOrg)), ~'y')
    p                <- ggplot(df2plot, aes(x = x, y = y, colour = y))
  } else {
    df2plot        <- df2plotOrg %>% rename_at(match(xCol, colnames(df2plotOrg)), ~'x')
    df2plot        <- df2plot %>% rename_at(match(yCol, colnames(df2plotOrg)), ~'y')
    df2plot        <- df2plot %>% rename_at(match(colCol, colnames(df2plot)), ~'col') %>% mutate(col = as.character(col))
    p                <- ggplot(df2plot, aes(x = x, y = y, colour = col))
  }
  p                <- p + geom_quasirandom(groupOnX = FALSE)
  p                <- p + scale_color_tableau() + theme_classic()
  p                <- p + labs(color = as.character(colLabel))
  p                <- p + xlab(as.character(xLabel)) + ylab(as.character(yLabel))
  # p                <- p + ggtitle("PCA") + theme(plot.title = element_text(hjust=0.5))
  # p                <- p + theme(legend.position = "none")
  print(p)
  ## ---
  return(p)
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
