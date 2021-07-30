## plotPseudotimeLineages(): make slingshot lineage on PCA dimensions                   ##
## Developed by Yan Li, Junly, 2021
##--------------------------------------------------------------------------------------##
#' plotPseudotime() Function
#' @details
#' This function is used to perform functional pseudo time analysis via PCA, Diffusion Map, and slingshot
#' @param pseudoRes required, returned results from function calcFnPseudo()$pseudoRes
#' @param plotname GMM slingshot lineage plot file name prefix, by default 'TEST_ptGMM_lineage.pdf'
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr rename
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom slingshot SlingshotDataSet
#' @importFrom graphics par
#'
#' @keywords plotPseudotimeLineages
#' @examples plotPseudotimeLineages()
#' @export
#'
#' @return
#' generate a PCA plot with slingshot calculated lineage in defined plotname.pdf
#'
#' ## ------------------------------------------------------------------------------------ ##
plotPseudotimeLineages    <- function(pseudoRes, plotname = 'TEST') {
  ## ---
  if (missing(pseudoRes)) stop("Please provide pseudotime analysis results in 'pseudoRes' for plotPseudotimeLineages() to make slingshot pseudotime lineage plots")
  ## ---
  df2plotOrg       <- as.data.frame(colData(pseudoRes$sceObj)) %>% dplyr::select(c('PCApc1', 'PCApc2'))
  df2plot          <- df2plotOrg  %>% dplyr::rename(PC1 = PCApc1) %>% dplyr::rename(PC2 = PCApc2) %>% as.data.frame()
  ## color options
  palettes         <- ggthemes::ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
  selectedcol      <- palettes$`Tableau 10`%>% pull(value)
  ## -
  pdf(file = sprintf('%s_ptGMM_lineage.pdf', plotname), width = 5.2, height = 5.2)
  par(mar = c(4, 4.5, 1, 1) + 0.1, cex.axis = 2, cex.lab = 2 )
  plot(df2plot, col = selectedcol[pseudoRes$sceObj$GMM], pch=16, asp = 1)
  lines(slingshot::SlingshotDataSet(pseudoRes$slingshotRes), lwd=2, type = 'lineages', col = 'black')
  dev.off()
  ## -

}
#' ## ------------------------------------------------------------------------------------ ##
#'
