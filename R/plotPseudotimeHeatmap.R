## plotPseudotimeHeatmap(): used to generate heatmap of GAM ranked genes on slingshot pseudo time cells ##
## Developed by Yan Li, June, 2021
##--------------------------------------------------------------------------------------##
#' plotPseudotimeHeatmap() Function
#' @details
#' This function is used to perform functional pseudo time analysis via PCA, Diffusion Map, and slingshot
#' @param pseudoRes required, returned results from function calcFnPseudo()$pseudoRes$sceObj
#' @param plotname GAM clustering heatmap plot file name prefix, by default 'TEST_ptGMM_GAMheatmap.pdf'
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom gplots heatmap.2
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom dplyr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics par
#'
#' @keywords plotPseudotimeHeatmap
#' @examples plotPseudotimeHeatmap()
#' @export
#'
#' @return
#' a plot with corresponding specified x, y-axis used for plot
#'
#' ## ------------------------------------------------------------------------------------ ##
plotPseudotimeHeatmap   <- function(pseudoRes, plotname = 'TEST') {
  ## ---
  # library(RColorBrewer)
  # library(dplyr)
  # library(SingleCellExperiment)
  ## ---
  if (missing(pseudoRes)) stop("Please provide pseudotime analysis results in 'pseudoRes' for plotPseudotimeHeatmap() to make heatmap plots")
  ## ---
  slingshotPtOrder  <- order(pseudoRes$sceObj$slingPseudotime_1, na.last = NA)
  heatdata          <- as.matrix(assays(pseudoRes$sceObj)$logcounts[pseudoRes$rankGene[1:100], slingshotPtOrder])
  # heatdata2         <- as.matrix(assays(pseudoRes$sceObj)$counts[pseudoRes$rankGene[1:100], slingshotPtOrder])
  heatclus          <- pseudoRes$sceObj$GMM[slingshotPtOrder]
  # heatmap(log1p(heatdata2), Colv = NA, ColSideColors = brewer.pal(length(levels(factor(heatclus))),"Set1")[heatclus])

  ## color options
  palettes         <- ggthemes::ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
  selectedcol      <- palettes$`Tableau 10`%>% pull(value)

  pdf(file = sprintf('%s_ptGMM_GAMheatmap.pdf', plotname), width = 8, height = 6)
  # heatmap(heatdata, Rowv = NA, Colv = NA, scale = 'row',
  #               ColSideColors = brewer.pal(length(levels(factor(heatclus))),"Set1")[heatclus])
  # heatmap(heatdata, Colv = NA, scale = 'row',
  #         ColSideColors = brewer.pal(length(levels(factor(heatclus))),"Set1")[heatclus])
  par(lwd=1.5)
  hmRes <- heatmap.2(as.matrix(heatdata), Rowv=T, Colv=F, distfun=dist,
                     scale="row", hclustfun =hclust, dendrogram="row",
                     cexCol=0.1, labCol = NULL,
                     key.title= NA, key.xlab = 'scaled log counts',
                     key=TRUE, symkey=FALSE, keysize = 1,
                     key.par = list(cex.lab=1, cex.axis = 1),
                     density.info="none", trace="none", cexRow=1,
                     ColSideColors= selectedcol[1:length(levels(factor(heatclus)))][heatclus])
  dev.off()
  write.table(x = rownames(as.matrix(heatdata))[hmRes$rowInd], file =  sprintf('%s_ptGMM_GAMheatmapGenenames.txt', plotname), quote = F, sep = '\t', row.names = F, col.names = F)
  ## ---
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
