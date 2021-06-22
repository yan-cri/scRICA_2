## seuratIntegrate(): integrate seurat objects with identified anchor via seurat method   ##
## Developed by Yan Li, June, 2021                                                        ##
##----------------------------------------------------------------------------------------##
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(gridExtra)
## -------------------------------------------------------------------------------------- ##
#' getExpCondClusterMarkers() Function
#' @details
#' This function is used to identify positively expressed cluster markers in specified expreimental conditions
#'
#' @param resDir full path of integration results analysis returned in getClusterMarkers()
#' @param newAnnotation logical value, whether to provide manual annotation
#' @param newAnnotationRscriptName if newAnnotation == T, this script is used to redefine the old clusters
#' @param expCondSepName character string, user defined name either to be 'org' or any character string
#' @param expCondName2change if above 'expCondSepName' is defined not as 'org', provide the name to be changed
#' @param pAdjValCutoff adjusted p-value cutoff of significant positively expressed cluster markers
#' @param topNo output this number of top up and down significant positively expressed cluster markers in each identified/annotated clusters
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @importFrom Seurat RenameIdents
#' @importFrom Seurat FindAllMarkers
#' @importFrom utils write.table
#' @importFrom Seurat DoHeatmap
#' @importFrom Seurat NoLegend
#' @importFrom dplyr %>%
#' @importFrom dplyr top_n
#' @importFrom dplyr group_by
#' @importFrom dplyr distinct
#' @importFrom dplyr summarise
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @keywords getExpCondClusterMarkers
#' @examples getExpCondClusterMarkers()
#' @export
#' @return
#' a list item including 3 elements:
## -------------------------------------------------------------------------------------- ##
getExpCondClusterMarkers <- function(resDir, newAnnotation, newAnnotationRscriptName, clusterReOrder, reorderedClusters, expCondSepName, expCondName2change, pAdjValCutoff, topNo) {
  options(java.parameters = "-Xmx32000m")
  if (missing(topNo)) topNo = 10
  if (missing(clusterReOrder)) clusterReOrder = as.logical(F)
  if (missing(reorderedClusters)) reorderedClusters <- ''
  ## --------------------------------------------
  rdsFname                <- paste(resDir, "RDS_Dir/analysis_results_integration_results.rds", sep = '/' )
  # print(sprintf('85858 rdsFname is %s', rdsFname))
  if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
  seuratObjFinal          <<- readRDS(file = as.character(rdsFname))
  print('Original identified clusters idents as below:')
  print(table(Idents(seuratObjFinal)))
  print('Done for RDS readin')
  if (clusterReOrder) {
    if (length(reorderedClusters)!=length(levels(Idents(seuratObjFinal)))) stop("Please provide corresponding reorderedClusters options")
  }
  ## ------
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir <- paste(resDir, 'clusterMarkerGenes_results_wNewAnnotation', sep = '/')
  } else {
    if (clusterReOrder) {
      resDir <- paste(resDir, 'clusterMarkerGenes_results_wOrgClusterAnnotation_Reorder', sep = '/')
    } else {
      resDir <- paste(resDir, 'clusterMarkerGenes_results_wOrgClusterAnnotation', sep = '/')
    }
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('Cluster marker genes identification on cell clusters with respect experimental condition %s will be saved at %s', expCondSepName, resDir))
  resDir                <- paste(sprintf('%s/expCond_%s', resDir, expCondSepName ))
  if (!dir.exists(resDir)) dir.create(resDir)
  ## --------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  } else {
    if (clusterReOrder) {
      # seuratObjFinal2@active.ident <- factor(seuratObjFinal2@active.ident, levels = reorderedClusters)
      Seurat::Idents(seuratObjFinal) <- factor(Seurat::Idents(seuratObjFinal), levels = as.character(reorderedClusters))
      print('Identified cell clusters are re-ordered as below:')
    } else {
      print('Identified cell clusters are sorted as orginal below:')
    }
    print(table(Idents(seuratObjFinal)))
  }
  # ## --------------------------------------------
  # if (expCondSepName == 'A') {
  #   expCondName2change = expCondName2change
  #   # expCondName2change = '3041|3043|3061|3203|3295|3296|3369|3391' ## replace with levels of condition B (patient No.) for FT
  #   # expCondName2change = '3041|3061|3203|3296|3391' ## replace with levels of condition B (patient No.) for Ovary
  # } else if (expCondSepName == 'B') {
  #   expCondName2change = expCondName2change
  #   # expCondName2change = 'A|F|I' ## replace with levels of condition A (tissue type) for FT
  #   # expCondName2change = 'O' ## replace with levels of condition A (tissue type) for Ovary, only 1 level no need to use, the same as 'org'
  # } else {
  #   expCondName2change = expCondName2change
  # }
  ## --------------------------------------------
  if (expCondSepName == 'org') {
    seuratObjFinal        <- seuratObjFinal
  } else {
    seuratObjFinal@meta.data$expCond <- gsub(pattern = as.character(expCondName2change), replacement = '', x = seuratObjFinal@meta.data$expCond)
  }
  ## updated expCond factor after above updating
  expCondLevels           <- levels(factor(seuratObjFinal@meta.data$expCond))
  ## --------------------------------------------
  Seurat::DefaultAssay(seuratObjFinal) <- "integrated"
  expCondPosMarkers       <- list()
  expCondSigPosMarkers    <- list()
  for (l in 1:length(expCondLevels)) {
    # Seurat::DefaultAssay(seuratObjFinal)   <- "RNA"
    ## ---
    systime1              <- Sys.time()
    ## update 'resDir' to ceate dir under 'result_wNewAnnotation' or 'results_wOrgClusterAnnotation'
    expCondLevel          <- expCondLevels[l]
    ## subsetting different experimental conditions
    seuratObjFinalexpCond <- subset(seuratObjFinal, expCond == expCondLevel)
    ## identify cluster positively expressed markers
    print(sprintf('Start %s: finding positive regulated cluster marker genes for experimental condition: %s', l, expCondLevel))
    allPosMarkers         <- FindAllMarkers(seuratObjFinalexpCond, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.table(x = allPosMarkers, file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_full.txt', resDir, expCondSepName, expCondLevel)), quote = F, sep = '\t', row.names = T, col.names = NA)
    print(sprintf('Maximum p_value is %s, Maximum adjusted p_value is %s', round(max(allPosMarkers$p_val), digits = 4), round(max(allPosMarkers$p_val_adj), digits = 4)))
    ## identify significant positively expressed cluster markers
    allPosMarkersAdjSig   <- allPosMarkers %>% dplyr::filter(p_val_adj <= pAdjValCutoff) %>% dplyr::mutate(perDiff = pct.1-pct.2)
    if (dim(allPosMarkersAdjSig)[1] > 0) write.table(x = allPosMarkersAdjSig, file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig.txt', resDir, expCondSepName, expCondLevel)), quote = F, sep = '\t', row.names = T, col.names = NA)
    print(sprintf('A total of %s positively expressed genes identified for experimental condition %s at %s, among them %s are significant up expressed at adjusted p-value significant level of %s', dim(allPosMarkers)[1], expCondSepName, expCondLevel, dim(allPosMarkersAdjSig)[1], pAdjValCutoff ))
    ## summarize the no. of significant positively expressed cluster markers
    allPosMarkersAdjSigNo <- allPosMarkersAdjSig %>% dplyr::group_by(cluster) %>% dplyr::distinct() %>% dplyr::summarise('geneNo' = n()) %>% as.data.frame()
    write.table(x = allPosMarkersAdjSigNo, file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig_NoSummary.txt', resDir, expCondSepName, expCondLevel)), quote = F, sep = '\t', row.names = F, col.names = T)
    ## -
    expCondPosMarkers[[l]]     <- allPosMarkers
    expCondSigPosMarkers[[l]]  <- allPosMarkersAdjSig
    systime2               <- Sys.time()
    print(sprintf('END %s: finding positive regulated cluster marker genes for experimental condition in %s with computation time: %s %s.', l, expCondLevel, round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
    print('Start: Step 6 making cluster marker genes heatmap plot')
    ## top N markers identification
    topMarkers             <- allPosMarkers %>% group_by(cluster) %>% top_n(n = topNo, wt = avg_log2FC) %>% as.data.frame()
    write.table(x = topMarkers, file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_top%s.txt', resDir, expCondSepName, expCondLevel, topNo )), quote = F, sep = '\t', row.names = T, col.names = NA)
    # Seurat::DefaultAssay(seuratObjFinalexpCond) <- "integrated"
    ## top N markers heatmap
    cluterTopMarkerheatmap <- DoHeatmap(seuratObjFinalexpCond, features = topMarkers$gene)
    pdf(file = file.path(sprintf('%s/expCond_%s_%s_allCluster_top%sPosMarkers_heatmap.pdf', resDir, expCondSepName, expCondLevel, topNo )), width = 25, height = 20)
    print(cluterTopMarkerheatmap)
    dev.off()
    print(sprintf('complete top %s cluster marker genes heatmap plot', topNo))
    ## all significant cluster markers heatmap
    cluterAllsigMarkerheatmap <- DoHeatmap(seuratObjFinalexpCond, features = allPosMarkersAdjSig$gene)
    pdf(file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig_heatmap.pdf', resDir, expCondSepName, expCondLevel )), width = 25, height = 20)
    print(cluterAllsigMarkerheatmap)
    dev.off()
    print('END: Step 6 making cluster marker genes heatmap plot')
    print('********************')
    ## ---
  }
  names(expCondPosMarkers)    <- expCondLevels
  names(expCondSigPosMarkers) <- expCondLevels
  ## --------------------------------------------
  return(list('expCondPosMarkers' = expCondPosMarkers, 'expCondSigPosMarkers' = expCondSigPosMarkers))
  ## -------------------------------------------------------------------------------------
}
## -------------------------------------------------------------------------------------- ##
