## Developed by Yan Li, June, 2021                                                        ##
##----------------------------------------------------------------------------------------##
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(gridExtra)
## -------------------------------------------------------------------------------------- ##
#' getExpCondClusterMarkers() Function
#' @details
#' This function is used to identify positively expressed cluster markers for all originally identified/annotated cell clusters from the experimental conditions specified in the metadata table.
#'
#' @param resDir full path of integration results analysis are saved, where RDS file is saved inside the 'RDS_Dir'. This path is also returned by getClusterMarkers() execution.
#' @param rdsFname User also can provide the full path of RDS file instead of 'resDir' where RDS file is saved in. If this option is used, please also provide 'resDir' to specify where the analysis results will be saved.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param pAdjValCutoff adjusted p-value cutoff for significant positively expressed cluster markers, by default = 0.05.
#' @param topNo specify the top number of up and down significant positively expressed cluster markers in each identified/annotated clusters, by default = 10.
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
getExpCondClusterMarkers <- function(resDir=NULL, rdsFname=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, expCondCheck='sample', expCondSepName = NULL, pAdjValCutoff = 0.05, topNo = 10) {
  options(java.parameters = "-Xmx32000m")
  ###--------------------------------------------------------------------------------------##
  pAdjValCutoff           <- as.numeric(pAdjValCutoff)
  topNo                   <- as.numeric(topNo)
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) print("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  ##--------------------------------------------------------------------------------------##
  if (is.null(resDir) & !is.null(rdsFname)) {
    rdsFname              <- rdsFname
    resDir                <- getwd()
  } else if (!is.null(resDir) & is.null(rdsFname)) {
    rdsFname              <- sprintf('%s/RDS_Dir/%s.rds', resDir, basename(resDir))
    resDir                <- resDir
  } else {
    stop("Error: please provide either option 'resDir' or 'rdsFname'. ")
  }
  ##--------------------------------------------------------------------------------------##
  if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
  seuratObjFinal          <<- readRDS(file = as.character(rdsFname))
  print('Done for RDS readin')
  # if (clusterReOrder) {
  #   if (length(reorderedClusters)!=length(levels(Idents(seuratObjFinal)))) stop("Please provide corresponding reorderedClusters options")
  # }
  ##--------------------------------------------------------------------------------------##
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir <- paste(resDir, 'clusterMarkerGenes_results_wNewAnnotation', sep = '/')
  } else {
    resDir <- paste(resDir, 'clusterMarkerGenes_results_wOrgClusterAnnotation', sep = '/')
    # if (clusterReOrder) {
    #   resDir <- paste(resDir, 'clusterMarkerGenes_results_wOrgClusterAnnotation_Reorder', sep = '/')
    # } else {
    #   resDir <- paste(resDir, 'clusterMarkerGenes_results_wOrgClusterAnnotation', sep = '/')
    # }
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  # ##--------------------------------------------------------------------------------------##
  # if (expCondCheck == 'sample') {
  #   if (is.null(expCondCheckFname)) {
  #     expCondCheckFname        <- 'expCond_sample'
  #   } else {
  #     expCondCheckFname        <- expCondCheckFname
  #   }
  # } else {
  #   if (is.null(expCondCheckFname)) {
  #     expCondCheckFname        <- as.character(expCondCheck)
  #   } else {
  #     expCondCheckFname        <- expCondCheckFname
  #   }
  # }
  ##--------------------------------------------------------------------------------------##
  print(sprintf('Cluster marker genes identification on cell clusters with respect experimental condition %s will be saved at %s', expCondSepName, resDir))
  resDir                <- paste(sprintf('%s/%s', resDir, expCondSepName ))
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  }
  print(table(Idents(seuratObjFinal)))
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond'
  if (expCondCheck == 'sample') {
    seuratObjFinal                     <- seuratObjFinal
  } else if (expCondCheck == 'idents') {
    seuratObjFinal@meta.data$expCond   <- Seurat::Idents(seuratObjFinal)
  } else {
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(as.character(expCondCheck), colnames(seuratObjFinal@meta.data))]
    }
  }
  ##--------------------------------------------------------------------------------------##
  ## updated expCond factor after above updating
  expCondLevels           <- levels(factor(seuratObjFinal@meta.data$expCond))
  ## ---
  # Seurat::DefaultAssay(seuratObjFinal) <- "integrated"
  expCondPosMarkers       <- list()
  expCondSigPosMarkers    <- list()
  ## -------------------------------------------------------------------------------------
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
    ## identify significant positively expressed cluster markers
    if (assay == 'integrated') {
      allPosMarkers <- FindAllMarkers(seuratObjFinalexpCond, assay = assay, slot = "scale.data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    } else {
      allPosMarkers <- FindAllMarkers(seuratObjFinalexpCond, assay = assay, slot = "data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    }
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
    ## all significant cluster markers heatmap
    cluterAllsigMarkerheatmap <- DoHeatmap(seuratObjFinalexpCond, assay = assay, slot = "scale.data", features = allPosMarkersAdjSig$gene)
    pdf(file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig_heatmap.pdf', resDir, expCondSepName, expCondLevel )), width = 25, height = 20)
    print(cluterAllsigMarkerheatmap)
    dev.off()
    print('END: Step 6 making cluster marker genes heatmap plot')
    print('********************')
    ## ---
  }
  names(expCondPosMarkers)    <- expCondLevels
  names(expCondSigPosMarkers) <- expCondLevels
  ##--------------------------------------------------------------------------------------##
  return(list('expCondPosMarkers' = expCondPosMarkers, 'expCondSigPosMarkers' = expCondSigPosMarkers))
  ##--------------------------------------------------------------------------------------##
}
## -------------------------------------------------------------------------------------- ##
