## getClusterPseudo(): used to perform functional pseudotime analysis on seurat object selected clusters##
## Developed by Yan Li, May, 2021
##--------------------------------------------------------------------------------------##
#' getClusterPseudo() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param resDir full path of integration results analysis returned in getClusterMarkers()
#' @param newAnnotation logical value, whether to provide manual annotation
#' @param newAnnotationRscriptName if newAnnotation == T, this script is used to redefine the old clusters
#' @param clusterName specified a single or several cluster names used for pseudotime analysis
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @importFrom Seurat RenameIdents
#'
#' @keywords getClusterPseudo
#' @examples getClusterPseudo()
#' @export
#'
#' @return
#' a SingleCellExperiment where 3 methods functional pseudotime analysis results are saved in colData
## ------------------------------------------------------------------------------------ ##
getClusterPseudo   <- function(resDir, newAnnotation, newAnnotationRscriptName, clusterName) {
  ## ------
  if (missing(newAnnotation)) newAnnotation <- as.logical(F)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Please provide corresponding 'newAnnotationRscriptName', becasue 'newAnnotation' == True.")
  ## ---
  rdsFname                <- paste(resDir, "RDS_Dir/analysis_results_integration_results.rds", sep = '/' )
  ## above is the same as other functions defined as below
  # resDirName            <- strsplit(x = resDir, split = '/')[[1]][length(strsplit(x = resDir, split = '/')[[1]])]
  # rdsFname              <- paste(resDir, sprintf("RDS_Dir/%s.rds", resDirName), sep = '/' )
  ## --
  # print(sprintf('85858 rdsFname is %s', rdsFname))
  if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
  seuratObjFinal          <<- readRDS(file = as.character(rdsFname))
  print('Done for RDS readin')
  ## ------
  ## ------
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir <- paste(resDir, 'results_wNewAnnotation_pseudoTime', sep = '/')
  } else {
    resDir <- paste(resDir, 'results_wOrgClusterAnnotation_pseudoTime', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('Pseudotime analysis results will be saved at %s', resDir))
  ## -------------------------------------------------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  }
  print('=====================================================')
  print(sprintf("START functional pseudo time trajectory analysis for subset '%s' at identified cell clusters shown as below", paste(clusterName, collapse = ', ')  ))
  print(table(Seurat::Idents(seuratObjFinal)))
  if (!all(clusterName %in% names(table(Seurat::Idents(seuratObjFinal))))) stop(sprintf("Provided 'clusterName' is/are not in above clusters, please edit it correspondingly."))
  ## -
  seuratObjFinalSubset                <- subset(seuratObjFinal, idents = clusterName)
  print(sprintf("Subsetted objects are shown as below:"))
  print(table(Seurat::Idents(seuratObjFinalSubset)))
  print('-=-=-=-=-')
  ##' if subclusers included, slingshot  used seurat subclusters
  if (length(levels(factor(seuratObjFinalSubset@meta.data$seurat_clusters))) > 1 ) {
    pseudoRes                         <- calcPCApseudo(obj = seuratObjFinalSubset, slingshotclusterLabels = 'seurat_clusters')
  } else {
    pseudoRes                         <- calcPCApseudo(obj = seuratObjFinalSubset)
  }
  # if (subsetCond == 'idents') {
  # above codes
  # } else if (subsetCond == 'expCond') {
  #   print(sprintf("START functional pseudo time trajectory analysis for subset ('%s') at experimental conditions shown as below", paste(clusterName, collapse = ', ') ))
  #   print(table(seuratObjFinal@meta.data$expCond))
  #   if (!all(clusterName %in% names(table(seuratObjFinal@meta.data$expCond)))) stop(sprintf("Provided 'clusterName' is/are not in above clusters, please edit it correspondingly."))
  #   seuratObjFinalSubset               <- subset(seuratObjFinal, subset = expCond == clusterName )
  #   print(sprintf("Subsetted objects are shown as below:"))
  #   print(table(seuratObjFinalSubset@meta.data$expCond))
  #   print('-=-=-=-=-')
  #   pseudoRes                          <- calcPCApseudo(obj = seuratObjFinalSubset, slingshotclusterLabels = 'ident'  )
  # } else {
  #   stop("Please provide 'subsetCond' in either 'idents' or 'expCond'")
  # }
  ## ---
  return(list(pseudoRes = pseudoRes, resDir = resDir, clusterName = clusterName))
  ## ------
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
