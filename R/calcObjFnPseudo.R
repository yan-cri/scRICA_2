## calcPCApseudo(): used to perform functional pseudotime analysis on seurat object     ##
## Developed by Yan Li, May, 2021
##--------------------------------------------------------------------------------------##
#' calcObjPCApseudo() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param resDir full path of integration results analysis returned in getClusterMarkers()
#' @param newAnnotation logical value, whether to provide manual annotation
#' @param newAnnotationRscriptName if newAnnotation == T, this script is used to redefine the old clusters
#' @param expCondSepName character string, user defined name either to be 'org' or any character string
#' @param expCondName2change if above 'expCondSepName' is defined not as 'org', provide the name to be changed
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @importFrom Seurat RenameIdents
#'
#' @keywords calcObjPCApseudo
#' @examples calcObjPCApseudo()
#' @export
#'
#' @return
#' a SingleCellExperiment where 3 methods functional pseudotime analysis results are saved in colData
## ------------------------------------------------------------------------------------ ##
calcObjPCApseudo   <- function(resDir, newAnnotation, newAnnotationRscriptName, expCondSepName, expCondName2change, subsetCond, subsetClusters) {
  ## ------
  rdsFname                <- paste(resDir, "RDS_Dir/analysis_results_integration_results.rds", sep = '/' )
  ## above is the same as other functions defined as below
  # resDirName              <- strsplit(x = resDir, split = '/')[[1]][length(strsplit(x = resDir, split = '/')[[1]])]
  # rdsFname                <- paste(resDir, sprintf("RDS_Dir/%s.rds", resDirName), sep = '/' )
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
  print(sprintf('cluster summary and new plots will be saved at %s', resDir))
  ## -------------------------------------------------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  }
  ## -------------------------------------------------------------------------------------
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding 'plotResDir' for dotplot to save
  print(sprintf("expCondSep is '%s'", expCondSepName))
  if (expCondSepName == 'org') {
    seuratObjFinal                     <- seuratObjFinal
  } else if (expCondSepName == 'comb') {
    seuratObjFinal@meta.data$expCond   <- Seurat::Idents(seuratObjFinal)
  } else {
    seuratObjFinal@meta.data$expCond   <- gsub(pattern = as.character(expCondName2change), replacement = '', x = seuratObjFinal@meta.data$expCond)
  }
  # print('97979799')
  # print(table(seuratObjFinal@meta.data$expCond))
  # print(sprintf('4444 plotResDir is %s', plotResDir))
  ## ------
  print('=====================================================')
  if (subsetCond == 'idents') {
    print(sprintf("START functional pseudo time trajectory analysis for subset '%s' at identified cell clusters shown as below", paste(subsetClusters, collapse = ', ')  ))
    print(table(Seurat::Idents(seuratObjFinal)))
    if (!all(subsetClusters %in% names(table(Seurat::Idents(seuratObjFinal))))) stop(sprintf("Provided 'subsetClusters' is/are not in above clusters, please edit it correspondingly."))
    ## -
    seuratObjFinalSubset               <- subset(seuratObjFinal, idents = subsetClusters)
    print(sprintf("Subsetted objects are shown as below:"))
    print(table(Seurat::Idents(seuratObjFinalSubset)))
    print('-=-=-=-=-')
    pseudoRes                          <- calcPCApseudo(obj = seuratObjFinalSubset, slingshotclusterLabels = 'expCond'  )
  } else if (subsetCond == 'expCond') {
    print(sprintf("START functional pseudo time trajectory analysis for subset ('%s') at experimental conditions shown as below", paste(subsetClusters, collapse = ', ') ))
    print(table(seuratObjFinal@meta.data$expCond))
    if (!all(subsetClusters %in% names(table(seuratObjFinal@meta.data$expCond)))) stop(sprintf("Provided 'subsetClusters' is/are not in above clusters, please edit it correspondingly."))
    seuratObjFinalSubset               <- subset(seuratObjFinal, subset = expCond == subsetClusters )
    print(sprintf("Subsetted objects are shown as below:"))
    print(table(seuratObjFinalSubset@meta.data$expCond))
    print('-=-=-=-=-')
    pseudoRes                          <- calcPCApseudo(obj = seuratObjFinalSubset, slingshotclusterLabels = 'ident'  )
  } else {
    stop("Please provide 'subsetCond' in either 'idents' or 'expCond'")
  }
  ## ---
  return(list(pseudoRes = pseudoRes, resDir = resDir))
  ## ------
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
