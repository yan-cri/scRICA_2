## compAnno(): compare the Idents annoation results from 2 seurat objects               ##
##--------------------------------------------------------------------------------------##
#' compAnno() Function
#' @details
#' This function is used to perform functional pseudo time analysis via PCA, Diffusion Map, and slingshot
#' @param obj either a Seurat or SingleCellExperiment object.
#' @param slingshotclusterLabels cluster used for slingshot with 3 options: NULL, GMM, Seurat_clusters, by default NULL.
#' @param resSave whether to save the analysis results into Rdata, by default 'FALSE'
#' @param resFnamePrefix if 'resSave = T', the Rdata file name prefix, by default 'TEST'
#'
#' @importFrom slingshot slingshot
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom tradeSeq fitGAM
#' @importFrom tradeSeq associationTest
#' @importFrom SingleCellExperiment reducedDims
#'
#'
#' @keywords compAnno
#' @examples compAnno()
#' @export
#'
#' @return
#' results in cell number and relative percentage of idents cells.
## ------------------------------------------------------------------------------------ ##
compAnno <- function(seuratObj, seuratObj2Comp) {
  ## ------
  cellName1       <- (sapply(strsplit(rownames(seuratObj@meta.data), split = '_'), '[[', 1))
  seuratObj       <- AddMetaData(object = seuratObj, metadata = as.vector(paste(cellName1, seuratObj$expCond, sep = '-')), col.name = 'cellName' )
  cellName2       <- (sapply(strsplit(rownames(seuratObj2Comp@meta.data), split = '_'), '[[', 1))
  seuratObj2Comp <- AddMetaData(object = seuratObj2Comp, metadata = as.vector(paste(cellName2, seuratObj2Comp$expCond, sep = '-')), col.name = 'cellName' )
  vennComp <- gplots::venn(data = list('org' = seuratObjFinal1@meta.data$cellName, 'comp' = seuratObjFinal2@meta.data$cellName))
  ## ------
  seuratObj.clusterMatch           <- seuratObj@meta.data
  seuratObj.clusterMatch$annoOrg   <- as.factor(Idents(seuratObj))
  seuratObj.clusterMatch$anno2Comp <- as.factor(Idents(seuratObj2Comp)[match(seuratObj@meta.data$cellName, seuratObj2Comp@meta.data$cellName)] )
  print(table(seuratObj.clusterMatch$annoOrg))
  print(table(seuratObj.clusterMatch$anno2Comp))
  ## ------
  print(sprintf("%s cells are used to make 2x2 comparison table", dim(seuratObj.clusterMatch)[1] ))
  ## RDS metadata 2 columns with annotation to compare the results with 2x2 table
  clusterLevels <- levels(factor(seuratObj.clusterMatch$anno2Comp))
  table(seuratObj.clusterMatch$annoOrg)
  table(seuratObj.clusterMatch$anno2Comp)
  clusterMatchDf    <- data.frame('annoOrg'= seuratObj.clusterMatch$annoOrg, 'anno2Comp' = seuratObj.clusterMatch$anno2Comp)
  # rownames(clusterMatchDf) <- rownames(seuratObj.clusterMatch@meta.data)
  ## -------
  clusterAnnoRes <- matrix(data = 0, nrow = length(clusterLevels), ncol = length(levels(factor(seuratObj.clusterMatch$annoOrg))) )
  clusterAnnoResPer <- clusterAnnoRes
  colnames(clusterAnnoRes)    <- levels(factor(seuratObj.clusterMatch$annoOrg))
  colnames(clusterAnnoResPer) <- levels(factor(seuratObj.clusterMatch$annoOrg))
  for (i in 1:length(clusterLevels)) {
    # print(sprintf("processing cluster level %s: %s", i, clusterLevels[i]))
    clusterREs <- clusterMatchDf %>% dplyr::filter(anno2Comp == as.character(clusterLevels[i])) %>% as.data.frame()
    # print(sort(table(clusterREs$annoOrg), decreasing = T))
    clusterAnnoSorted     <- sort(table(clusterREs$annoOrg), decreasing = T)
    clusterAnnoSorted2    <- clusterAnnoSorted[sort(names(clusterAnnoSorted))]
    clusterAnnoSorted2Per <- round(clusterAnnoSorted2*100/sum(clusterAnnoSorted2), digits = 1)
    clusterAnnoRes[i,]    <- clusterAnnoSorted2[match(colnames(clusterAnnoRes), names(clusterAnnoSorted2))]
    clusterAnnoResPer[i,] <- clusterAnnoSorted2Per[match(colnames(clusterAnnoRes), names(clusterAnnoSorted2Per))]
    # print('*******************************************************')
  }
  rownames(clusterAnnoRes)   <- clusterLevels
  rownames(clusterAnnoResPer)<- clusterLevels
  clusterAnnoRes[is.na(clusterAnnoRes)] <- 0
  clusterAnnoResPer[is.na(clusterAnnoResPer)] <- 0
  ## ---
  clusterAnnoResSorted    <- clusterAnnoRes[names(sort(apply(clusterAnnoRes, 1, which.max))),]
  clusterAnnoResPerSorted <- clusterAnnoResPer[names(sort(apply(clusterAnnoResPer, 1, which.max))),]
  return(list('count' = clusterAnnoResSorted, 'per' = clusterAnnoResPerSorted))
}

## ------------------------------------------------------------------------------------ ##
