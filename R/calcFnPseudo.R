## calcPCApseudo(): used to perform functional pseudotime analysis                      ##
## Developed by Yan Li, May, 2021
##--------------------------------------------------------------------------------------##
#' calcPCApseudo() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param obj required, either a seurat or SingleCellExperiment object
#' @param slingshotclusterLabels default no cluster used for slingshot, 3 options: NULL, GMM, seurat_clusters
#' @param resSave optional, default is F
#' @param resFnamePrefix optional, default is 'TEST', when resSave is on
#'
#' @importFrom slingshot slingshot
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom tradeSeq fitGAM
#' @importFrom tradeSeq associationTest
#' @importFrom SingleCellExperiment reducedDims
#'
#'
#' @keywords calcPCApseudo
#' @examples calcPCApseudo()
#' @export
#'
#' @return
#' a SingleCellExperiment where 3 methods functional pseudotime analysis results are saved in colData
## ------------------------------------------------------------------------------------ ##
calcPCApseudo <- function(obj, slingshotclusterLabels = NULL, resSave = 'F', resFnamePrefix = 'TEST') {
  ## ------
  ## This analysis is based on RNA assay regardless whether integrative assay has been implemented
  ## 0. decide the obj input and change into 'SingleCellExperiment' object
  if (class(obj) == "Seurat") {
    Seurat::DefaultAssay(obj) <- "RNA"
    print(sprintf("currently is operating on %s assay", Seurat::DefaultAssay(obj) ))
    seuratObj     <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seuratObj     <- FindVariableFeatures(seuratObj, selection.method = 'vst', nfeatures = 2000)
    seuratObj     <- ScaleData(object = seuratObj)
    seuratVarFea  <- subset(seuratObj, features = VariableFeatures(object = seuratObj))
    sceObj        <- Seurat::as.SingleCellExperiment( seuratVarFea )
    # head(assays(sceObj)$logcounts[,1:3])
  } else if (class(obj) == "SingleCellExperiment") {
    sceObj         <- obj
  } else {
    stop("please provide either Seurat or SingleCellExperiment object")
  }
  # print(sceObj)
  ## ---
  ## 1. perform BiocSingular PCA analysis
  print('-=-=-=-=-=-=-=-=-=-=')
  print(sprintf('Start step1: BiocSingular PCA pseudotime calculation at %s',  Sys.time() ))
  systime1         <- Sys.time()
  # print('8494934934893843984')
  sceObj           <- scater::runPCA(sceObj, ncomponents = 50)
  # print('8494934934893843984')
  pca              <- SingleCellExperiment::reducedDim(sceObj, "PCA")
  # dim(pca)
  # head(pca[,1:3])
  sceObj$PCApc1    <- pca[, 1]
  sceObj$PCApc2    <- pca[, 2]
  sceObj$ptPC1     <- rank(sceObj$PCApc1)  # rank cells by their PC1 score
  # print(head(as.data.frame(colData(sceObj))))
  systime2         <- Sys.time()
  print(sprintf('End step1: BiocSingular PCA pseudotime calculation at %s',  systime2 ))
  print(sprintf("Step1 PCA analysis used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  ## ---
  ## 2. perform Diffusion map pseudotime
  print('-=-=-=-=-=-=-=-=-=-=')
  print(sprintf('Start step2: Diffusion map pseudotime calculation at %s',  Sys.time() ))
  systime1         <- Sys.time()
  dmPca            <- destiny::DiffusionMap(pca)
  dpt              <- destiny::DPT(dmPca)
  # dpt1 <- DPT(dmPca, tips = 268)
  ## -
  sceObj$dmPc1     <- destiny::eigenvectors(dmPca)[, 1]
  sceObj$dmPc2     <- destiny::eigenvectors(dmPca)[, 2]
  sceObj$dmapDpt1  <- dpt$dpt
  sceObj$dmapDptRank <- rank(dpt$dpt)
  systime2         <- Sys.time()
  print(sprintf('End step2: Diffusion map pseudotime calculation at %s',  systime2 ))
  print(sprintf("Step2 Diffsion Map analysis used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  ## ---
  ## 3. perform Slingshot map pseudotime
  print('-=-=-=-=-=-=-=-=-=-=')
  print(sprintf('Start step3: slingshot pseudotime calculation at %s', Sys.time() ))
  systime1         <- Sys.time()
  if (is.null(slingshotclusterLabels)) {
    print(sprintf('No cluster labels will be used in slingshot calculation' ))
    slingshotRes   <- slingshot::slingshot(sceObj, reducedDim = 'PCA')
  } else if (slingshotclusterLabels == 'GMM') {
    print(sprintf('%s cluster label will be used in slingshot calculation', slingshotclusterLabels ))
    ## Note: expCond level cannot be used with 'seurat_clusters'
    # sce2seurat     <- as.Seurat(sceObj, counts = "counts", data = "logcounts")
    # library(mclust, quietly = TRUE)
    mc               <- mclust::Mclust(reducedDims(sceObj)$PCA[,1:2])$classification
    sceObj$GMM       <- mc
    slingshotRes     <- sceObj
    slingshotRes     <- slingshot::slingshot(sceObj, clusterLabels = 'GMM', reducedDim = 'PCA')
  } else if (slingshotclusterLabels == 'seurat_clusters') {
    print(sprintf('%s cluster label will be used in slingshot calculation', slingshotclusterLabels ))
    slingshotRes     <- slingshot::slingshot(sceObj, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')
  }
  ## ---
  sceObj$slingPseudotime_1 <-  slingshotRes$slingPseudotime_1
  systime2           <- Sys.time()
  print(sprintf('End step3: slingshot pseudotime calculation at %s', systime2 ))
  print(sprintf("Step3 slingshot pseudotime analysis used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  ## ---
  ## 4. identify temporally dynamic genes based on slingshot pseudotime
  print('-=-=-=-=-=-=-=-=-=-=')
  print(sprintf('Start step4: slingshot pseudotime GAM gene calculation at %s', Sys.time() ))
  systime1           <- Sys.time()
  # library(tradeSeq)
  ## fit negative binomial GAM
  sceObjfitGam       <- tradeSeq::fitGAM(slingshotRes)
  # test for dynamic expression
  sceObjfitGamATres  <- tradeSeq::associationTest(sceObjfitGam)
  gamSortedGenes     <- rownames(sceObjfitGamATres[order(sceObjfitGamATres$pvalue), ])
  systime2           <- Sys.time()
  print(sprintf('End step4: slingshot pseudotime GAM gene calculation at %s', systime2 ))
  print(sprintf("Step4 GAM slingshot pseudotime analysis used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  ## ---
  if (resSave) save(sceObj, slingshotRes, sceObjfitGam, sceObjfitGamATres, file = sprintf('%s_fnPseudoTimeRes.Rdata', resFnamePrefix) )
  ## ------
  return(list('sceObj'=sceObj, 'slingshotRes' = slingshotRes, 'rankGene'=gamSortedGenes))
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
