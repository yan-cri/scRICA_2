## calcPCApseudo(): used to perform functional pseudotime analysis                      ##
## Developed by Yan Li, May, 2021
##--------------------------------------------------------------------------------------##
#' calcPCApseudo() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param obj required, either a seurat or SingleCellExperiment object
#' @param slingshotclusterLabels optional, default no cluster used for slingshot
#' @param resSave optional, default is F
#' @param resFnamePrefix optional, default is 'TEST', when resSave is on
#'
#' @keywords calcPCApseudo
#' @examples calcPCApseudo()
#' @export
#'
#' @return
#' a SingleCellExperiment where 3 methods functional pseudotime analysis results are saved in colData
## ------------------------------------------------------------------------------------ ##
calcPCApseudo <- function(obj, slingshotclusterLabels, resSave, resFnamePrefix) {
  if ( missing(slingshotclusterLabels) ) slingshotclusterLabels <- NA
  if (missing(resSave)) resSave <- as.logical('F')
  if (missing(resFnamePrefix)) resFnamePrefix <- as.character('TEST')
  ## ------
  ## 0. decide the obj input and change into 'SingleCellExperiment' object
  if (class(obj) == "Seurat") {
    sceObj         <- Seurat::as.SingleCellExperiment( obj )
  } else if (class(obj) == "SingleCellExperiment") {
    sceObj         <- obj
  } else {
    stop("please provide either Seurat or SingleCellExperiment object")
  }
  # print(sceObj)
  ## ---
  ## 1. perform BiocSingular PCA analysis
  print('-=-=-=-=-=-=-=-=-=-=')
  print(sprintf('Start 1: BiocSingular PCA pseudotime calculation at %s',  Sys.time() ))
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
  print(sprintf('End 1: BiocSingular PCA pseudotime calculation at %s',  systime2 ))
  print(sprintf("This step's analysis used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  ## ---
  ## 2. perform Diffusion map pseudotime
  print('-=-=-=-=-=-=-=-=-=-=')
  print(sprintf('Start 2: Diffusion map pseudotime calculation at %s',  Sys.time() ))
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
  print(sprintf('End 2: Diffusion map pseudotime calculation at %s',  systime2 ))
  print(sprintf("This step's analysis used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  ## ---
  ## 3. perform Slingshot map pseudotime
  print('-=-=-=-=-=-=-=-=-=-=')
  print(sprintf('Start 3: slingshot pseudotime calculation at %s', Sys.time() ))
  systime1         <- Sys.time()
  if (is.na(slingshotclusterLabels)) {
    print(sprintf('No cluster labels will be used in slingshot calculation' ))
    slingshotRes   <- slingshot::slingshot(sceObj, reducedDim = 'PCA')
  } else {
    ## Note: expCond level cannot be used with 'clusterLabels'
    print(sprintf('%s cluster labels will be used in slingshot calculation', slingshotclusterLabels ))
    slingshotRes   <- slingshot::slingshot(sceObj, clusterLabels = slingshotclusterLabels, reducedDim = 'PCA')
  }
  sceObj$slingPseudotime_1 <-  slingshotRes$slingPseudotime_1
  sceObj$slingPseudotime_2 <-  slingshotRes$slingPseudotime_2
  systime2         <- Sys.time()
  print(sprintf('End 3: slingshot pseudotime calculation at %s', systime2 ))
  print(sprintf("This step's analysis used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  ## ---
  if (resSave) save(sceObj, file = sprintf('%s_fnPseudoTimeRes.Rdata') )
  ## ------
  return(sceObj)
}
## -----------END-----------END-----------END-----------END-----------END-------------- ##
