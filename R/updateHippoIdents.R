#' updateHippoIdents() Function
#' @details
#' This function is used to update seurat object's idents based on provided hippo clustering results.
#'
#' @param resDir full path of integration results analysis are saved, where RDS file is saved inside the 'RDS_Dir'. This path is also returned by getClusterMarkers() execution.
#' @param rds User also can provide the full path of RDS file instead of 'resDir' where RDS file is saved in. If this option is used, please also provide 'resDir' to specify where the analysis results will be saved.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param resSaveFname whether to save the seurat object with the updated idents, by default it is 'NULL' without saving idents updated seurat object.
#'
#' @importFrom Seurat SingleRasterMap
#' @importFrom Seurat Idents
#' @importFrom SeuratObject DefaultAssay
#'
#' @keywords updateHippoIdents
#' @examples updateHippoIdents(rds, hippoResList)
#' @export
#' @return seurat object with updated idents from hippo results
#'
## ------------------------------------------------------------------------------------ ##
updateHippoIdents <- function(resDir = NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, resSaveFname = NULL, hippoResList, hippoResNames = NULL, hippoResK, lightHippo = T) {
  ## ----------------------------------------------------- ##
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  ## ---
  if (is.null(resDir) & !is.null(rds)) {
    if (class(rds)=='Seurat') {
      seuratObjFinal      <<- rds
      print('RDS is provided with rds option')
    } else {
      rdsFname            <- rds
      ## ---
      if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
      seuratObjFinal      <<- readRDS(file = as.character(rdsFname))
      print('Done for RDS read in')
    }
    resDir                <- getwd()
  } else if (!is.null(resDir) & is.null(rds)) {
    rdsFname              <- sprintf('%s/RDS_Dir/%s.rds', resDir, basename(resDir))
    resDir                <- resDir
    ## ---
    if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
    seuratObjFinal          <<- readRDS(file = as.character(rdsFname))
    print('Done for RDS read in')
  } else {
    stop("Error: please provide either option 'resDir' or 'rdsFname'. ")
  }
  rdsDir                <- sprintf('%s/RDS_Dir', resDir)
  if (!dir.exists(rdsDir)) dir.create(rdsDir)
  ## ----------------------------------------------------- ##
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
    print('-=-=-=-')
    print('updated Idents are as below:')
    print(table(Idents(seuratObjFinal)))
    print('-=-=-=-')
  }
  ## ----------------------------------------------------- ##
  seuratObjHippo                  <- seuratObjFinal
  seuratObjHippo$seurat_clusters2 = as.character(Seurat::Idents(seuratObjHippo) )
  print(table(seuratObjHippo$seurat_clusters2))
  print(table(seuratObjHippo$seurat_clusters))
  ## ---
  if (is.null(names(hippoResList)) & is.null(hippoResNames)) {
    hippoResNames <- paste('hippoRes', 1:length(hippoResList), sep = '')
  } else if (!is.null(names(hippoResList)) & is.null(hippoResNames)){
    hippoResNames <- names(hippoResList)
  }
  ## ---
  for (h in 1:length(hippoResList)) {
    print(sprintf("%s. processing '%s' at '%s'.", h, names(hippoResList)[h], hippoResList[[h]]))
    seuratObjHippo$seurat_clusters2 <- as.character(seuratObjHippo$seurat_clusters2)
    load(hippoResList[[h]])
    if (lightHippo) {
      lightHippoCluserRes  <- lightHippo::cut_hierarchy(cluster_res = lightHippoRes, K = hippoResK[h])
      hippo.cluster.update <- paste(as.character(hippoResNames[h]), lightHippoCluserRes, sep = '')
      print(table(hippo.cluster.update))
      if (length(hippo.cluster.update) != length(colnames(inputData))) stop('updated hippo cluster does not match original clusters')
      seuratObjHippo$seurat_clusters2[match(colnames(inputData), rownames(seuratObjHippo@meta.data))] = hippo.cluster.update
    } else {
      seuratObjHippo$seurat_clusters2[match(colnames(inputData@int_metadata$hippo$X), rownames(seuratObjHippo@meta.data))] = paste(as.character(hippoResNames), inputData@int_metadata$hippo$labelmatrix[,dim(inputData@int_metadata$hippo$labelmatrix)[2]], sep = '')
    }
    seuratObjHippo$seurat_clusters2 <- as.factor(seuratObjHippo$seurat_clusters2)
    seuratObjHippo$seurat_clusters  = seuratObjHippo$seurat_clusters2
    Seurat::Idents(seuratObjHippo) = seuratObjHippo$seurat_clusters
    print(table(Seurat::Idents(seuratObjHippo)))
    print('-=-=-=-=-=-=-')
  }
  if (!is.null(resSaveFname)) saveRDS(object = seuratObjHippo, file = file.path(rdsDir, sprintf('%s.rds', resSaveFname)))
  ## ----------------------------------------------------- ##
  return(seuratObjHippo)
  ## ----------------------------------------------------- ##
}
## ------------------------------------------------------------------------------------ ##
