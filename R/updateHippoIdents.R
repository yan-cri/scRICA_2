#' updateHippoIdents() Function
#' @details
#' This function is used to update seurat object's idents based on provided hippo clustering results.
#'
#' @param resDir specify a exiting full path of integration results analysis are saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param resSaveFname whether to save the seurat object with the updated idents, by default it is 'NULL' without saving idents updated seurat object.
#' @param hippoResK corresponding cluster k number with respect to list input in 'hippoResList'
#' @param hippoResList list of Rdata inputs for cell types to be updated with HIPPO results
#' @param hippoResNames corresponding updated HIPPO results names with respect to list input in 'hippoResList'
#' @param debug whether to turn on for debug check, by default FALSE (turned off).
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
updateHippoIdents <- function(resDir = NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, resSaveFname = NULL, hippoResList, hippoResNames = NULL, hippoResK, lightHippo = T, debug = F) {
  ## ----------------------------------------------------- ##
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  ##--------------------------------------------------------------------------------------##
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
  } else if (is.null(resDir) & is.null(rds)){
    stop("Error: please provide either option 'resDir' or 'rds', or both. ")
  } else if (!is.null(resDir) & !is.null(rds)){
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
    resDir                <- resDir
  }
  ##--------------------------------------------------------------------------------------##
  rdsDir                <- sprintf('%s/RDS_Dir', resDir)
  if (!dir.exists(rdsDir)) dir.create(rdsDir)
  hippoResDir <- sprintf('%s/hippo_results', resDir)
  if (!dir.exists(hippoResDir)) dir.create(hippoResDir)
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
  if (debug) {
    print(table(seuratObjHippo$seurat_clusters2))
    print(table(seuratObjHippo$seurat_clusters))
  }
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
    if(debug) print(table(Seurat::Idents(seuratObjHippo)))
    ## ---------
    # if (lightHippo) {
    #   pdf(file = file.path(hippoResDir, sprintf('%s_hierarchy_k%s.pdf', names(hippoResList)[h], hippoResK[h])), width = 8, height = 5)
    #   newcut.lightHippoRes <- lightHippo::cut_hierarchy(lightHippoRes, K = hippoResK[h], cut_sequence = T)
    #   lightHippo::visualize_hippo_hierarchy(newcut.lightHippoRes)
    #   dev.off()
    # }
    ## ---------
    print('-=-=-=-=-=-=-')
  }
  if (!is.null(resSaveFname)) saveRDS(object = seuratObjHippo, file = file.path(rdsDir, sprintf('%s.rds', resSaveFname)))
  ## ----------------------------------------------------- ##
  return(seuratObjHippo)
  ## ----------------------------------------------------- ##
}
## ------------------------------------------------------------------------------------ ##
