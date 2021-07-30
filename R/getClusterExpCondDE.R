#' getClusterExpCondDe() Function
#' @details
#' This function is used to identify DEGs in each identified/annotated cell cluster based on the experimental conditions from metadata table column 'sample', 'expCond1', or 'expCond2'
#'
#' @param resDir full path of integration results analysis are saved, where RDS file is saved inside the 'RDS_Dir'. This path is also returned by getClusterMarkers() execution.
#' @param rdsFname User also can provide the full path of RDS file instead of 'resDir' where RDS file is saved in. If this option is used, please also provide 'resDir' to specify where the analysis results will be saved.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param expCondCheck 3 options: 'sample', 'expCond1', or 'expCond2' to specify which experimental conditions to be explored with this function.
#' @param expCondSepName part of file name string to specify the analysis results folder name.
#' @param expCondName2change character string to indicate part of characters specified here can be removed from sample name defined in the metadata table, if additional samples combination needs to be explored which has not been specified in the column of 'expCond1' or 'expCond2'.
#' @param compGroup specify 2 group names for comparisons, these group names is corresponding to experimental condition factor levels.
#' @param deMethod DE test method with options: 'wilcox', 't', 'negbinom', 'poisson', 'MAST', 'DESeq2', by default = 'wilcox'.
#' @param pAdjValCutoff adjusted p-value cutoff for significant DEGs detection, by default = 0.05.
#' @param topNo specify the top number of up and down significant DEGs in each identified/annotated cell clusters, by default = 10.
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @importFrom Seurat RenameIdents
#' @importFrom Seurat FindMarkers
#' @importFrom utils write.table
#' @importFrom xlsx write.xlsx
#'
#' @keywords FindMarkers
#' @examples getClusterExpCondDe()
#' @export
#' @return
#' a list item including 3 elements:
#' 1: 'clusterDeMarkers'
#' 2. 'clusterDeResSummary'
#' 3. 'clusterTopDeMarkers'
## ---------------------------------------------------------------------------------------
getClusterExpCondDe <- function(resDir=NULL, rdsFname=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, expCondCheck='sample', expCondSepName = NULL, expCondName2change=NULL, compGroup, deMethod = 'wilcox', pAdjValCutoff = 0.05, topNo = 10) {
  options(java.parameters = "-Xmx32000m")
  ## ---
  if (missing(compGroup)) stop("Please provide option 'compGroup' to specify which 2 groups in your updated experimental condition levels with 'expCondCheck' for comparision.")
  pAdjValCutoff           <- as.numeric(pAdjValCutoff)
  topNo                   <- as.numeric(topNo)
  deMethod                <- as.character(deMethod)
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) print("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  ## ---
  if (is.null(resDir) & !is.null(rdsFname)) {
    rdsFname              <- rdsFname
    resDir                <- getwd()
  } else if (!is.null(resDir) & is.null(rdsFname)) {
    rdsFname              <- sprintf('%s/RDS_Dir/%s.rds', resDir, basename(resDir))
    resDir                <- resDir
  } else {
    stop("Error: please provide either option 'resDir' or 'rdsFname'. ")
  }
  ## ---
  if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
  seuratObjFinal          <<- readRDS(file = as.character(rdsFname))
  print('Done for RDS readin')
  ## ------
  ## ------
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir <- paste(resDir, 'results_wNewAnnotation', sep = '/')
  } else {
    resDir <- paste(resDir, 'results_wOrgClusterAnnotation', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('DEG identification on cell clusters will be saved at %s', resDir))
  ## -------------------------------------------------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  }
  ## -------------------------------------------------------------------------------------
  if (expCondCheck == 'sample') {
    if (is.null(expCondSepName)) {
      expCondSepName        <- 'expCond_sample'
    } else {
      expCondSepName        <- expCondSepName
    }
  } else {
    if (is.null(expCondSepName)) {
      expCondSepName        <- as.character(expCondCheck)
    } else {
      expCondSepName        <- expCondSepName
    }
  }
  ## update 'resDir' to ceate dir under 'result_wNewAnnotation' or 'results_wOrgClusterAnnotation'
  resDir                  <- paste(sprintf('%s/DEG_%s_comp_%s', resDir, expCondSepName, deMethod))
  if (!dir.exists(resDir)) dir.create(resDir)
  ## -------------------------------------------------------------------------------------
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding updated 'resDir' for new tSNE/UMAP plots to save
  resDir                  <- paste(resDir, sprintf('DEG_%s_comp_%s_%s', expCondSepName, gsub('/', '-', compGroup), deMethod ), sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  ## -------------------------------------------------------------------------------------
  if (expCondCheck == 'sample') {
    seuratObjFinal                     <- seuratObjFinal
  } else if (expCondCheck == 'expCond1') {
    if (!'expCond1' %in% colnames(seuratObjFinal@meta.data)){
      print("Error: 'expCond1' has not been included in the original integration analysis.")
      seuratObjFinal@meta.data$expCond <- gsub(pattern = as.character(expCondName2change), replacement = '', x = seuratObjFinal@meta.data$expCond)
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data$expCond1
    }
  } else if (expCondCheck == 'expCond2') {
    if (!'expCond2' %in% colnames(seuratObjFinal@meta.data)){
      print("Error: 'expCond2' has not been included in the original integration analysis.")
      seuratObjFinal@meta.data$expCond <- gsub(pattern = as.character(expCondName2change), replacement = '', x = seuratObjFinal@meta.data$expCond)
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data$expCond2
    }
  }
  ## -------------------------------------------------------------------------------------
  # Seurat::DefaultAssay(seuratObjFinal)   <- "RNA"
  Seurat::DefaultAssay(seuratObjFinal) <- "integrated"
  ## ---
  noClusters                             <- levels(Idents(seuratObjFinal))
  print(sprintf('Start step 1: identifing DEGs for each identified clusters based on sample name experimental conditions'))
  print(sprintf("A total of %s clusters to process", length(noClusters)))
  ## -
  clusterDeMarkers                       <- list()
  clusterDeResSummary                    <- matrix(data = NA, nrow = length(noClusters), ncol = 4)
  clusterTopDeMarkers                    <- list()
  clusterTopDeMarkersUp                  <- NA
  clusterTopDeMarkersUpCluster           <- NA
  clusterTopDeMarkersDown                <- NA
  clusterTopDeMarkersDownCluster         <- NA
  clusterDeMarkersFname                  <- paste(sprintf('%s/expCondCompDeMarkers_%s', resDir, gsub('/', '-', compGroup) ))
  ## -
  for (i in 1:length(noClusters)) {
  # for (i in c(1,3,4)) {
    print(sprintf('Start 1.%s processing cluster %s for DE markers identification', i, noClusters[i]))
    seuratObjFinalSubet                  <- subset(seuratObjFinal, idents = as.character(noClusters[i]) )
    print("-=-=-=")
    print("experimental conditions:")
    print(table(seuratObjFinalSubet$expCond))
    print("-=-=-=")
    if ( all(table(seuratObjFinalSubet$expCond) >= 3) ) {
      compGroupName                      <- strsplit(compGroup, split = '/')[[1]]
      if ( all(compGroupName %in% names(table(seuratObjFinalSubet$expCond))) ) {
        ## ---
        if ( length(compGroupName)>2 ) {
          stop("Too many comparision groups provided in 'compGroup', only 2 names shoulb be provided ")
        } else {
          if (length(compGroupName) == 1) {
            deMarkers                    <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName, group.by = 'expCond', test.use = deMethod, min.cells.group = 3, verbose = T)
          } else if (length(compGroupName) == 2) {
            deMarkers                    <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName[1], ident.2 = compGroupName[2], group.by = 'expCond', test.use = deMethod, min.cells.group = 3, verbose = T)
          }
          clusterDeMarkers[[i]]          <- deMarkers
          print(sprintf('Maximum p_value is %s, Maximum adjusted p_value is %s', round(max(deMarkers$p_val), digits = 4), round(max(deMarkers$p_val_adj), digits = 4)))
          deMarkersAdjSig                <- deMarkers %>% dplyr::filter(p_val_adj <= pAdjValCutoff) %>% dplyr::mutate(perDiff = pct.1-pct.2)
          deMarkersAdjSigUp              <- deMarkersAdjSig %>% dplyr::filter(avg_log2FC > 0) %>% dplyr::arrange(desc(perDiff))
          deMarkersAdjSigDown            <- deMarkersAdjSig %>% dplyr::filter(avg_log2FC < 0) %>% dplyr::arrange(desc(perDiff))
          clusterTopDeMarkers[[i]]       <- list('up' = rownames(deMarkersAdjSigUp)[1:topNo],
                                                 'down' = rownames(deMarkersAdjSigDown)[1:topNo])
          clusterTopDeMarkersUp          <- c(clusterTopDeMarkersUp, rownames(deMarkersAdjSigUp)[1:topNo])
          clusterTopDeMarkersUpCluster   <- c(clusterTopDeMarkersUpCluster, rep(noClusters[i], length(rownames(deMarkersAdjSigUp)[1:topNo])) )
          clusterTopDeMarkersDown        <- c(clusterTopDeMarkersDown, rownames(deMarkersAdjSigDown)[1:topNo])
          clusterTopDeMarkersDownCluster <- c(clusterTopDeMarkersDownCluster, rep(noClusters[i], length(rownames(deMarkersAdjSigDown)[1:topNo])))
          clusterDeResSummary[i,]        <- c(dim(deMarkers)[1], dim(deMarkersAdjSig)[1], dim(deMarkersAdjSigUp)[1], dim(deMarkersAdjSigDown)[1] )
          print(sprintf('out of %s DE markers, %s are significantly DE at adjusted p-value of %s', dim(deMarkers)[1], dim(deMarkersAdjSig)[1], pAdjValCutoff))
          print(sprintf('out of %s significantly DE markers, %s are positively expressed, %s are negative expressed.', dim(deMarkersAdjSig)[1], dim(deMarkersAdjSigUp)[1], dim(deMarkersAdjSigDown)[1] ))
          ## -
          if (newAnnotation) {
            resFname1                  <- paste(clusterDeMarkersFname, '_full_wNewAnnotation.xlsx', sep = '')
            resFname2                  <- paste(clusterDeMarkersFname, '_adjSig_wNewAnnotation.xlsx', sep = '')
            resFname3                  <- paste(clusterDeMarkersFname, '_adjSig_up_wNewAnnotation.xlsx', sep = '')
            resFname4                  <- paste(clusterDeMarkersFname, '_adjSig_down_wNewAnnotation.xlsx', sep = '')
          } else {
            resFname1                  <- paste(clusterDeMarkersFname, '_full.xlsx', sep = '')
            resFname2                  <- paste(clusterDeMarkersFname, '_adjSig.xlsx', sep = '')
            resFname3                  <- paste(clusterDeMarkersFname, '_adjSig_up.xlsx', sep = '')
            resFname4                  <- paste(clusterDeMarkersFname, '_adjSig_down.xlsx', sep = '')
          }
          ## -
          if (i ==1) {
            if (dim(deMarkers)[1] > 0) write.xlsx(x = deMarkers, file = resFname1, sheetName = paste('cluster', gsub('/|:', ' ', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSig)[1] > 0) write.xlsx(x = deMarkersAdjSig, file = resFname2, sheetName = paste('cluster', gsub('/|:', ' ', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSigUp)[1] > 0) write.xlsx(x = deMarkersAdjSigUp, file = resFname3, sheetName = paste('cluster', gsub('/|:', ' ', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSigDown)[1] > 0) write.xlsx(x = deMarkersAdjSigDown, file = resFname4, sheetName = paste('cluster', gsub('/|:', ' ', noClusters[i]), sep = '_'), row.names = T, append = F )
          } else {
            if (dim(deMarkers)[1] > 0) write.xlsx(x = deMarkers, file = resFname1, sheetName = paste('cluster', gsub('/|:', ' ', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSig)[1] > 0)write.xlsx(x = deMarkersAdjSig, file = resFname2, sheetName = paste('cluster', gsub('/|:', ' ', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSigUp)[1] > 0) write.xlsx(x = deMarkersAdjSigUp, file = resFname3, sheetName = paste('cluster', gsub('/|:', ' ', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSigDown)[1] > 0) write.xlsx(x = deMarkersAdjSigDown, file = resFname4, sheetName = paste('cluster', gsub('/|:', ' ', noClusters[i]), sep = '_'), row.names = T, append = T )
          }
        }
        ## ---
      } else {
        stop("Provided comparision group names in 'compGroup' does not match experimental conditions ")
      }
    } else {
      print(sprintf('cluster %s has less than 3 cells in one experimental condiction, no DE markers can be identfied for this cluster', noClusters[i]))
    }
    print(sprintf('End 1.%s processing cluster %s for DE markers identification', i, noClusters[i]))
    print('=========')
  }
  ## -
  clusterTopDeMarkersUpComb          <- data.frame('geneName' = clusterTopDeMarkersUp[-1], 'cluster'= clusterTopDeMarkersUpCluster[-1])
  clusterTopDeMarkersDownComb        <- data.frame('geneName' = clusterTopDeMarkersDown[-1], 'cluster' = clusterTopDeMarkersDownCluster[-1])
  write.xlsx( x = clusterTopDeMarkersUpComb[!is.na(clusterTopDeMarkersUpComb$geneName),], file = sprintf('%s_top%s_upDe.xlsx', clusterDeMarkersFname, topNo), sheetName = 'up', row.names = F, col.names = T, append = F)
  write.xlsx( x = clusterTopDeMarkersDownComb[!is.na(clusterTopDeMarkersDownComb$geneName),], file = sprintf('%s_top%s_downDe.xlsx', clusterDeMarkersFname, topNo), sheetName = 'down', row.names = F, col.names = T, append = F)
  ## -
  rownames(clusterDeResSummary)      <- noClusters
  colnames(clusterDeResSummary)      <- c('Total', 'DE', 'Up', 'Down')
  clusterDeResSummaryFname           <- paste(clusterDeMarkersFname, '_NoDEmarkers_summary.xlsx', sep = '')
  write.xlsx( x = clusterDeResSummary, file = clusterDeResSummaryFname, sheetName = 'No. summary', row.names = T, col.names = T, append = F)
  ## ---
  print(sprintf('End step 1: identifing DEGs for each identified clusters based on sample name experimental conditions'))
  print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
  return(list('clusterDeMarkers' = clusterDeMarkers, 'clusterDeResSummary' = clusterDeResSummary, 'clusterTopDeMarkers' = clusterTopDeMarkers ))
  ## -------------------------------------------------------------------------------------
}
