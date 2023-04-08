## getClusterMarkers(): integrate Seurat objects with identified anchor via Seurat method   ##
## Developed by Yan Li, Jan, 2021                                                         ##
##----------------------------------------------------------------------------------------##
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(gridExtra)
## -------------------------------------------------------------------------------------- ##
#' getClusterMarkers() Function
#' @details
#' This function is used to integrate Seurat objects with identified anchor via Seurat method,
#' followed with clustering analysis, and cluster markers identification.
#' @param qcProcessedResults pre-processed QC results object by processQC().
#' @param integrationMethod specify integration method applied for anchor identification for integration, 2 options, 'CCA' and 'RPCA', by default = 'CCA'.
#' @param topN specify the number of identified over expressed genes from each identified cell clusters, by default = 10.
#' @param resDirName specify the folder/directory name where integration analysis results will be saved, by default: the save directory where QC results are saved.
#' @param ribo logical option (True or False), if True, rRNA genes are removed from the variable features for integration analysis.
#' @param int.k.weight specify the number of neighbors to consider when weighting anchors in the step of integration, by default = 100, if too few number cells used in the study, this number can be reduced.
#' If processQC() full results used for option 'qcProcessedResults', it will use the same 'resDirName' used in processQC()
#'
#' @importFrom ggplot2 theme
#' @importFrom gridExtra grid.arrange
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat FindIntegrationAnchors
#' @importFrom Seurat IntegrateData
#' @importFrom Seurat SelectIntegrationFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat ElbowPlot
#' @importFrom Seurat RunUMAP
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat RunTSNE
#' @importFrom Seurat DimPlot
#' @importFrom Seurat FindAllMarkers
#' @importFrom Seurat DoHeatmap
#' @importFrom Seurat NoLegend
#' @importFrom dplyr %>%
#' @importFrom dplyr top_n
#' @importFrom dplyr group_by
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom stats filter
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom utils write.table
#'
#' @keywords getClusterMarkers, seuratIntegrate
#' @examples getClusterMarkers(qcProcessedResults)
#' @export
#' @return
#' the entire integration analysis results in the defined 'resDirName' inherited from 'processQC()'.
#' Additionally, a list item including 3 elements:
#' 1. 'integratedObj': integrated Seurat object;
#' 2. 'posMarkers': identified positively expressed cluster gene markers;
#' 3. 'resDir': full path of results directory, where the entire integration analysis results are saved.
#'
##----------------------------------------------------------------------------------------
getClusterMarkers <- function(qcProcessedResults, integrationMethod = 'CCA', nfeatures = 2000, resDirName = NULL, topN = 10, ribo = F, int.k.weight = 100) {
  ## ---
  topN                          <- as.numeric(topN)
  if ('resDir' %in% names(qcProcessedResults) & 'countReadInOjb' %in% names(qcProcessedResults) ){
    qcProcessedSeuratObjList    <- qcProcessedResults$qcProcessObj
    resDir                      <- qcProcessedResults$resDir
    if (!is.null(resDirName)) stop("Full processQC() results provided in 'qcProcessedResults', no need for option 'resDirName', please remove it")
    if (!dir.exists(resDir)) stop('Error: please make sure QC analysis has been processed correctly and successfully with function processQC()')
  } else {
    qcProcessedSeuratObjList    <- qcProcessedResults
    if (is.null(resDirName)) resDirName = 'integration_results'
    resDir                      <- paste(getwd(), resDirName, sep = '/')
    if (!dir.exists(resDir)) dir.create(resDir)
  }
  resDirName                  <- basename(resDir)
  ##--------------------------------------------------------------------------------------##
  ## intermediate RDS result dir
  rdsDir                        <- paste(resDir, 'RDS_Dir', sep = '/')
  if (!dir.exists(rdsDir)) dir.create(rdsDir)
  print(sprintf('Integration results will be saved in %s', resDir))
  ##--------------------------------------------------------------------------------------##
  theme1noLegend       <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                                # legend.key.size = unit(0.7, "cm"),
                                axis.title = element_text(size = 20),
                                axis.text = element_text(size = 25),
                                legend.position="bottom",
                                legend.text = element_text(size = 14) ) + NoLegend()
  theme1wLegend        <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                                # legend.key.size = unit(0.7, "cm"),
                                axis.title = element_text(size = 20),
                                axis.text = element_text(size = 25),
                                legend.position="bottom",
                                legend.text = element_text(size = 14) )
  ##--------------------------------------------------------------------------------------##
  if (length(qcProcessedSeuratObjList) == 1) {
    ## ---
    print("No integration is needed, only 1 Seurat object is provided.")
    print(sprintf('Start: Step 2 tSNE and UMP clustering at %s.', Sys.time()))
    seuratObjIntegrated         <-  qcProcessedSeuratObjList[[1]]
  } else {
    ## 3. Integration
    ##--------------------------------------------------------------------------------------##
    ## 3.1 Finding anchors, by default running CCA
    print(sprintf('Step 1: Process data integration at %s', Sys.time()))
    print(sprintf('%s samples will be integrated', length(qcProcessedSeuratObjList)))
    print(sprintf("%s anchor integration is implemented on %s features.", integrationMethod, nfeatures))
    ## if provided 'qcProcessedSeuratObjList' has not implemented with 'FindVariableFeatures', do it here with normalization again
    if (any(unlist(lapply(qcProcessedSeuratObjList, function(x) length(Seurat::VariableFeatures(x))==0)))) {
      print("Note: some input data in 'qcProcessedResults' seems not to be normalized, before integration, conducing normalization again")
      qcProcessedSeuratObjList <- lapply(X = qcProcessedSeuratObjList, FUN = function(x) {
        x <- Seurat::NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
        x <- Seurat::FindVariableFeatures(x, selection.method = 'vst', nfeatures = nfeatures)
        if (ribo == TRUE) {
          Seurat::VariableFeatures(x) <- Seurat::VariableFeatures(x)[-grep("^RP[SL]", Seurat::VariableFeatures(x))]
        }
      })
    } else if (any(unlist(lapply(qcProcessedSeuratObjList, function(x) length(Seurat::VariableFeatures(x))<nfeatures)))) {
      print("Note: some input data in 'qcProcessedResults' has fewer number of detected varialbe features requested in 'nfeatures' for integration, before integration, conducing normalization again")
      qcProcessedSeuratObjList <- lapply(X = qcProcessedSeuratObjList, FUN = function(x) {
        x <- Seurat::NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
        x <- Seurat::FindVariableFeatures(x, selection.method = 'vst', nfeatures = nfeatures)
        if (ribo == TRUE) {
          Seurat::VariableFeatures(x) <- Seurat::VariableFeatures(x)[-grep("^RP[SL]", Seurat::VariableFeatures(x))]
        }
      })
    } else {
      print("Note: input data in 'qcProcessedResults' has been appropriately normalized, no need to normalize before integration.")
    }
    ##--------------------------------------------------------------------------------------##
    features                    <- SelectIntegrationFeatures(object.list = qcProcessedSeuratObjList, nfeatures = nfeatures)

    ##Remove rRNA content from the variable features prior to integration, if desired
    if (ribo == TRUE) {
        features = features[-grep("^RP[SL]", features)]
      }

    if (integrationMethod == 'RPCA') {
      # normalize and identify variable features for each dataset independently
      qcProcessedSeuratObjList <- lapply(X = qcProcessedSeuratObjList, FUN = function(x) {
        x <- ScaleData(object = x)
        x <- RunPCA(object = x)
      })

      anchors                     <- FindIntegrationAnchors(object.list = qcProcessedSeuratObjList,
                                                            anchor.features = features,
                                                            reduction = 'rpca' )

    } else {
      anchors                     <- FindIntegrationAnchors(object.list = qcProcessedSeuratObjList,
                                                            anchor.features = features )
    }
    ##--------------------------------------------------------------------------------------##
    ## 3.2 Integrating data based on found anchors above
    seuratObjIntegrated         <- IntegrateData(anchorset = anchors, k.weight = int.k.weight)
    print(sprintf('End Step 1: data integration at %s', Sys.time()))
    print('---===---===---===---===---===---')
    ##--------------------------------------------------------------------------------------##
    print(sprintf('Start: Step 2 tSNE and UMP clustering at %s.', Sys.time()))
    ## 4. dimensional reduction with PCA, KNN/clusters, tsne/umap
    ## 4.1 scale top 2000 identified variable features with default liner model in 'model.use' option
    ##     The results of this are stored in seuratObjFinal[["RNA"]]@scale.data for sep & seuratObjFinal[["integrated"]]@scale.data for integrated data
    Seurat::DefaultAssay(seuratObjIntegrated) <- "integrated"
  }
  ##--------------------------------------------------------------------------------------##
  seuratObjFinal                <- ScaleData(object = seuratObjIntegrated)
  ##--------------------------------------------------------------------------------------##
  ## 4.2 PCA
  seuratObjFinal                <- RunPCA(object = seuratObjFinal)
  ## 4.2.2 determine the 'dimensinality' of the dataset
  # seuratObjFinal <- JackStraw(seuratObjFinal, num.replicate = 100)
  # seuratObjFinal <- ScoreJackStraw(seuratObjFinal, dims = 1:20)
  elbowPlot                     <- ElbowPlot(seuratObjFinal, ndims = 30)
  ## 4.3 umap, tsne, knn/clusters
  seuratObjFinal                <- RunUMAP(seuratObjFinal, reduction = "pca", dims = 1:20)
  seuratObjFinal                <- FindNeighbors(seuratObjFinal, reduction = "pca", dims = 1:20)
  seuratObjFinal                <- FindClusters(seuratObjFinal, resolution = 0.5)
  seuratObjFinal                <- RunTSNE(object = seuratObjFinal, dims = 1:20)
  ##--------------------------------------------------------------------------------------##
  ## PCA plot
  pcaCluster  <- DimPlot(seuratObjFinal, reduction = "pca") + labs(title = 'PCA clustering', x = "PC 1", y = 'PC 2')
  pcaExpCond  <- DimPlot(seuratObjFinal, reduction = "pca", group.by = 'expCond') + labs(title = 'PCA clustering', x = "PC 1", y = 'PC2')
  ## -
  pdf(file = file.path(resDir, 'pca_plot.pdf'), width = 10, height = 6)
  grid.arrange(pcaCluster + theme1wLegend, pcaExpCond + theme1wLegend, ncol = 2, nrow = 1)
  dev.off()
  ##--------------------------------------------------------------------------------------##
  ## UMAP plot
  umapCluster <- DimPlot(seuratObjFinal, reduction = "umap", label = T, repel = T) + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
  umapExpCond <- DimPlot(seuratObjFinal, reduction = "umap", group.by = 'expCond') + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
  umapSplit   <- DimPlot(seuratObjFinal, reduction = "umap", split.by = 'expCond') + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
  ## -
  pdf(file = file.path(resDir, 'umap_plot.pdf'), width = 10, height = 6)
  grid.arrange(umapCluster + theme1wLegend, umapExpCond + theme1wLegend, nrow = 1, ncol = 2)
  # grid.arrange(umapSplit + theme1noLegend,
  #              arrangeGrob(umapCluster + theme1wLegend, umapExpCond + theme1wLegend, nrow = 1, ncol = 2),
  #              nrow = 2, ncol = 1)
  dev.off()
  pdf(file = file.path(resDir, 'umap_plot_samplSep.pdf'), width = 4.5*length(qcProcessedSeuratObjList), height = 6)
  grid.arrange(umapSplit + theme1noLegend)
  dev.off()
  ##--------------------------------------------------------------------------------------##
  ## tsne plot
  tsneExpCond <- DimPlot(seuratObjFinal, reduction = "tsne", group.by = 'expCond') + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
  tsneCluster <- DimPlot(seuratObjFinal, reduction = "tsne", label = T, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
  tsneSplit   <- DimPlot(seuratObjFinal, reduction = "tsne", split.by = 'expCond') + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
  ## -
  pdf(file = file.path(resDir, 'tsne_plot.pdf'), width = 10, height = 6)
  grid.arrange(tsneCluster + theme1wLegend, tsneExpCond + theme1wLegend, nrow = 1, ncol = 2)
  # grid.arrange(tsneSplit + theme1noLegend,
  #              arrangeGrob(tsneCluster + theme1wLegend, tsneExpCond + theme1wLegend, nrow = 1, ncol = 2),
  #              nrow = 2, ncol = 1)
  dev.off()
  pdf(file = file.path(resDir, 'tsne_plot_samplSep.pdf'), width = 4.5*length(qcProcessedSeuratObjList), height = 6)
  print(tsneSplit + theme1noLegend)
  dev.off()
  print(sprintf('END: Step 2 tSNE and UMP clustering at %s.', Sys.time()))
  print('---===---')
  print(seuratObjFinal)
  saveRDS(seuratObjFinal, file = file.path(rdsDir, sprintf("%s.rds", resDirName)) )
  print('Complete save RDS')
  print('---===---')
  ##--------------------------------------------------------------------------------------##
  ## 5. 1. Find all marker
  print('Start: Step 3 finding positive regulated cluster marker genes')
  allPosMarkers        <- FindAllMarkers(seuratObjFinal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(x = allPosMarkers, file = file.path(resDir, 'allCluster_pos_markers.txt'), quote = F, sep = '\t', row.names = T, col.names = NA)
  Sys.time()
  ## -
  noClusterMarkers     <- as.data.frame(table(allPosMarkers$cluster))
  write.table(x = noClusterMarkers, file = file.path(resDir, 'allCluster_pos_markers_no.txt'), quote = F, sep = '\t', row.names = T, col.names = NA)
  print('END: Step 3 finding positive regulated cluster marker genes')
  print('---===---')
  ## -
  print('Start: Step 4 making cluster marker genes heatmap plot')
  top10                 <- allPosMarkers %>% group_by(cluster) %>% top_n(n = topN, wt = avg_log2FC) %>% as.data.frame()
  write.table(x = top10, file = file.path(sprintf('%s/allCluster_pos_markers_top%s.txt', resDir, topN)), quote = F, sep = '\t', row.names = T, col.names = NA)
  cluterTop10heatmap    <- DoHeatmap(seuratObjFinal, features = top10$gene) + NoLegend()
  pdf(file = file.path(resDir, sprintf('cluster_heatmap_top%sPosMarkers.pdf', topN)), width = 25, height = 20)
  print(cluterTop10heatmap)
  dev.off()
  ##--------------------------------------------------------------------------------------##
  print('END: Step 4 making cluster marker genes heatmap plot')
  print('END===END===END')
  ##--------------------------------------------------------------------------------------##
  return(list('resDir' = resDir, 'integratedObj' = seuratObjFinal, 'posMarkers' = allPosMarkers))
}
##----------------------------------------------------------------------------------------
