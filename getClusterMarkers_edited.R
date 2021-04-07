## seuratIntegrate(): integrate seurat objects with identified anchor via seurat method   ##
## Developed by Yan Li, Jan, 2021 
## Edited by Michiko Ryu, April, 2021
##----------------------------------------------------------------------------------------##
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(gridExtra)
##--------------------------------------------------------------------------------------##
#' getClusterMarkers() Function
#' @details
#' This function is used to integrate seurat objects with identified anchor via seurat method,
#' followed with conducing cluster analysis, and identifying cluster markers
#' @param qcProcessedSeuratObjList returned 'qcProcessObj' from countReadin() function.
#' @param anchorIntegrate whehter to implement anchorIntegrate, default = T.
#' @param resDirName optional, define folder/directory name where integration analysis results will be saved
#'
#' @importFrom ggplot2 theme
#' @importFrom gridExtra grid.arrange
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat FindIntegrationAnchors
#' @importFrom Seurat IntegrateData
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
#' @examples getClusterMarkers()
#' @export
#' @return
#' a list item including 3 elements:
#' 1. 'integratedObj': finalized integrated seurat object;
#' 2. 'posMarkers': identified positively expressed cluster markers;
#' 3. 'resDir': full path of results directory, where includes Integration results, pca plot, umap plot, tsne plot, positive regulated cluster marker genes and their  heatmap plot.
##----------------------------------------------------------------------------------------
getClusterMarkers <- function(qcProcessedSeuratObjList, anchorIntegrate, resDirName) {
  if (missing(anchorIntegrate)) anchorIntegrate <- as.logical(T)
  if (length(qcProcessedSeuratObjList) == 1 & anchorIntegrate == as.logical(T)) stop("ERROR: only 1 item in input 'qcProcessedSeuratObjList', no integration can be performed ")
  if (missing(resDirName)) resDirName <- 'integration_analysis_results'
  ## 0. create 'resDir' based on provided 'resDirName' under current workDir
  resDir               <- paste(getwd(), resDirName, sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  ## intermediate RDS result dir
  rdsDir               <- paste(resDir, 'RDS_Dir', sep = '/')
  if (!dir.exists(rdsDir)) dir.create(rdsDir)
  ## setup cutome theme for plotting
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
  ## ------------------
  print(sprintf('Integration results will be saved in %s', resDir))
  ## wether integrate input objecct
  if (anchorIntegrate) {
    ## 3. Integration
    ## 3.1 Finding anchors, by default runing CCA
    print(sprintf('Step 3: Process data integration at %s', Sys.time()))
    print(sprintf('%s samples will be integrated', length(qcProcessedSeuratObjList)))
    anchors                     <- FindIntegrationAnchors(object.list = qcProcessedSeuratObjList)
    ## 3.2 Integrating data based on found anchors above
    seuratObjIntegrated         <- IntegrateData(anchorset = anchors)
    print(sprintf('End Step 3: data integration at %s', Sys.time()))
    print('---===---===---===---===---===---')
    ## ---
    print(sprintf('Start: Step 4 tSNE and UMP clustering at %s.', Sys.time()))
    ## 4. dimentional reduction with PCA, KNN/clusters, tsne/umap
    ## 4.1 scale top 2000 identified variable features with default liner model in 'model.use' option
    ##     The results of this are stored in seuratObjFinal[["RNA"]]@scale.data for sep & seuratObjFinal[["integrated"]]@scale.data for integrated data
    Seurat::DefaultAssay(seuratObjIntegrated) <- "integrated"
  } else {
    ## ---
    print("No integration is needed, only 1 object in the input 'qcProcessedSeuratObjList'.")
    print(sprintf('Start: Step 4 tSNE and UMP clustering at %s.', Sys.time()))
    seuratObjIntegrated         <-  qcProcessedSeuratObjList[[1]]
  }
  ## ---
  seuratObjFinal                <- ScaleData(object = seuratObjIntegrated)
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
  ## PCA plot
  pcaCluster  <- DimPlot(seuratObjFinal, reduction = "pca") + labs(title = 'PCA clustering', x = "PC 1", y = 'PC 2')
  pcaExpCond  <- DimPlot(seuratObjFinal, reduction = "pca", group.by = 'expCond') + labs(title = 'PCA clustering', x = "PC 1", y = 'PC2')
  ## -
  pdf(file = file.path(resDir, 'pca_plot.pdf'), width = 10, height = 6)
  grid.arrange(pcaCluster + theme1wLegend, pcaExpCond + theme1wLegend, ncol = 2, nrow = 1)
  dev.off()
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
  print(sprintf('END: Step 4 tSNE and UMP clustering at %s.', Sys.time()))
  print('---===---')
  print(seuratObjFinal)
  saveRDS(seuratObjFinal, file = file.path(rdsDir, sprintf("%s.rds", resDirName)) )
  print('Complete save RDS')
  print('---===---')
  ## 5. 1. Find all marker
  ## -
  print('Start: Step 5 finding positive regulated cluster marker genes')
  allPosMarkers <- FindAllMarkers(seuratObjFinal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(x = allPosMarkers, file = file.path(resDir, 'allCluster_pos_markers.txt'), quote = F, sep = '\t', row.names = T, col.names = NA)
  Sys.time()
  ## -
  noClusterMarkers <- as.data.frame(table(allPosMarkers$cluster))
  write.table(x = noClusterMarkers, file = file.path(resDir, 'allCluster_pos_markers_no.txt'), quote = F, sep = '\t', row.names = T, col.names = NA)
  print('END: Step 5.1 finding positive regulated cluster marker genes')
  print('---===---')
  ## -
  print('Start: Step 6 making cluster marker genes heatmap plot')
  top10                 <- allPosMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% as.data.frame()
  write.table(x = top10, file = file.path(resDir, 'allCluster_pos_markers_top10.txt'), quote = F, sep = '\t', row.names = T, col.names = NA)
  cluterTop10heatmap    <- DoHeatmap(seuratObjFinal, features = top10$gene) + NoLegend()
  pdf(file = file.path(resDir, 'cluster_heatmap_top10PosMarkers.pdf'), width = 25, height = 20)
  print(cluterTop10heatmap)
  dev.off()
  ## -
  top20                 <- allPosMarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% as.data.frame()
  write.table(x = top20, file = file.path(resDir, 'allCluster_pos_markers_top20.txt'), quote = F, sep = '\t', row.names = T, col.names = NA)
  cluterTop20heatmap    <- DoHeatmap(seuratObjFinal, features = top20$gene) + NoLegend()
  pdf(file = file.path(resDir, 'cluster_heatmap_top20PosMarkers.pdf'), width = 25, height = 20)
  print(cluterTop20heatmap)
  dev.off()
  ## -
  top50                 <- allPosMarkers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) %>% as.data.frame()
  write.table(x = top50, file = file.path(resDir, 'allCluster_pos_markers_top50.txt'), quote = F, sep = '\t', row.names = T, col.names = NA)
  cluterTop50heatmap    <- DoHeatmap(seuratObjFinal, features = top50$gene) + NoLegend()
  pdf(file = file.path(resDir, 'cluster_heatmap_top50PosMarkers.pdf'), width = 25, height = 20)
  print(cluterTop50heatmap)
  dev.off()
  ## -
  top100                <- allPosMarkers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% as.data.frame()
  write.table(x = top100, file = file.path(resDir, 'allCluster_pos_markers_top100.txt'), quote = F, sep = '\t', row.names = T, col.names = NA)
  cluterTop100heatmap   <- DoHeatmap(seuratObjFinal, features = top100$gene) + NoLegend()
  pdf(file = file.path(resDir, 'cluster_heatmap_top100PosMarkers.pdf'), width = 30, height = 25)
  print(cluterTop100heatmap)
  dev.off()
  print('END: Step 6 making cluster marker genes heatmap plot')
  print('END===END===END')
  return(list('resDir' = resDir, 'integratedObj' = seuratObjFinal, 'posMarkers' = allPosMarkers))
}
##----------------------------------------------------------------------------------------
