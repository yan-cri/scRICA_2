## This script is used to 1) annotate/assign and merge identified clusters,
## 2) subset certain types of cells to conduct subclustering
## ---------------------------------------------------------------------------------------
# library(Seurat)
# library(ggplot2)
# library(dplyr)
##--------------------------------------------------------------------------------------##
#' clusterSummary() Function
#' @details
#' This function is used to annotate/assign and merge identified clusters,
#'
#' @param resDir to be added
#' @param rdsFname to be added
#' @param clusteringPlotRemake To be added
#' @param clusterCellsNoSummary To be added
#'
#' @keywords clusterSummary
#' @examples clusterSummary()
#' @export
#' @return
#' None
## ---------------------------------------------------------------------------------------

clusterSummary <- function(resDir, rdsFname, clusteringPlotRemake, clusterCellsNoSummary) {
  # resDir                <- paste(workDir, resDirName, sep = '/')
  # seuratObjFinal        <- readRDS(file = file.path('/Users/yanli/Desktop/757_scRNA-seq/', sprintf("%s.rds", resDirName)))
  resDirName            <- strsplit(x = resDir, split = '/')[[1]][length(strsplit(x = resDir, split = '/')[[1]])]
  if (missing(rdsFname)) rdsFname <- paste(resDir, sprintf("%s.rds", resDirName), sep = '/' )
  if (!file.exists(rdsFname)) stop("Please provide 'rdsFname' with full path.")
  seuratObjFinal        <- readRDS(file = as.character(rdsFname))
  if (missing(clusteringPlotRemake)) clusteringPlotRemake  <- as.logical(T)
  if (missing(clusterCellsNoSummary)) clusterCellsNoSummary <- as.logical(T)
  ##  below scripts copied over to 'afterClustering_downstream_dirSetup.R'
  source(paste(orgDir, 'R/afterClustering_downstream_dirSetup.R', sep = '/'))
  ## ---------------------------------------------------------------------------------------
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding updated 'resDir' for new tSNE/UMAP plots to save
  resDir              <- paste(resDir, sprintf('expCond_%s', expCondSep), sep = '/')
  if (expCondSep == 'org') {
    seuratObjFinal    <- seuratObjFinal
  } else {
    seuratObjFinal@meta.data$expCond <- gsub(pattern = as.character(expCondName2change), replacement = '', x = seuratObjFinal@meta.data$expCond)
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  ## ---------------------------------------------------------------------------------------
  ## 1. re-make tSNE plot
  if (clusteringPlotRemake) {
    newResDir          <- paste(resDir, sprintf('new_tSNE_plot_%s', expCondSep), sep = '/')
    if(!dir.exists(newResDir)) dir.create(newResDir)
    ## tsne plot
    tsneCluster        <- DimPlot(seuratObjFinal, reduction = "tsne", label = T, label.size = 6, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
    tsneClusterNolabel <- DimPlot(seuratObjFinal, reduction = "tsne", label = F, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
    if (length(levels(as.factor(seuratObjFinal$expCond)))>1) {
      ## relevel the 'expCond' for split.by= ordering
      tsneSplit          <- DimPlot(seuratObjFinal, reduction = "tsne", label = T, label.size = 4, repel = T, split.by = 'expCond') + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
    }
    ## -
    if (newAnnotation) {
      plotName1 = paste(newResDir, 'tsne_plot_noLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'tsne_plot_wLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('tsne_plot_wLabel_newAnnotation_expCondSep%s.pdf', expCondSep), sep = '/')
    } else {
      plotName1 = paste(newResDir, 'tsne_plot_noLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'tsne_plot_wLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('tsne_plot_wLabel_orgAnnotation_expCondSep%s.pdf', expCondSep), sep = '/')
    }
    ## -
    pdf(file = plotName1, width = 5.7, height = 6.7)
    print(tsneClusterNolabel + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))))
    dev.off()
    ## -
    pdf(file = plotName2, width = 5.7, height = 6.7)
    print(tsneCluster + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
    dev.off()
    ## -
    if ( length(levels(as.factor(seuratObjFinal$expCond))) > 1 ) {
      if ( length(levels(as.factor(seuratObjFinal$expCond))) == 2) {
        pdf(file = plotName3, width = 11, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) == 3) {
        pdf(file = plotName3, width = 13, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) == 4) {
        pdf(file = plotName3, width = 21, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) > 4) {
        pdf(file = plotName3, width = 5.5*length(levels(as.factor(seuratObjFinal$expCond))), height = 7)
      }
      # print(tsneSplit + theme1noLegend)
      print(tsneSplit + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
      dev.off()
    }
    ## ---
    ## 2. re-make UMAP plot
    newResDir          <- paste(resDir, sprintf('new_UMAP_plot_%s', expCondSep), sep = '/')
    if(!dir.exists(newResDir)) dir.create(newResDir)
    ## tsne plot
    umapCluster        <- DimPlot(seuratObjFinal, reduction = "umap", label = T, label.size = 6, repel = T) + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    umapClusterNolabel <- DimPlot(seuratObjFinal, reduction = "umap", label = F, repel = T) + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    if ( length(levels(as.factor(seuratObjFinal$expCond))) > 1 ) {
      ## relevel the 'expCond' for split.by= ordering
      tsneSplit          <- DimPlot(seuratObjFinal, reduction = "umap", label = T, label.size = 4, repel = T, split.by = 'expCond') + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    }
    ## -
    if (newAnnotation) {
      plotName1 = paste(newResDir, 'UMAP_plot_noLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'UMAP_plot_wLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('UMAP_plot_wLabel_newAnnotation_expCondSep%s.pdf', expCondSep), sep = '/')
    } else {
      plotName1 = paste(newResDir, 'UMAP_plot_noLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'UMAP_plot_wLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('UMAP_plot_wLabel_orgAnnotation_expCondSep%s.pdf', expCondSep), sep = '/')
    }
    ## -
    pdf(file = plotName1, width = 5.7, height = 8)
    print(umapClusterNolabel + theme1wLegend + guides(colour = guide_legend(nrow=7, byrow=TRUE, override.aes = list(size=6))))
    dev.off()
    ## -
    pdf(file = plotName2, width = 5.7, height = 8)
    print(umapCluster + theme1wLegend + guides(colour = guide_legend(nrow=7, byrow=TRUE, override.aes = list(size=6))) )
    dev.off()
    ## -
    if ( length(levels(as.factor(seuratObjFinal$expCond))) > 1 ) {
      if ( length(levels(as.factor(seuratObjFinal$expCond))) == 2 ) {
        pdf(file = plotName3, width = 11, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) == 3 ) {
        pdf(file = plotName3, width = 13, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) == 4 ) {
        pdf(file = plotName3, width = 21, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) > 4) {
        pdf(file = plotName3, width = 5.5*length(levels(as.factor(seuratObjFinal$expCond))), height = 7)
      }
      # print(tsneSplit + theme1noLegend)
      print(tsneSplit + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
      dev.off()
    }
    ## -
  }
  ## ---------------------------------------------------------------------------------------
  ## 2. summarize cell no in each idetified clusters, if clleNo summary will change automately based on above whether to update on 'expCondSep'
  if (clusterCellsNoSummary) {
    DefaultAssay(seuratObjFinal)   <- "RNA"
    ## ---
    print(sprintf('Start step1: summarizing on identified cell no. in each cluster'))
    ## 1. output cell no. for each identified cluster, and exp conditions within each cluster
    clusterCellNo                  <- as.data.frame(table(Idents(seuratObjFinal)))
    seuratObjFinal$clusterExpCond  <- paste(Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_')
    clusterCellExpNo               <- as.data.frame(table(seuratObjFinal@meta.data$clusterExpCond))
    clusterCellExpNo$cluster       <- sapply(strsplit(as.character(clusterCellExpNo$Var1), split = '_'), '[[', 1)
    # clusterCellExpNo$exp           <- sapply(strsplit(as.character(clusterCellExpNo$Var1), split = '-'), '[[', 2)
    clusterCellExpNo$exp           <- sapply(strsplit(as.character(clusterCellExpNo$Var1), split = '_'), tail, 1)
    library(reshape2)
    library(xlsx)
    clusterCellExpNoWide           <- dcast(data = clusterCellExpNo, cluster ~ exp, value.var = 'Freq')
    clusterCellExpNoWide[is.na(clusterCellExpNoWide)] <- 0
    clusterCellExpNoWidePer        <- clusterCellExpNoWide
    clusterCellExpNoWideColSum     <- colSums(clusterCellExpNoWidePer %>% dplyr::select(-cluster))
    for (i in 2:dim(clusterCellExpNoWidePer)[2]) {
      clusterCellExpNoWidePer[,i]      <- round(clusterCellExpNoWidePer[,i]*100/clusterCellExpNoWideColSum[i-1], digits = 2)
    }
    clusterCellNoComb1             <- merge(clusterCellNo, clusterCellExpNoWide, by.x = 'Var1', by.y = 'cluster' )
    clusterCellNoComb              <- merge(clusterCellNoComb1, clusterCellExpNoWidePer, by.x = 'Var1', by.y = 'cluster' )
    colnames(clusterCellNoComb)    <- c('clusters', 'cellNo', paste('cellNo', colnames(clusterCellNoComb)[-c(1,2)], sep = '_'))
    colnames(clusterCellNoComb)    <- gsub(pattern = '.x', replacement = '', x = colnames(clusterCellNoComb))
    colnames(clusterCellNoComb)    <- gsub(pattern = '.y', replacement = '_Per', x = colnames(clusterCellNoComb))
    if (newAnnotation) {
      cellnoFname                  <- paste(resDir, sprintf('cellNo_summary_newAnnotation_expCond_%s.txt', expCondSep), sep = '/')
    } else {
      cellnoFname                  <- paste(resDir, sprintf('cellNo_summary_orgClusterAnnotation_expCond_%s.txt', expCondSep), sep = '/')
    }
    write.table(x = clusterCellNoComb, file = cellnoFname, quote = F, row.names = F, col.names = T, sep = '\t')
    print(sprintf('END step1: summarizing on identified cell no. in each cluster'))
    ## -
  }
}

## ---------------------------------------------------------------------------------------
