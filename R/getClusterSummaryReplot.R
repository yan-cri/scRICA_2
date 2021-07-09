#' getClusterSummaryReplot() Function
#' @details
#' This function is used to annotate/assign and merge identified clusters,
#'
#' @param resDir full path of integration results analysis returned in getClusterMarkers()
#' @param newAnnotation logical value, whether to provide manual annotation
#' @param newAnnotationRscriptName if newAnnotation == T, this script is used to redefine the old clusters
#' @param expCondSepName character string, user defined name either to be 'org' or any character string
#' @param expCondName2change if above 'expCondSepName' is defined not as 'org', provide the name to be changed
#' @param expCondNameReorder whether to change the experimental conditions order in summarized feature plots and percentage bar plots
#'
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom Seurat NoLegend
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @importFrom Seurat RenameIdents
#' @importFrom Seurat DimPlot
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom utils write.table
#' @importFrom ggplot2 element_text
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#'
#' @keywords clusterSummary
#' @examples clusterSummary()
#' @export
#' @return
## ---------------------------------------------------------------------------------------
getClusterSummaryReplot <- function(resDir, newAnnotation, newAnnotationRscriptName, expCondSepName, expCondName2change, clusterLevelReorder = T, expCondNameReorder = NULL) {
  rdsFname                <- paste(resDir, "RDS_Dir/analysis_results_integration_results.rds", sep = '/' )
  # print(sprintf('85858 rdsFname is %s', rdsFname))
  if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
  seuratObjFinal          <<- readRDS(file = as.character(rdsFname))
  print('Done for RDS readin')
  ## ------
  clusteringPlotRemake    <- as.logical(T)
  clusterCellsNoSummary   <- as.logical(T)
  ## ------
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir <- paste(resDir, 'results_wNewAnnotation', sep = '/')
  } else {
    resDir <- paste(resDir, 'results_wOrgClusterAnnotation', sep = '/')
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
  # setup cutome theme for plotting
  theme1noLegend          <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                                   # legend.key.size = unit(0.7, "cm"),
                                   axis.title = element_text(size = 20),
                                   axis.text = element_text(size = 25),
                                   legend.position="bottom",
                                   legend.text = element_text(size = 14) ) + NoLegend()
  theme1wLegend           <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                                   # legend.key.size = unit(0.7, "cm"),
                                   axis.title = element_text(size = 15),
                                   axis.text = element_text(size = 20),
                                   legend.position="bottom",
                                   legend.text = element_text(size = 15) )
  theme1wLegendRight      <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                                   # legend.key.size = unit(0.7, "cm"),
                                   axis.title = element_text(size = 20),
                                   axis.text = element_text(size = 25),
                                   legend.position="right",
                                   legend.text = element_text(size = 14) )
  ## -------------------------------------------------------------------------------------
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding updated 'resDir' for new tSNE/UMAP plots to save
  resDir                  <- paste(resDir, sprintf('expCond_%s', expCondSepName), sep = '/')
  ## ---------------------------------------------------------------------------------------
  if (expCondSepName == 'org') {
    seuratObjFinal        <- seuratObjFinal
  } else {
    seuratObjFinal@meta.data$expCond <- gsub(pattern = as.character(expCondName2change), replacement = '', x = seuratObjFinal@meta.data$expCond)
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  ## ---
  if (!is.null(expCondNameReorder)) seuratObjFinal$expCond <- factor(seuratObjFinal$expCond, levels = expCondNameReorder )
  ## -------------------------------------------------------------------------------------
  ## 1. re-make tSNE plot
  if (clusteringPlotRemake) {
    if (clusterLevelReorder) {
      Idents(seuratObjFinal) <- factor(Idents(seuratObjFinal), levels = fpClusterOrder )
    }
    ## color options
    ggplotColours <- function(n = 6, h = c(0, 360) + 15){
      if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
      hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
    }
    ## -
    # selectedCol  <- DiscretePalette(n = length(levels(Idents(seuratObjFinal))), palette = 'alphabet')
    selectedCol        <- ggplotColours(n=length(levels(Idents(seuratObjFinal))))
    ## ---
    print(sprintf('Start step1: remake tSNE/UMAP plots'))
    newResDir             <- paste(resDir, sprintf('new_tSNE_plot_%s', expCondSepName), sep = '/')
    if(!dir.exists(newResDir)) dir.create(newResDir)
    ## tsne plot
    tsneCluster           <- DimPlot(seuratObjFinal, reduction = "tsne", cols = selectedCol, label = T, label.size = 6, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
    tsneClusterNolabel    <- DimPlot(seuratObjFinal, reduction = "tsne", cols = selectedCol, label = F, repel = T) + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
    if (length(levels(as.factor(seuratObjFinal$expCond)))>1) {
      ## relevel the 'expCond' for split.by= ordering
      tsneSplit           <- DimPlot(seuratObjFinal, reduction = "tsne", cols = selectedCol, label = T, label.size = 4, repel = T, split.by = 'expCond') + labs(title = 'tSNE clustering', x = "tSNE 1", y = 'tSNE 2')
    }
    ## -
    if (newAnnotation) {
      plotName1 = paste(newResDir, 'tsne_plot_noLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'tsne_plot_wLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('tsne_plot_wLabel_newAnnotation_expCondSep%s.pdf', expCondSepName), sep = '/')
    } else {
      plotName1 = paste(newResDir, 'tsne_plot_noLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'tsne_plot_wLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('tsne_plot_wLabel_orgAnnotation_expCondSep%s.pdf', expCondSepName), sep = '/')
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
    newResDir          <- paste(resDir, sprintf('new_UMAP_plot_%s', expCondSepName), sep = '/')
    if(!dir.exists(newResDir)) dir.create(newResDir)
    ## umap plot
    umapCluster        <- DimPlot(seuratObjFinal, reduction = "umap", cols = selectedCol, label = T, label.size = 6, repel = T) + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    umapClusterNolabel <- DimPlot(seuratObjFinal, reduction = "umap", cols = selectedCol, label = F, repel = T) + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    if ( length(levels(as.factor(seuratObjFinal$expCond))) > 1 ) {
      ## relevel the 'expCond' for split.by= ordering
      tsneSplit          <- DimPlot(seuratObjFinal, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T, split.by = 'expCond') + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    }
    ## -
    if (newAnnotation) {
      plotName1 = paste(newResDir, 'UMAP_plot_noLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'UMAP_plot_wLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('UMAP_plot_wLabel_newAnnotation_expCondSep_%s.pdf', expCondSepName), sep = '/')
    } else {
      plotName1 = paste(newResDir, 'UMAP_plot_noLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'UMAP_plot_wLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('UMAP_plot_wLabel_orgAnnotation_expCondSep_%s.pdf', expCondSepName), sep = '/')
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
    print(sprintf('END step1: remake tSNE/UMAP plots'))
    ## -
  }
  ## -------------------------------------------------------------------------------------
  ## 2. summarize cell no in each idetified clusters, if clleNo summary will change automately based on above whether to update on 'expCondSepName'
  if (clusterCellsNoSummary) {
    Seurat::DefaultAssay(seuratObjFinal)   <- "RNA"
    ## ---
    print(sprintf('Start step2.1: summarizing on identified cell no. in each cluster'))
    ## 1. output cell no. for each identified cluster, and exp conditions within each cluster
    clusterCellNo                  <- as.data.frame(table(Idents(seuratObjFinal)))
    seuratObjFinal$clusterExpCond  <- paste(Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_')
    clusterCellExpNo               <- as.data.frame(table(seuratObjFinal@meta.data$clusterExpCond))
    clusterCellExpNo$cluster       <- sapply(strsplit(as.character(clusterCellExpNo$Var1), split = '_'), '[[', 1)
    # clusterCellExpNo$exp           <- sapply(strsplit(as.character(clusterCellExpNo$Var1), split = '-'), '[[', 2)
    clusterCellExpNo$exp           <- sapply(strsplit(as.character(clusterCellExpNo$Var1), split = '_'), tail, 1)
    # library(reshape2)
    # library(xlsx)
    clusterCellExpNoWide           <- reshape2::dcast(data = clusterCellExpNo, cluster ~ exp, value.var = 'Freq')
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
      cellnoFname                  <- paste(resDir, sprintf('cellNo_summary_newAnnotation_expCond_%s.txt', expCondSepName), sep = '/')
    } else {
      cellnoFname                  <- paste(resDir, sprintf('cellNo_summary_orgClusterAnnotation_expCond_%s.txt', expCondSepName), sep = '/')
    }
    write.table(x = clusterCellNoComb, file = cellnoFname, quote = F, row.names = F, col.names = T, sep = '\t')
    print(sprintf('END step2.1: summarizing on identified cell no. in each cluster'))
    ## -
    print(sprintf('Start step2.2: plotting identified cell no. percentage in each cluster'))
    # print(head(clusterCellExpNoWidePer))
    # print(dim(clusterCellExpNoWidePer))
    perData2plotLong               <- reshape2::melt(clusterCellExpNoWidePer, id.vars = c('cluster'))
    if (newAnnotation) {
      perData2plotLong$cluster     <- factor(perData2plotLong$cluster, levels = perClusterOrder )
    }
    ## ---------
    if (!is.null(expCondNameReorder)){
      perData2plotLong$variable    <- factor(perData2plotLong$variable, levels = rev(expCondNameReorder) )
    }
    ## ---------
    selectedCol2 <- selectedCol[match( perClusterOrder, fpClusterOrder)]
    # print(sprintf('dimention of clusterCellExpNoWidePer is %s, %s', dim(clusterCellExpNoWidePer)[1], dim(clusterCellExpNoWidePer)[2] ))
    if (dim(clusterCellExpNoWidePer)[2] > 2) {
      plotSizeHeight               <- round( (0.5*dim(clusterCellExpNoWidePer)[2]), digits = 0)
      plotSize <- c( round(dim(clusterCellExpNoWidePer)[1], digits = 0), round( (0.5*dim(clusterCellExpNoWidePer)[2]), digits = 0) )
    } else {
      plotSizeHeight               <- round( (0.8*dim(clusterCellExpNoWidePer)[2]), digits = 0)
      plotSize <- c( round(dim(clusterCellExpNoWidePer)[1], digits = 0), round( (0.8*dim(clusterCellExpNoWidePer)[2]), digits = 0) )
    }
    ## -
    if (dim(clusterCellExpNoWidePer)[1] < 11) {
      plotSizeWidth                <- round(dim(clusterCellExpNoWidePer)[1], digits = 0)
    } else if (dim(clusterCellExpNoWidePer)[1] > 10 & dim(clusterCellExpNoWidePer)[1] < 17 ) {
      plotSizeWidth                <- round(0.7*dim(clusterCellExpNoWidePer)[1], digits = 0)
    }  else if (dim(clusterCellExpNoWidePer)[1] > 16) {
      plotSizeWidth                <- round(0.5*dim(clusterCellExpNoWidePer)[1], digits = 0)
    }
    ## ---------
    plotSize <- c( plotSizeWidth, plotSizeHeight )
    pdf(file = gsub('.txt', '.pdf', cellnoFname), width = plotSizeWidth, height = plotSizeHeight )
    g1 <- ggplot2::ggplot(perData2plotLong, ggplot2::aes(x = value, y = factor(variable), fill = factor(cluster) )) + ggplot2::geom_bar(stat="identity")
    g1 <- g1 + ggplot2::scale_fill_manual(values=selectedCol2)
    g1 <- g1 + ggplot2::labs(title='', x = '', y = '')
    g1 <- g1 + ggplot2::labs(fill = '')
    g1 <- g1 + theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
                     axis.title.x = element_text(color="black", size=16, face="bold"),
                     axis.title.y = element_text(color="black", size=16, face="bold"))
    g1 <- g1 + theme(axis.text.x = element_text(face="plain", color="black", size=20, angle=0),
                     axis.text.y = element_text(face="plain", color="#000000", size=20, angle=0))
    g1 <- g1 + theme(legend.title = element_text(color = "black", size = 18),
                     legend.text = element_text(color = "black", size = 18) )
    if (dim(clusterCellExpNoWidePer)[1] > 10) {
      g1 <- g1 + guides(fill=guide_legend(ncol=2))
    }
    print(g1)
    dev.off()
    print(sprintf('END step2.2: plotting identified cell no. percentage in each cluster'))
    ## -
  }
}

## ---------------------------------------------------------------------------------------
