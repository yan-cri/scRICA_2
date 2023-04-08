#' getClusterSummaryReplot() Function
#' @details
#' This function is used to 1). summarize the number of cells in each originally identified/annotated cells clusters; 2). re-make clustering reduction UMAP/tSNE plots; 3). make the cells percentage plot based on experimental condition levels specified in the metadata table.
#'
#' @param resDir specify an exiting full path of directory, where results will be saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved with above rdsDir option.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param expCondCheckFname suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param fpRemake whether to re-make feature plots, default is True
#' @param perRemake whether to re-calculate cell number percentage and related percentage plots, default is True
#' @param colors by default it is null, colors will be automatically selected, otherwise, a vector string of colors can be specified with this option for colors used in the plots.
#' @param cellcluster by default it is null, all cell clusters will be plot, otherwise, cell clusters presented in the plots (UMAP/tSNE and percentage plots) can be specified with this option.
#' @param fpClusterOrder character string of cell clusters orders to specify the cell clusters orders corresponding to the color scheme, if not defined, sorted numerically, alphabetically, or the orders shown in the annotation R script.
#' @param perClusterOrder character string of cell clusters orders to specify the cell clusters orders presented in the percentage bar plot, if not defined, sorted numerically or alphabetically.
#' @param expCondNameOrder character string of the corresponding experimental condition factor levels' orders presented in the percentage bar plot, if not defined, sorted numerically or alphabetically.
#' @param stack whether to make stacked bar-plot
#' @param gap space of un-stacked bar-plot between different conditions, 0.1-1.
#' @param debug logical option, turn on to have extra analysis information printing out for debug, by defalut 'off'.
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
#' @examples clusterSummary(rds, expCondCheck='sample/idents/expCond*')
#' @export
## ---------------------------------------------------------------------------------------
getClusterSummaryReplot <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL,
                                    fpRemake = T, perRemake = T, colors = NULL, cellcluster = NULL,
                                    expCondCheck='sample', expCondCheckFname = NULL, fpClusterOrder = NULL, perClusterOrder = NULL,
                                    perPlotVertical = F, perPlotHeight = NULL, perPlotWidth = NULL, expCondNameOrder = NULL,
                                    stack = as.logical('T'), gap = 0.8, debug = as.logical(F)) {
  ## ---
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) print("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
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
  ## ------
  clusteringPlotRemake    <- as.logical(fpRemake)
  clusterCellsNoSummary   <- as.logical(perRemake)
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
    if (debug) {
      print("Before adding new annotation")
      print(table(Idents(seuratObjFinal)))
    }
    source(newAnnotationRscriptName)
    if(debug){
      print("After adding new annotation")
      print(table(Idents(seuratObjFinal)))
    }
  }
  ## -------------------------------------------------------------------------------------
  # setup custom theme for plotting
  theme1noLegend          <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                                   # legend.key.size = unit(0.7, "cm"),
                                   axis.title = element_text(size = 20),
                                   axis.text = element_text(size = 25),
                                   legend.position="bottom",
                                   legend.text = element_text(size = 14) ) + Seurat::NoLegend()
  theme1wLegend           <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                                   # legend.key.size = unit(0.7, "cm"),
                                   axis.title = element_text(size = 15),
                                   axis.text = element_text(size = 25),
                                   legend.position="bottom",
                                   legend.text = element_text(size = 15) )
  theme1wLegendRight      <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                                   # legend.key.size = unit(0.7, "cm"),
                                   axis.title = element_text(size = 20),
                                   axis.text = element_text(size = 25),
                                   legend.position="right",
                                   legend.text = element_text(size = 14) )
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond'
  if (expCondCheck == 'sample') {
    seuratObjFinal                     <- seuratObjFinal
  } else if (expCondCheck == 'idents') {
    seuratObjFinal@meta.data$expCond   <- Seurat::Idents(seuratObjFinal)
  } else {
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  seuratObjFinal$expCond <- gsub(pattern = '_', replacement = '-', x = seuratObjFinal$expCond)
  if(debug){
    print('6666666666666')
    print(table(Seurat::Idents(seuratObjFinal)))
    print('6666666666666')
  }

  ##--------------------------------------------------------------------------------------##
  clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% clusterLevels ) ) stop('Please provide the corresponding cell clusters ')
    print(sprintf('Subsetting %s specific cell clusters: %s', length(cellcluster), paste(cellcluster, collapse = ',')))
    seuratObjFinal                          <- subset(seuratObjFinal, idents = cellcluster )
  }
  if(debug){
    print('8787878787878')
    print(table(Idents(seuratObjFinal)))
    print('8787878787878')
  }
  ##--------------------------------------------------------------------------------------##
  if (expCondCheck == 'sample') {
    if (is.null(expCondCheckFname)) {
      expCondCheckFname        <- 'expCond_sample'
    } else {
      expCondCheckFname        <- expCondCheckFname
    }
  } else {
    if (is.null(expCondCheckFname)) {
      expCondCheckFname        <- as.character(expCondCheck)
    } else {
      expCondCheckFname        <- expCondCheckFname
    }
  }
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding updated 'resDir' for new tSNE/UMAP plots to save
  resDir                <- sprintf('%s/expCond_%s', resDir, expCondCheckFname)
  if (!dir.exists(resDir)) dir.create(resDir)
  ## ---------------------------------------------------------------------------------------
  if (!is.null(expCondNameOrder)) {
    if (all((expCondNameOrder %in% levels(factor(seuratObjFinal$expCond))))) {
      seuratObjFinal$expCond  <- factor(seuratObjFinal$expCond, levels = expCondNameOrder )
    } else {
      stop("Error: please provide corresponding option 'expCondNameOrder' consistent with updated experimental conditions using option 'expCondCheck'. ")
    }
  }
  ## -------------------------------------------------------------------------------------
  if (is.null(fpClusterOrder)) fpClusterOrder <- levels(factor(Seurat::Idents(seuratObjFinal)))
  Seurat::Idents(seuratObjFinal) <- factor(Seurat::Idents(seuratObjFinal), levels = fpClusterOrder )
  ## -
  if (is.null(colors)) {
    # selectedCol  <- DiscretePalette(n = length(levels(Seurat::Idents(seuratObjFinal))), palette = 'alphabet')
    selectedCol        <- ggplotColours(n=length(levels(Seurat::Idents(seuratObjFinal))))
  } else {
    selectedCol        <- colors
  }
  print(sprintf("used colors are '%s' for cells: '%s'.", paste(selectedCol, collapse = ','), paste(fpClusterOrder, collapse = ',') ))
  ## -------------------------------------------------------------------------------------
  ## 1. re-make tSNE plot
  if (clusteringPlotRemake) {
    print('***************************************************')
    ## ---
    print(sprintf('Start step1: remake tSNE/UMAP plots'))
    newResDir             <- paste(resDir, sprintf('new_tSNE_plot_%s', expCondCheckFname), sep = '/')
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
      plotName3 = paste(newResDir, sprintf('tsne_plot_wLabel_newAnnotation_%s.pdf', expCondCheckFname), sep = '/')
    } else {
      plotName1 = paste(newResDir, 'tsne_plot_noLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'tsne_plot_wLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('tsne_plot_wLabel_orgAnnotation_%s.pdf', expCondCheckFname), sep = '/')
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
    ## --------------------------------------------------------
    ## 2. re-make UMAP plot
    newResDir          <- paste(resDir, sprintf('new_UMAP_plot_%s', expCondCheckFname), sep = '/')
    if(!dir.exists(newResDir)) dir.create(newResDir)
    ## umap plot
    umapCluster        <- DimPlot(seuratObjFinal, reduction = "umap", cols = selectedCol, label = T, label.size = 6, repel = T) + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    umapClusterNolabel <- DimPlot(seuratObjFinal, reduction = "umap", cols = selectedCol, label = F, repel = T) + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    if ( length(levels(as.factor(seuratObjFinal$expCond))) > 1 ) {
      ## relevel the 'expCond' for split.by= ordering
      umapSplit          <- DimPlot(seuratObjFinal, reduction = "umap", cols = selectedCol, label = T, label.size = 4, repel = T, split.by = 'expCond') + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
      umapSplitNolabel   <- DimPlot(seuratObjFinal, reduction = "umap", cols = selectedCol, label = F, label.size = 4, repel = T, split.by = 'expCond') + labs(title = 'UMAP clustering', x = "UMAP 1", y = 'UMAP 2')
    }
    ## -
    if (newAnnotation) {
      plotName1 = paste(newResDir, 'UMAP_plot_noLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'UMAP_plot_wLabel_integrate_newAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('UMAP_plot_wLabel_newAnnotation_%s.pdf', expCondCheckFname), sep = '/')
      plotName4 = paste(newResDir, sprintf('UMAP_plot_noLabel_newAnnotation_%s.pdf', expCondCheckFname), sep = '/')
    } else {
      plotName1 = paste(newResDir, 'UMAP_plot_noLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName2 = paste(newResDir, 'UMAP_plot_wLabel_integrate_orgAnnotation.pdf', sep = '/')
      plotName3 = paste(newResDir, sprintf('UMAP_plot_wLabel_orgAnnotation_%s.pdf', expCondCheckFname), sep = '/')
      plotName4 = paste(newResDir, sprintf('UMAP_plot_noLabel_orgAnnotation_%s.pdf', expCondCheckFname), sep = '/')
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
      print(umapSplit + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
      dev.off()
    }
    ## -
    if ( length(levels(as.factor(seuratObjFinal$expCond))) > 1 ) {
      if ( length(levels(as.factor(seuratObjFinal$expCond))) == 2 ) {
        pdf(file = plotName4, width = 11, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) == 3 ) {
        pdf(file = plotName4, width = 13, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) == 4 ) {
        pdf(file = plotName4, width = 21, height = 7)
      } else if ( length(levels(as.factor(seuratObjFinal$expCond))) > 4) {
        pdf(file = plotName4, width = 5.5*length(levels(as.factor(seuratObjFinal$expCond))), height = 7)
      }
      print(umapSplitNolabel + theme1wLegend + guides(colour = guide_legend(nrow=4, byrow=TRUE, override.aes = list(size=6))) )
      dev.off()
    }
    print(sprintf('END step1: remake tSNE/UMAP plots'))
    ## -
  }
  ## -------------------------------------------------------------------------------------
  ## 2. summarize cell no in each identified clusters, if cellNo summary will change automatically based on above whether to update on 'expCondCheckFname'
  ## '_' is used here to combine idents() with expCond' factor levels.
  if (clusterCellsNoSummary) {
    Seurat::DefaultAssay(seuratObjFinal)   <- "RNA"
    ## ---
    print(sprintf('Start step2.1: summarizing on identified cell no. in each cluster'))
    ## 1. output cell no. for each identified cluster, and exp conditions within each cluster
    clusterCellNo                  <- as.data.frame(table(Seurat::Idents(seuratObjFinal)))
    seuratObjFinal$clusterExpCond  <- paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_')
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
      cellnoFname                  <- paste(resDir, sprintf('cellNo_summary_newAnnotation_%s.txt', expCondCheckFname), sep = '/')
    } else {
      cellnoFname                  <- paste(resDir, sprintf('cellNo_summary_orgClusterAnnotation_%s.txt', expCondCheckFname), sep = '/')
    }
    write.table(x = clusterCellNoComb, file = cellnoFname, quote = F, row.names = F, col.names = T, sep = '\t')
    print(sprintf('END step2.1: summarizing on identified cell no. in each cluster'))
    ## -
    print(sprintf('Start step2.2: plotting identified cell no. percentage in each cluster'))
    if(debug){
      print(head(clusterCellExpNoWidePer))
      print(dim(clusterCellExpNoWidePer))
    }
    perData2plotLong               <- reshape2::melt(clusterCellExpNoWidePer, id.vars = c('cluster'))
    ## ---------
    ## updating annotated cell clusters levels order
    if (is.null(perClusterOrder))  perClusterOrder <- levels(factor(perData2plotLong$cluster))
    ## ---------
    ## updating 'expCond' levels order
    if (!is.null(expCondNameOrder)){
      if (perPlotVertical) {
        perData2plotLong$variable    <- factor(perData2plotLong$variable, levels = expCondNameOrder )
      } else {
        perData2plotLong$variable    <- factor(perData2plotLong$variable, levels = rev(expCondNameOrder) )
      }
    }
    ## ---------
    if (perPlotVertical) {
      perData2plotLong$cluster     <- factor(perData2plotLong$cluster, levels = perClusterOrder )
      selectedCol2 <- selectedCol[match( perClusterOrder, fpClusterOrder)]
    } else {
      perData2plotLong$cluster     <- factor(perData2plotLong$cluster, levels = rev(perClusterOrder) )
      selectedCol2 <- selectedCol[match( rev(perClusterOrder), fpClusterOrder)]
    }

    if (is.null(perPlotHeight)) {
      # print(sprintf('dimention of clusterCellExpNoWidePer is %s, %s', dim(clusterCellExpNoWidePer)[1], dim(clusterCellExpNoWidePer)[2] ))
      if (dim(clusterCellExpNoWidePer)[2] > 2) {
        plotSizeHeight               <- round( (0.5*dim(clusterCellExpNoWidePer)[2]), digits = 0)
        # plotSize <- c( round(dim(clusterCellExpNoWidePer)[1], digits = 0), round( (0.5*dim(clusterCellExpNoWidePer)[2]), digits = 0) )
      } else {
        plotSizeHeight               <- round( (0.8*dim(clusterCellExpNoWidePer)[2]), digits = 0)
        # plotSize <- c( round(dim(clusterCellExpNoWidePer)[1], digits = 0), round( (0.8*dim(clusterCellExpNoWidePer)[2]), digits = 0) )
      }
      if (perPlotVertical) plotSizeWidth <- plotSizeHeight+3
    } else {
      plotSizeHeight <- perPlotHeight
    }

    ## -
    if (is.null(perPlotWidth)) {
      if (dim(clusterCellExpNoWidePer)[1] < 11) {
        plotSizeWidth                <- round(dim(clusterCellExpNoWidePer)[1], digits = 0)
      } else if (dim(clusterCellExpNoWidePer)[1] > 10 & dim(clusterCellExpNoWidePer)[1] < 17 ) {
        plotSizeWidth                <- round(0.7*dim(clusterCellExpNoWidePer)[1], digits = 0)
      }  else if (dim(clusterCellExpNoWidePer)[1] > 16) {
        plotSizeWidth                <- round(0.5*dim(clusterCellExpNoWidePer)[1], digits = 0)
      }
      if (perPlotVertical) plotSizeHeight <- plotSizeWidth-2
    } else {
      plotSizeWidth <- perPlotWidth
    }
    ## -
    ## ---------
    # print(head(perData2plotLong))
    # print('777777777')

    if (perPlotVertical) {
      pdf(file = gsub('.txt', '_vertical.pdf', cellnoFname), width = plotSizeWidth, height = plotSizeHeight  )

      if (stack) {
        g1 <- ggplot2::ggplot(perData2plotLong, ggplot2::aes(x = factor(variable), y = value, fill = factor(cluster) )) + ggplot2::geom_bar(stat="identity")
        g1 <- g1 + ggplot2::scale_fill_manual(values=selectedCol2)
        g1 <- g1 + ggplot2::labs(title='', x = '', y = '')
        g1 <- g1 + ggplot2::labs(fill = '')
        g1 <- g1 + theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
                         axis.title.x = element_text(color="black", size=16, face="bold"),
                         axis.title.y = element_text(color="black", size=16, face="bold"))
        g1 <- g1 + theme(axis.text.x = element_text(face="plain", color="black", size=20, angle=90, vjust = 0, hjust = 1),
                         axis.text.y = element_text(face="plain", color="#000000", size=20, angle=0))
        g1 <- g1 + theme(legend.title = element_text(color = "black", size = 18),
                         legend.text = element_text(color = "black", size = 18) )
        if (dim(clusterCellExpNoWidePer)[1] > 10) {
          g1 <- g1 + guides(fill=guide_legend(ncol=1))
        }
      } else {
        g1 <- ggplot2::ggplot(perData2plotLong, ggplot2::aes(x = value, y = factor(cluster, levels = rev(perClusterOrder)), fill = factor(variable) ))+ ggplot2::geom_bar(stat="identity", position="dodge", colour="black", width = gap) + ggplot2::coord_flip()
        # g1 <- g1 + ggplot2::scale_fill_manual(values=selectedCol2)
        g1 <- g1 + ggplot2::labs(title='', y = 'Relative cell %', x = '')
        g1 <- g1 + ggplot2::labs(fill = '')
        g1 <- g1 + ggplot2::theme_void() + theme(panel.grid.major = element_line(colour = "grey85"))
        g1 <- g1 + ggplot2::theme(plot.margin = unit(c(0.5,0.1,0.5,0.1), "cm"))
        g1 <- g1 + ggplot2::labs(fill = '')
        g1 <- g1 + theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
                         axis.title.x = element_text(color="black", size=16, face="bold"),
                         axis.title.y = element_text(color="black", size=16, face="bold"))
        g1 <- g1 + theme(axis.text.x = element_text(face="plain", color="black", size=14, angle=90, vjust = 0, hjust = 1),
                         axis.text.y = element_text(face="plain", color="#000000", size=20, angle=0))
        g1 <- g1 + theme(legend.title = element_text(color = "black", size = 18),
                         legend.text = element_text(color = "black", size = 18) )
        if (dim(clusterCellExpNoWidePer)[1] > 10) {
          g1 <- g1 + guides(fill=guide_legend(ncol=1))
        }
      }
      print(g1)
      dev.off()
    } else {
      pdf(file = gsub('.txt', '.pdf', cellnoFname), width = plotSizeWidth, height = plotSizeHeight )
      if (stack) {
        g1 <- ggplot2::ggplot(perData2plotLong, ggplot2::aes(x = value, y = factor(variable), fill = factor(cluster) )) + ggplot2::geom_bar(stat="identity")
        g1 <- g1 + ggplot2::scale_fill_manual(values=selectedCol2)
        g1 <- g1 + ggplot2::labs(title='', x = '', y = '')
        g1 <- g1 + ggplot2::labs(fill = '')
        g1 <- g1 + theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
                         axis.title.x = element_text(color="black", size=14, face="bold"),
                         axis.title.y = element_text(color="black", size=20, face="bold"))
        g1 <- g1 + theme(axis.text.x = element_text(face="plain", color="black", size=20, angle=0 ),
                         axis.text.y = element_text(face="plain", color="#000000", size=20, angle=0, hjust=0))
        g1 <- g1 + theme(legend.title = element_text(color = "black", size = 18),
                         legend.text = element_text(color = "black", size = 18) )
        if (dim(clusterCellExpNoWidePer)[1] > 10) {
          g1 <- g1 + guides(fill=guide_legend(ncol=2, reverse = T))
        }
      } else {
        g1 <- ggplot2::ggplot(perData2plotLong, ggplot2::aes(x = value, y = factor(cluster, levels = rev(perClusterOrder)), fill = factor(variable) ))
        g1 <- g1 + ggplot2::geom_bar(stat="identity", colour="black", position = 'dodge', width = gap )
        # g1 <- g1 + ggplot2::scale_fill_manual(values=selectedCol2)
        g1 <- g1 + ggplot2::labs(title='', x = 'Relative cell %', y = '')
        g1 <- g1 + ggplot2::labs(fill = '')
        g1 <- g1 + ggplot2::theme_void() + theme(panel.grid.major = element_line(colour = "grey85"))
        g1 <- g1 + ggplot2::theme(plot.margin = unit(c(0.1,0.5,0.1,0.5), "cm"))
        g1 <- g1 + theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
                         axis.title.x = element_text(color="black", size=20, face="bold"),
                         axis.title.y = element_text(color="black", size=14, face="bold"))
        g1 <- g1 + theme(axis.text.x = element_text(face="plain", color="black", size=20, angle=0 ),
                         axis.text.y = element_text(face="plain", color="#000000", size=14, angle=0, hjust=0))
        g1 <- g1 + theme(legend.title = element_text(color = "black", size = 20),
                         legend.text = element_text(color = "black", size = 24) )
        if (dim(clusterCellExpNoWidePer)[1] > 10) {
          g1 <- g1 + guides(fill=guide_legend(ncol=1, reverse = T))
        }
      }
      print(g1)
      dev.off()
    }
    print(sprintf('END step2.2: plotting identified cell no. percentage in each cluster'))
    ## -
  }
}

## color options
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
## ---------------------------------------------------------------------------------------
