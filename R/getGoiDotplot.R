#' getGoiDotplot() Function
#' @details
#' This function is used to make dotplot of marker/features genes in provided 'goiFname'.
#'
#' @param resDir full path of integration results analysis are saved, where RDS file is saved inside the 'RDS_Dir'. This path is also returned by getClusterMarkers() execution.
#' @param rdsFname User also can provide the full path of RDS file instead of 'resDir' where RDS file is saved in. If this option is used, please also provide 'resDir' to specify where the analysis results will be saved.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param goiFname path to file, where a list of marker/features genes are provided in column 'Gene', if column 'Cell Type' is also provided, option 'geneCellTypeOrder' can be used to adjust orders.
#' @param geneCellTypeOrder if column 'Cell Type' is provided in the 'goiFname' file, this option can be used to adjust the orders of marker gene's cell types, if not provided, marker gene's cell types will be sorted alphabetically.
#' @param expCondCheck 3 options: 'sample', 'expCond1', or 'expCond2' to specify which experimental conditions to be explored with this function.
#' @param expCondSepName suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param expCondName2change a character string to indicate part of characters specified here can be removed from sample name defined in the metadata table, if additional samples combination needs to be explored which has not been specified in the column of 'expCond1' or 'expCond2'.
#' @param expCondReorderLevels a character string of the corresponding experimental condition factor levels' orders presented on the y-axis of the dot-plot from bottom to top, if not defined, sorted numerically or alphabetically.
#' @param dotPlotFnamePrefix prefix of the dot plot file name, if not defined, by default = 'goiDotplots'.
#' @param dotPlotMinExpCutoff minimum expression value threshold presented in the dot plot, if not defined, by default = 0.3.
#' @param dotPlotWidth dot plot width, if not defined, will be decided automatically based on the number of marker genes presented in 'goiFname'
#' @param dotPlotHeight dot plot height, if not defined, will be decided automatically based on the number of experimental condition or sample's cell clusters.
#' @param legendPer specify the legend proportion in the dot plot, if not specified, by default = 0.1
#'
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 ggsave
#' @importFrom Seurat Idents
#' @importFrom Seurat NoLegend
#' @importFrom Seurat RotatedAxis
#' @importFrom SeuratObject DefaultAssay
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom utils write.table
#' @importFrom tools file_ext
#' @importFrom utils read.delim
#' @importFrom xlsx read.xlsx
#' @importFrom cowplot plot_grid
#' @importFrom cowplot get_legend
#' @importFrom grid unit
#'
#' @keywords GoiDotplot
#' @examples getGoiDotplot()
#' @export
#' @return the dotplots of provided GOI(gene of interest) saved in '' inside the provided 'resDir'
## ---------------------------------------------------------------------------------------
# library(ggplot2)
# library(Seurat)
# library(grDevices)
# library(tools)
# library(xlsx)
getGoiDotplot <- function(resDir=NULL, rdsFname=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, goiFname, geneCellTypeOrder=NULL, expCondCheck='sample', expCondSepName = NULL, expCondName2change=NULL, expCondReorderLevels=NULL, dotPlotFnamePrefix='goiDotplots', dotPlotMinExpCutoff=0.3, dotPlotWidth=NULL, dotPlotHeight=NULL, legendPer=NULL ){
  ## ---
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
  print('Done for RDS read in')
  ## ------
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir                <- paste(resDir, 'results_wNewAnnotation', sep = '/')
  } else {
    resDir                <- paste(resDir, 'results_wOrgClusterAnnotation', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
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
  ## -------------------------------------------------------------------------------------
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding 'plotResDir' for feature-plot to save
  plotResDir            <- paste(resDir, sprintf('dotplot_selected_markers_expCond_%s', expCondSepName), sep = '/')
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
  ## ------
  if (expCondCheck == 'sample') {
    seuratObjFinal                     <- seuratObjFinal
  } else if (expCondCheck == 'comb') {
    seuratObjFinal@meta.data$expCond   <- Seurat::Idents(seuratObjFinal)
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
  # print('97979799')
  # print(table(seuratObjFinal@meta.data$expCond))
  # print(sprintf('4444 plotResDir is %s', plotResDir))
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
  print(sprintf('GOI dot plots will be saved in %s', plotResDir))
  ## ------
  print(table(Seurat::Idents(seuratObjFinal)))
  ## ------
  ## using column 'gene' as marker genes, if column 'celltype' exist, will be used to categorize the genes
  if (file_ext(basename(goiFname)) == 'xlsx') {
    markerGenesPrep         <- read.xlsx(file = as.character(goiFname), sheetIndex = 1, header = T)
  } else if (file_ext(basename(goiFname)) == 'txt') {
    markerGenesPrep         <- read.delim(file = as.character(goiFname), header = T, sep = '\t')
  }
  colnames(markerGenesPrep) <- tolower(colnames(markerGenesPrep))
  # print(head(markerGenesPrep))
  print("----------------")
  if (sum(duplicated(markerGenesPrep$gene))>0) print(sprintf("%s genes are duplicated genes, they are: %s", sum(duplicated(markerGenesPrep$gene)), paste(markerGenesPrep$gene[duplicated(markerGenesPrep$gene)], collapse = ', ' ) ))
  if (dim(markerGenesPrep)[2] == 1) {
    markerGenesPrep       <- data.frame('gene'= markerGenesPrep[!duplicated(markerGenesPrep$gene),])
  } else {
    markerGenesPrep       <- markerGenesPrep[!duplicated(markerGenesPrep$gene),]
  }
  print(sprintf("A total of %s genes will be searched for GOI dot-plot", length(unique(markerGenesPrep$gene))))
  print("----------------")
  ## marker gene for dot plot is based on either 'gene' column or 1st column
  if ("gene" %in% colnames(markerGenesPrep)) {
    markerGenes           <- as.character(unique(markerGenesPrep$gene))
  } else {
    markerGenes           <- as.character(unique(markerGenesPrep[,1]))
  }
  ## -
  if ('celltype' %in% colnames(markerGenesPrep) ){
    markerGenesCat        <- as.factor(markerGenesPrep$celltype)
  } else {
    markerGenesCat        <- NULL
  }
  ## -
  if (!is.null(markerGenesCat)){
    markerGenesDf         <- data.frame('gene' = markerGenes, 'celltype'=markerGenesCat)
    if (!is.null(geneCellTypeOrder)) {
      markerGenesDf$celltype <- factor(markerGenesDf$celltype, levels = geneCellTypeOrder)
      markerGenesDf          <- markerGenesDf[order(markerGenesDf$celltype),]
    }
  }
  ## Note: if column 'celltype' does not exist, markerGenesCat=NULL, NO 'markerGenesDf' exist
  ## -------------------------------------------------------------------------------------
  # dotplotFname           <- paste(plotResDir, 'selected_markerGenesV2_dotplot.pdf', sep = '/')
  dotplotFname           <- file.path(plotResDir, sprintf('%s_markerGenes_dotplot_expCond_%s.pdf', dotPlotFnamePrefix, expCondSepName))
  ## -
  print(sprintf("A total of %s marker genes will be used for downstream dotplot at '%s'. ", length(markerGenes), basename(dotplotFname) ))
  ## ---
  ## 2. make dotplot with provided gene markers; dot plot of all selected marker genes presented on x-axis
  SeuratObject::DefaultAssay(seuratObjFinal)   <- "RNA" ## suggested by Seurat tutorial at 'https://satijalab.org/seurat/articles/integration_introduction.html' using RNA slot for feature plot and dot plot and downstream analysis
  # SeuratObject::DefaultAssay(seuratObjFinal)   <- "integrated"
  ## adding expCond to the idents of identified clusters
  if (all( names(Seurat::Idents(seuratObjFinal)) == names(seuratObjFinal$expCond) )) {
    ## re-define Seurat::Idents(seuratObjFinal) with level name: cluster + expCond if expCondCheck != 'comb'
    ## ---
    if (is.null(expCondReorderLevels)) {
      expCondReorderLevels          <- levels(factor(seuratObjFinal@meta.data$expCond))
    }
    if ( !all(expCondReorderLevels %in% levels(factor(seuratObjFinal$expCond))) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort y-axis dot plot")
    if (expCondCheck != 'comb') {
      Seurat::Idents(seuratObjFinal)        <- factor( paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'),
                                                       levels = paste(rep(levels(Seurat::Idents(seuratObjFinal)), each = length(levels(factor(seuratObjFinal$expCond))) ), levels(factor(seuratObjFinal$expCond, levels = expCondReorderLevels)), sep = '_') )
    } else {
      Seurat::Idents(seuratObjFinal)        <- factor(Seurat::Idents(seuratObjFinal), levels = expCondReorderLevels)
    }

    # if (!is.null(expCondReorderLevels)) {
    #   if (cellclusterNameSort) {
    # if ( !all(expCondReorderLevels %in% levels(factor(seuratObjFinal$expCond))) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort y-axis dot plot")
    # Seurat::Idents(seuratObjFinal) <- factor(paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'), levels = paste(rep(levels(Seurat::Idents(seuratObjFinal)), each = length(levels(factor(seuratObjFinal$expCond))) ), levels(factor(seuratObjFinal$expCond, levels = expCondReorderLevels)), sep = '_') )
    #   } else {
    # Seurat::Idents(seuratObjFinal) <- factor(paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'), levels = unlist(lapply(levels(factor(seuratObjFinal$expCond, levels = expCondReorderLevels)), function(x) paste(levels(Seurat::Idents(seuratObjFinal)), x, sep = '_')) ) )
    #   }
    # } else {
    #   if (cellclusterNameSort) {
    #     Seurat::Idents(seuratObjFinal) <- factor(paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'), levels = paste(rep(levels(Seurat::Idents(seuratObjFinal)), each = length(levels(factor(seuratObjFinal$expCond))) ), levels(factor(seuratObjFinal$expCond)), sep = '_') )
    #   } else {
    #     Seurat::Idents(seuratObjFinal) <- factor(paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'), levels = unlist(lapply(levels(factor(seuratObjFinal$expCond)), function(x) paste(levels(Seurat::Idents(seuratObjFinal)), x, sep = '_')) ) )
    #   }
    # }
    ## ---
  }
  print(sprintf('Updated idents information with experimental condition are as below:'))
  print(table(Seurat::Idents(seuratObjFinal) ))
  print('"-=-=-=-=-=')
  ## -----------------------------------------------------------------------------
  ## Currently Seurat::DotPlot() x-axis is based on Seurat::DotPlot()$data$features.plot levels after dropping off unused levels by droplevels()
  if (is.null(markerGenesCat)) {
    g1                <- Seurat::DotPlot(seuratObjFinal, features = markerGenes, cols = c('#D3D3D3', '#CC0000'), col.min = dotPlotMinExpCutoff, scale = T, scale.by = 'size', dot.min = 0.01 ) + Seurat::RotatedAxis()
    # print(Seurat::DotPlot(seuratObjFinal, features = markerGenes, cols = c('#D3D3D3', '#CC0000'), col.min = 0.3, scale = T, scale.by = 'size', dot.min = 0.01, idents = levels(Seurat::Idents(seuratObjFinal))[-c(7,14)] ) + Seurat::RotatedAxis())
    # print(Seurat::DotPlot(seuratObjFinal, features = markerGenes, cols = c('#D3D3D3', '#CC0000'), col.min = dotPlotMinExpCutoff, scale = T, scale.by = 'size', dot.min = 0.01 ) + Seurat::RotatedAxis())
    # print(Seurat::DotPlot(seuratObjFinal, features = markerGenes, col.min = 0.3 ) + Seurat::RotatedAxis())
    # print(Seurat::DotPlot(seuratObjFinal, features = markerGenes ) + Seurat::RotatedAxis())

  } else {
    g1                <- Seurat::DotPlot(seuratObjFinal, features = markerGenesDf$gene, cols = c('#D3D3D3', '#CC0000'), col.min = dotPlotMinExpCutoff, scale = T, scale.by = 'size', dot.min = 0.01 ) + Seurat::RotatedAxis()
  }
  # ggsave(filename = 'TEST1.pdf', plot = g1, width = 40, height = 15)
  ## ------
  print(sprintf('A total of %s genes displayed in GOI dot plot.', length(unique(g1$data$features.plot))))
  ## ---
  if (is.null(dotPlotWidth)) {
    dotPlotWidth  = round(length(unique(g1$data$features.plot))*0.5)
  }
  if (is.null(dotPlotHeight)) {
    if (length(levels(factor(Seurat::Idents(seuratObjFinal) ))) < 15 ) {
      dotPlotHeight = round(length( levels(factor(Seurat::Idents(seuratObjFinal))))*0.6 )
    } else if (length(levels(factor(Seurat::Idents(seuratObjFinal) ))) > 14 &  length(levels(factor(Seurat::Idents(seuratObjFinal) ))) < 25 ) {
      dotPlotHeight = round(length( levels(factor(Seurat::Idents(seuratObjFinal))))*0.4 )
    } else {
      dotPlotHeight = round(length( levels(factor(Seurat::Idents(seuratObjFinal))))*0.25 )
    }
  }
  ## ---
  if (is.null(markerGenesCat)) {
    pdf(file = dotplotFname, width = dotPlotWidth, height = dotPlotHeight )
    print(g1)
    dev.off()
  } else {
    ## 'markerGenesDf' exist with 2 columns 'gene' & 'celltype'
    markerGenesDf2      <- markerGenesDf %>% dplyr::filter(gene %in% levels(droplevels(g1$data$features.plot)) )
    markerGenesDf2$gene <- factor(markerGenesDf2$gene, levels = levels(droplevels(g1$data$features.plot)) )
    ## below DiscretePalette() in Seurat with 4 color scheme options, color scheme display seen in https://kwstat.github.io/pals/
    colschemeg2         <- Seurat::DiscretePalette(n = length(levels(markerGenesDf2$celltype)), palette = 'glasbey')
    g2     <- ggplot2::ggplot(markerGenesDf2) + ggplot2::geom_bar(mapping = aes(x = gene, y = 1, fill = celltype), stat = "identity", width = 1)
    g2     <- g2 + ggplot2::scale_fill_manual(values=colschemeg2)
    # ggsave(filename = 'TEST11.pdf', plot = g2, width = 40, height = 2)
    # g2     <- g2 + theme(panel.spacing.x = grid::unit(1, "mm"))
    g2     <- g2 + ggplot2::theme_void() + theme(panel.spacing.x = grid::unit(1, "mm")) ##+ facet_grid(.~colorBar, scales = "free_x") ##testing to make sure gene name aligned
    g2     <- g2 + theme(legend.text = element_text(color = "black", size = 20), legend.title = element_blank() )
    g1     <- g1 + theme(legend.text = element_text(color = "black", size = 20), axis.text.x = element_text(color = "black", size = 24), axis.text.y = element_text(color = "black", size = 18), axis.title = element_blank() )
    legend              <- cowplot::plot_grid(cowplot::get_legend(g2), cowplot::get_legend(g1), ncol = 1, align = 'h', axis = 'l')
    g1woLegend          <- g1 + theme(legend.position = "none")
    g2woLegend          <- g2 + theme(legend.position = "none")

    if (length(levels(factor(Seurat::Idents(seuratObjFinal) ))) < 151 ){
      plot              <- cowplot::plot_grid(g2woLegend, g1woLegend, align = "v", ncol = 1, axis = 'lr', rel_heights = c(0.015*dotPlotHeight, 0.95*dotPlotHeight))
    } else {
      plot              <- cowplot::plot_grid(g2woLegend, g1woLegend, align = "v", ncol = 1, axis = 'lr', rel_heights = c(0.03*dotPlotHeight, 0.95*dotPlotHeight))
    }
    if(is.null(legendPer)) legendPer <- 0.1
    plotWlegend         <- cowplot::plot_grid(plot, legend, nrow = 1, align = 'h', axis = 'none', rel_widths = c((1-legendPer)*dotPlotWidth, legendPer*dotPlotWidth))
    ggplot2::ggsave(filename = dotplotFname, plot = plotWlegend, width = dotPlotWidth, height = dotPlotHeight, limitsize = FALSE)
  }
}
## ---------------------------------------------------------------------------------------
