#' getGoiDotplot() Function
#' @details
#' This function is used to make feature genes dotplot of a set of marker/feature genes provided with 'goiFname' in an excel file.
#' @param resDir specify an exiting full path of directory, where results will be saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved with above rdsDir option.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param goiFname path to file, where a list of marker/features genes are provided in column 'Gene', if column 'Cell Type' is also provided, option 'geneTypeOrder' can be used to adjust orders.
#' @param geneTypeOrder if column 'Cell Type' is provided in the 'goiFname' file, this option can be used to adjust the orders of marker gene's cell types, if not provided, marker gene's cell types will be sorted alphabetically.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param expCondCheckFname suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param expCondReorderLevels a character string of the corresponding experimental condition factor levels' orders presented on the y-axis of the dot-plot from bottom to top, if not defined, sorted numerically or alphabetically.
#' @param expCond subset the specified experimental condition corresponding to 'expCondCheck' levels to make the corresponding dot plot.
#' @param cellcluster specify cell clusters to be displayed on the dot plot.
#' @param dotPlotFnamePrefix prefix of the dot plot file name, if not defined, by default = 'goiDotplots'.
#' @param dotPlotMinExpCutoff minimum expression value threshold presented in the dot plot, if not defined, by default = 0.3.
#' @param dotPlotMaxExpCutoff maximum expression value threshold presented in the dot plot, if not defined, by default = 2.5.
#' @param dotPlotWidth dot plot width, if not defined, will be decided automatically based on the number of marker genes presented in 'goiFname'
#' @param dotPlotHeight dot plot height, if not defined, will be decided automatically based on the number of experimental condition or sample's cell clusters.
#' @param legendPer specify the legend proportion in the dot plot, if not specified, by default = 0.1
#' @param genetypebarPer specify the proportion top gene types color bar, by default 0.01.
#' @param fontsize.x specify the font size of gene names on x-axis, by default 24.
#' @param fontsize.y specify the font size of samples attributes on y-axis, by default 18.
#' @param fontangle.x change the x-axis text label position, by default vertically aligned with x-axis with 90 degree.
#' @param fontangle.y change the y-axis text label position, by default horizontally displayed with y-axis with 90 degree.
#' @param fontsize.legend1 specify the font size of expression color legend bar.
#' @param fontsize.legend2 specify the font size of input gene list's types legend bar.
#' @param gridOn logical option, whether to have grid on the plot background, by default it is on.
#' @param geneTypeLegendOn logical option, whether to include feature genes type legend bar in the plot, by default 'on'.
#' @param debug logical option, turn on to have extra analysis information printing out for debug, by defalut 'off'.
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
#' @importFrom ggplot2 element_line
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
#' @examples getGoiDotplot(rds, goiFname, expCondCheck='sample/idents/expCond*')
#' @export
#' @return the dotplots of provided GOI(gene of interest) saved in '' inside the provided 'resDir'
## ---------------------------------------------------------------------------------------
# library(ggplot2)
# library(Seurat)
# library(grDevices)
# library(tools)
# library(xlsx)
## personal notes: turn on scale in Seurat::DotPlot(), normalized count will be scale() on each group.by conditions. cluster.idents=F do not scale across gene features to avoid high expressed gene strench the visualization.
getGoiDotplot <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, goiFname, geneTypeOrder=NULL,
                          expCondCheck='sample', expCondCheckFname = NULL, expCondReorderLevels=NULL,
                          expCond = NULL, cellcluster = NULL, dotPlotFnamePrefix='goiDotplots',
                          dotPlotMinExpCutoff=0.3, dotPlotMaxExpCutoff = 2.5, dotPlotWidth=NULL,
                          scale.min = NA, scale.max = NA,
                          dotPlotHeight=NULL, legendPer=NULL, genetypebarPer=NULL, fontsize.x = 24, fontsize.y = 18, fontangle.x = 90, fontangle.y = 0, fontsize.legend1 = 20, fontsize.legend2 = NULL, gridOn = as.logical(T), geneTypeLegendOn = as.logical(T), debug = F ){
  ## ---
  newAnnotation           <- as.logical(newAnnotation)
  expCondName             <- expCond ## for subsetting on expCond
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
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir                <- paste(resDir, 'results_wNewAnnotation', sep = '/')
  } else {
    resDir                <- paste(resDir, 'results_wOrgClusterAnnotation', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    if(debug){
      print("Before adding new annotation")
      print(table(Idents(seuratObjFinal)))
    }
    source(newAnnotationRscriptName)
    if (debug) {
      print("After adding new annotation")
      print(table(Idents(seuratObjFinal)))
    }
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
  ## create corresponding 'plotResDir' for feature-plot to save
  plotResDir            <- paste(resDir, sprintf('dotplot_selected_markers_%s', expCondCheckFname), sep = '/')
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
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
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(as.character(expCondCheck), colnames(seuratObjFinal@meta.data))]
    }
  }
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
  print(sprintf('GOI dot plots will be saved in %s', plotResDir))
  ##--------------------------------------------------------------------------------------##
  if (debug) print(table(seuratObjFinal@meta.data$expCond))
  expCondLevels = levels(factor(seuratObjFinal@meta.data$expCond))
  if (!is.null(expCondName)) {
    if (any(!expCondName %in% expCondLevels ) ) stop('Please provide the corresponding experimental condition.')
    if (length(expCondName)==1) {
      print(sprintf('Subsetting a specific experimental condition level: %s', expCondName))
      seuratObjFinal                          <- subset(seuratObjFinal, expCond == expCondName )
    } else {
      print(sprintf('Subsetting %s specific experimental condition levels: %s', length(expCondName), paste(expCondName, collapse = ', ')))
      for (e in 1:length(expCondName)) {
        seuratObjFinalPrep                    <- subset(seuratObjFinal, subset = expCond == expCondName[e] )
        # print(seuratObjFinalPrep)
        if (e == 1) {
          seuratObjFinal2                     <- seuratObjFinalPrep
        } else {
          seuratObjFinal2                     <- merge(seuratObjFinal2, seuratObjFinalPrep)
        }
      }
      seuratObjFinal <- seuratObjFinal2
    }
    if(debug) print('-=-=-=-=-=-=-=-=-=-=-')
    if(debug) print('subsetting on expCond')
    if(debug) print("expCond table is shown below")
    if(debug) print(table(seuratObjFinal@meta.data$expCond))
    if(debug) print('******')
    if(debug) print(head(seuratObjFinal))
    if(debug) print("Idents table is shown below")
    if(debug) print(table(Seurat::Idents(seuratObjFinal)))
    if(debug) print('-=-=-=-=-=-=-=-=-=-=-')
  }

  ## ------
  clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% clusterLevels ) ) stop('Please provide the corresponding cell clusters ')
    print(sprintf('Subsetting %s specific cell clusters: %s', length(cellcluster), paste(cellcluster, collapse = ',')))
    seuratObjFinal                          <- subset(seuratObjFinal, idents = cellcluster )
    if (debug) print("subset cell clusters")
    if(debug) print("subsetted idents table is shown below")
    if(debug) print(table(Seurat::Idents(seuratObjFinal)))
    if(debug) print('-=-=-=-=-=-=-=-=-=-=-')
  }
  ##--------------------------------------------------------------------------------------##
  ## using column 'gene' as marker genes, if column 'geneType' exist, will be used to categorize the genes
  if (file_ext(basename(goiFname)) == 'xlsx') {
    markerGenesPrep         <- read.xlsx(file = as.character(goiFname), sheetIndex = 1, header = T)
  } else if (file_ext(basename(goiFname)) == 'txt') {
    markerGenesPrep         <- read.delim(file = as.character(goiFname), header = T, sep = '\t')
  }
  colnames(markerGenesPrep) <- tolower(colnames(markerGenesPrep))
  if (debug) print(head(markerGenesPrep))
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
  if ('genetype' %in% colnames(markerGenesPrep) ){
    markerGenesCat        <- as.factor(markerGenesPrep$genetype)
  } else {
    markerGenesCat        <- NULL
  }
  ## -
  if (!is.null(markerGenesCat)){
    markerGenesDf         <- data.frame('gene' = markerGenes, 'genetype'=markerGenesCat)
    if (!is.null(geneTypeOrder)) {
      markerGenesDf$genetype <- factor(markerGenesDf$genetype, levels = geneTypeOrder)
      markerGenesDf          <- markerGenesDf[order(markerGenesDf$genetype),]
    }
  }
  ## Note: if column 'geneType' does not exist, markerGenesCat=NULL, NO 'markerGenesDf' exist
  ##--------------------------------------------------------------------------------------##
  dotplotFname           <- file.path(plotResDir, sprintf('%s_markerGenes_dotplot_%s.pdf', dotPlotFnamePrefix, expCondCheckFname))
  ## -
  print(sprintf("A total of %s marker genes will be used for downstream dotplot at '%s'. ", length(markerGenes), basename(dotplotFname) ))
  ## ---
  ## 2. make dotplot with provided gene markers; dot plot of all selected marker genes presented on x-axis
  SeuratObject::DefaultAssay(seuratObjFinal)   <- "RNA" ## suggested by Seurat tutorial at 'https://satijalab.org/seurat/articles/integration_introduction.html' using RNA slot for feature plot and dot plot and downstream analysis
  ## ---------
  ## adding expCond to the idents of identified clusters
  ## '_' is used here to combine idents() with expCond' factor levels.
  if (all( names(Seurat::Idents(seuratObjFinal)) == names(seuratObjFinal$expCond) )) {
    ## re-define Seurat::Idents(seuratObjFinal) with level name: cluster + expCond if expCondCheck != 'comb'
    if (expCondCheck != 'idents') {
      # Seurat::Idents(seuratObjFinal)        <- factor( paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'),
      #                                                  levels = paste(rep(levels(Seurat::Idents(seuratObjFinal)), each = length(levels(factor(seuratObjFinal$expCond))) ), levels(factor(seuratObjFinal$expCond, levels = expCondReorderLevels)), sep = '_') )
      if (is.null(expCondReorderLevels)) {
        expCondReorderLevels                <- levels(factor(paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_')))
      }
      ## ---
      if ( !all(levels(factor(paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'))) %in% expCondReorderLevels ) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort y-axis dot plot experimental conditions.")
      ## ---
      Seurat::Idents(seuratObjFinal)        <- factor(paste(Seurat::Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'), levels =expCondReorderLevels )
    } else {
      if (is.null(expCondReorderLevels)) {
        expCondReorderLevels                <- levels(factor(Seurat::Idents(seuratObjFinal)))
      }
      ## ---
      if ( !all(levels(factor(Seurat::Idents(seuratObjFinal)))  %in% expCondReorderLevels) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort y-axis dot plot experimental conditions.")
      ## ---
      Seurat::Idents(seuratObjFinal)        <- factor(Seurat::Idents(seuratObjFinal), levels = expCondReorderLevels)

    }
    ## ---

       ## ---
  } else{
    print("Error: cell name does not match.")
  }
  if(debug) {
    print(sprintf('Updated idents information with experimental condition are as below:'))
    print(table(Seurat::Idents(seuratObjFinal) ))
    print('"-=-=-=-=-=')
  }
  ##--------------------------------------------------------------------------------------##
  ## Currently Seurat::DotPlot() x-axis is based on Seurat::DotPlot()$data$features.plot levels after dropping off unused levels by droplevels()
  if (is.null(markerGenesCat)) {
    g1                <- Seurat::DotPlot(seuratObjFinal, features = markerGenes, cols = c('#D3D3D3', '#CC0000'), col.min = dotPlotMinExpCutoff, col.max = dotPlotMaxExpCutoff, scale = T, scale.by = 'size', dot.min = 0.01, scale.min = scale.min, scale.max = scale.max ) + Seurat::RotatedAxis()
    # print(Seurat::DotPlot(seuratObjFinal, features = markerGenes, cols = c('#D3D3D3', '#CC0000'), col.min = 0.3, scale = T, scale.by = 'size', dot.min = 0.01, idents = levels(Seurat::Idents(seuratObjFinal))[-c(7,14)] ) + Seurat::RotatedAxis())
    # print(Seurat::DotPlot(seuratObjFinal, features = markerGenes, cols = c('#D3D3D3', '#CC0000'), col.min = dotPlotMinExpCutoff, scale = T, scale.by = 'size', dot.min = 0.01 ) + Seurat::RotatedAxis())
    # print(Seurat::DotPlot(seuratObjFinal, features = markerGenes, col.min = 0.3 ) + Seurat::RotatedAxis())
    # print(Seurat::DotPlot(seuratObjFinal, features = markerGenes ) + Seurat::RotatedAxis())
  } else {
    g1                <- Seurat::DotPlot(seuratObjFinal, features = markerGenesDf$gene, cols = c('#D3D3D3', '#CC0000'), col.min = dotPlotMinExpCutoff, col.max = dotPlotMaxExpCutoff, scale = T, scale.by = 'size', dot.min = 0.01, scale.min = scale.min, scale.max = scale.max ) + Seurat::RotatedAxis()
    if (debug) {
      print(g1)
      print('0909090909')
    }

  }
  # print(g1)
  # ggsave(filename = 'TEST1.pdf', plot = g1, width = 40, height = 15)
  ##--------------------------------------------------------------------------------------##
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
  if (fontangle.y < 90) {
    hjustVal = 1 ##right align
    vjustVal.y = 0.5
  } else if (fontangle.y > 90 & fontangle.y < 180) {
    hjustVal = 0 ##left align
    vjustVal.y = 1
  } else {
    hjustVal = 0 ##left align
    vjustVal.y = 0.5
  }
  ##--------------------------------------------------------------------------------------##
  if (is.null(markerGenesCat)) {
    pdf(file = dotplotFname, width = dotPlotWidth, height = dotPlotHeight )
    g1     <- g1 + theme(legend.text = element_text(color = "black", size = fontsize.legend1), axis.text.x = element_text(color = "black", size = fontsize.x, angle = fontangle.x, vjust = 0.5), axis.text.y = element_text(color = "black", size = fontsize.y, angle = fontangle.y, hjust = hjustVal, vjust = vjustVal.y), axis.title = element_blank() )
    if (gridOn) {
      g1 <- g1 + theme(panel.grid.major = element_line(colour = "grey85"))
    }
    print(g1)
    dev.off()
  } else {
    ## 'markerGenesDf' exist with 2 columns 'gene' & 'geneType'
    markerGenesDf2      <- markerGenesDf %>% dplyr::filter(gene %in% levels(droplevels(g1$data$features.plot)) )
    markerGenesDf2$gene <- factor(markerGenesDf2$gene, levels = levels(droplevels(g1$data$features.plot)) )
    ## below DiscretePalette() in Seurat with 4 color scheme options, color scheme display seen in https://kwstat.github.io/pals/
    colschemeg2         <- Seurat::DiscretePalette(n = length(levels(markerGenesDf2$genetype)), palette = 'glasbey')
    g2     <- ggplot2::ggplot(markerGenesDf2) + ggplot2::geom_bar(mapping = aes(x = gene, y = 1, fill = genetype), stat = "identity", width = 1)
    g2     <- g2 + ggplot2::scale_fill_manual(values=colschemeg2)
    # ggsave(filename = 'TEST11.pdf', plot = g2, width = 40, height = 2)
    # g2     <- g2 + theme(panel.spacing.x = grid::unit(1, "mm"))
    g2     <- g2 + ggplot2::theme_void() + theme(panel.spacing.x = grid::unit(1, "mm")) ##+ facet_grid(.~colorBar, scales = "free_x") ##testing to make sure gene name aligned
    if (is.null(fontsize.legend2)) fontsize.legend2 = fontsize.legend1
    g2     <- g2 + theme(legend.text = element_text(color = "black", size = fontsize.legend2), legend.title = element_blank() )
    if (gridOn) {
      g1 <- g1 + theme(panel.grid.major = element_line(colour = "grey80"))
    }
    g1     <- g1 + theme(legend.text = element_text(color = "black", size = fontsize.legend1), axis.text.x = element_text(color = "black", size = fontsize.x, angle = fontangle.x, vjust = 0.5), axis.text.y = element_text(color = "black", size = fontsize.y, angle = fontangle.y, hjust = hjustVal, vjust = vjustVal.y), axis.title = element_blank() )
    legend              <- cowplot::plot_grid(cowplot::get_legend(g2), cowplot::get_legend(g1), ncol = 1, align = 'h', axis = 'l')
    legend1             <- cowplot::plot_grid(cowplot::get_legend(g1), ncol = 1, align = 'h', axis = 'l')
    g1woLegend          <- g1 + theme(legend.position = "none")
    g2woLegend          <- g2 + theme(legend.position = "none")
    if (is.null(genetypebarPer)) genetypebarPer <- 0.01
    plot                <- cowplot::plot_grid(g2woLegend, g1woLegend, align = "v", ncol = 1, axis = 'lr', rel_heights = c(genetypebarPer*dotPlotHeight, (1-genetypebarPer)*dotPlotHeight))
    # if (length(levels(factor(Seurat::Idents(seuratObjFinal) ))) < 151 ){
    #   plot              <- cowplot::plot_grid(g2woLegend, g1woLegend, align = "v", ncol = 1, axis = 'lr', rel_heights = c(0.015*dotPlotHeight, 0.95*dotPlotHeight))
    # } else {
    #   plot              <- cowplot::plot_grid(g2woLegend, g1woLegend, align = "v", ncol = 1, axis = 'lr', rel_heights = c(0.03*dotPlotHeight, 0.95*dotPlotHeight))
    # }
    if(is.null(legendPer)) legendPer <- 0.1
    if (geneTypeLegendOn) {
      plotWlegend         <- cowplot::plot_grid(plot, legend, nrow = 1, align = 'h', axis = 'none', rel_widths = c((1-legendPer)*dotPlotWidth, legendPer*dotPlotWidth))
      ggplot2::ggsave(filename = dotplotFname, plot = plotWlegend, width = dotPlotWidth, height = dotPlotHeight, limitsize = FALSE)
    } else {
      plotWlegend2        <- cowplot::plot_grid(plot, legend1, nrow = 1, align = 'h', axis = 'none', rel_widths = c((1-legendPer)*dotPlotWidth, legendPer*dotPlotWidth))
      ggplot2::ggsave(filename = dotplotFname, plot = plotWlegend2, width = dotPlotWidth, height = dotPlotHeight, limitsize = FALSE)
    }
  }
  ##--------------------------------------------------------------------------------------##
}
## ---------------------------------------------------------------------------------------
