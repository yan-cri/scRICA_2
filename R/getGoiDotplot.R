#' getGoiDotplot() Function
#' @details
#' This function is used to make dotplot of marker/features genes in provided 'goiFname'.
#'
#' @param resDir full path of integration results analysis returned in getClusterMarkers()
#' @param newAnnotation logical value, whether to provide manual annotation
#' @param newAnnotationRscriptName if newAnnotation == T, this script is used to redefine the old clusters
#' @param expCondSepName character string, user defined name either to be 'org' or any character string
#' @param expCondName2change if above 'expCondSepName' is defined not as 'org', provide the name to be changed
#' @param goiFname full path of a file name, where a list of marker/features genes provided
#' @param dotPlotFnamePrefix dotPlot file name prefix
#' @param dotPlotMinExpCutoff dotPlot miniumn expression threshold
#' @param dotPlotWidth dotPlot width
#' @param dotPlotHeight dotPlot height
#'
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 labs
#' @importFrom Seurat NoLegend
#' @importFrom Seurat DefaultAssay
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom utils write.table
#' @importFrom tools file_ext
#' @importFrom utils read.delim
#' @importFrom xlsx read.xlsx
#'
#' @keywords GoiDotplot
#' @examples getGoiDotplot()
#' @export
#' @return the dotplots of provided GOI(gene of interest) saved in '' inside the provided 'resDir'
## ---------------------------------------------------------------------------------------
getGoiDotplot <- function(resDir, newAnnotation, newAnnotationRscriptName, expCondSepName, expCondName2change, goiFname, dotPlotFnamePrefix, dotPlotMinExpCutoff, dotPlotWidth, dotPlotHeight ){
  if (missing(expCondName2change)) expCondName2change <- NA
  if (expCondSepName == 'org' | expCondSepName == 'comb' & !is.na(expCondName2change)) print("'expCondName2change' will not be applied")
  ## ------
  resDirName              <- strsplit(x = resDir, split = '/')[[1]][length(strsplit(x = resDir, split = '/')[[1]])]
  # print(sprintf('TEST TETST  resDirName is %s', resDirName))
  rdsFname                <- paste(resDir, sprintf("RDS_Dir/%s.rds", resDirName), sep = '/' )
  # print(sprintf('85858 rdsFname is %s', rdsFname))
  if (!file.exists(rdsFname)) stop("Please execute getClusterMarker() to conduct integration analysis before running getClusterSummaryReplot().")
  seuratObjFinal          <<- readRDS(file = as.character(rdsFname))
  print('Done for RDS readin')
  ## ------
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir <- paste(resDir, 'results_wNewAnnotation', sep = '/')
  } else {
    resDir <- paste(resDir, 'results_wOrgClusterAnnotation', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  ## -------------------------------------------------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  }
  ## -------------------------------------------------------------------------------------
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding 'plotResDir' for dotplot to save
  print(sprintf("expCondSep is '%s'", expCondSepName))
  plotResDir             <- paste(resDir, sprintf('dotplot_selected_markers_ExpCond_%s', expCondSepName), sep = '/')
  if (expCondSepName == 'org') {
    seuratObjFinal                     <- seuratObjFinal
  } else if (expCondSepName == 'comb') {
    seuratObjFinal@meta.data$expCond   <- Idents(seuratObjFinal)
  } else {
    seuratObjFinal@meta.data$expCond   <- gsub(pattern = as.character(expCondName2change), replacement = '', x = seuratObjFinal@meta.data$expCond)
  }
  # print('97979799')
  # print(table(seuratObjFinal@meta.data$expCond))
  # print(sprintf('4444 plotResDir is %s', plotResDir))
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
  print(sprintf('GOI dot plots will be saved in %s', plotResDir))
  ## ------
  ## 1. define/input selected markers/features
  ## below 'topUpDEmarkers' were obtained based on file 'paste(resDir, 'clusterDeMarkers_adjSig_up_wNewAnnotation.xlsx', sep = '/')'
  ## below 'topDownDEmarkers' were obtained based on file 'paste(resDir, 'clusterDeMarkers_adjSig_down_wNewAnnotation.xlsx', sep = '/')'
  print(table(Idents(seuratObjFinal)))
  if (file_ext(basename(goiFname)) == 'xlsx') {
    markerGenesPrep       <- read.xlsx(file = as.character(goiFname), sheetIndex = 1, header = T)
  } else if (file_ext(basename(goiFname)) == 'txt') {
    markerGenesPrep       <- read.delim(file = as.character(goiFname), header = T, sep = '\t')
  }
  ## ---
  if ("Gene" %in% colnames(markerGenesPrep)) {
    markerGenes           <- as.character(unique(markerGenesPrep$Gene))
  } else {
    markerGenes           <- as.character(unique(markerGenesPrep[,1]))
  }
  ## ---
  # dotplotFname           <- paste(plotResDir, 'selected_markerGenesV2_dotplot.pdf', sep = '/')
  dotplotFname           <- file.path(plotResDir, sprintf('%s_markerGenes_dotplot_expCondSep_%s.pdf', dotPlotFnamePrefix, expCondSepName))
  ## -
  print(sprintf("A total of %s marker genes will be used for downstream dotplot at '%s'. ", length(markerGenes), basename(dotplotFname) ))
  ## ---
  ## 2. make dotplot with provided gene markers; dot plot of all selected marker genes presented on x-axis
  DefaultAssay(seuratObjFinal)   <- "RNA"
  ## adding expCond to the idents of identified clusters
  if (all( names(Idents(seuratObjFinal)) == names(seuratObjFinal$expCond) )) {
    ## level name: cluster in front of expCond
    # Idents(seuratObjFinal) <- factor(paste(seuratObjFinal$seurat_clusters, seuratObjFinal$expCond, sep = '_'), levels = unlist(lapply(levels(seuratObjFinal$seurat_clusters), function(x) paste(x, levels(factor(seuratObjFinal$expCond)), sep = '_')) ) )
    ## level name: cluster after expCond
    if (expCondSepName != 'comb') {
      Idents(seuratObjFinal) <- factor(paste(Idents(seuratObjFinal), seuratObjFinal$expCond, sep = '_'), levels = unlist(lapply(levels(factor(seuratObjFinal$expCond)), function(x) paste(levels(Idents(seuratObjFinal)), x, sep = '_')) ) )
    }
  }
  print(sprintf('Updated idents information with experimental condition are as below:'))
  print(table(Idents(seuratObjFinal) ))
  print('"-=-=-=-=-=')
  ## ---
  pdf(file = dotplotFname, width = dotPlotWidth, height = dotPlotHeight )
  # print(DotPlot(seuratObjFinal, features = markerGenes, cols = c('#D3D3D3', '#CC0000'), col.min = 0.3, scale = T, scale.by = 'size', dot.min = 0.01, idents = levels(Idents(seuratObjFinal))[-c(7,14)] ) + RotatedAxis())
  ## for stem and proliferting cell to use
  print(DotPlot(seuratObjFinal, features = markerGenes, cols = c('#D3D3D3', '#CC0000'), col.min = dotPlotMinExpCutoff, scale = T, scale.by = 'size', dot.min = 0.01 ) + RotatedAxis())
  # print(DotPlot(seuratObjFinal, features = markerGenes, col.min = 0.3 ) + RotatedAxis())
  # print(DotPlot(seuratObjFinal, features = markerGenes ) + RotatedAxis())
  dev.off()
}
## ---------------------------------------------------------------------------------------
