#' getGoiHeatmap() Function
#'
#' @details
#' This function is used to make heatmap of marker/features genes in provided 'goiFname' or geneNames with 3 different visualization options on 1). identified cell clusters ('idents'); 2). specific samples attributes from an experimental condition levels ('expCond'); 3). combined levels of cell clusters and samples attributes ('expCond.idents).
#'
#' @param heatmap.view specify heatmap visualization options, including 3 options: 'idents', 'expCond', and 'expCond.idents' to make heatmap, by default 'idents'.
#' @param resDir specify an exiting full path of directory, where results will be saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved with above rdsDir option.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param goiFname path to file, where a list of marker/features genes are provided in column 'Gene', if column 'Cell Type' is also provided, option 'geneTypeOrder' can be used to adjust orders.
#' @param geneNames provide gene names directly with this option, if do not want to use option 'goiFname'
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param expCondCheckFname suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param cellcluster specify cell clusters to be displayed on the dot plot
#' @param expCond specify the specific experimental conditions to be plotted.
#' @param expCondReorderLevels a character string of the corresponding experimental condition factor levels' orders shown on top color bar.
#' @param plotFnamePrefix prefix of the heatmap plot file name, if not defined, by default = 'goiHeatmap'.
#' @param plotWidth dot plot width, if not defined, will be decided automatically based on the number of marker genes presented in 'goiFname'
#' @param plotHeight dot plot height, if not defined, will be decided automatically based on the number of experimental condition or sample's cell clusters.
#' @param fontangle.top modify the angle of the text on the top legend color bar.
#' @param fontsize.top modify the font size of top legend color bar.
#' @param fontsize.y modify the font size of y-axis (gene names).
#' @param fontangle.y modify the y-axis (gene names) angle position.
#' @param display 3 options: 'cell', 'sample.merge' or 'expCond.merge'.
#' @param draw.lines whether to have an empty cell space on heatmap to separate different groups.
#' @param debug if on, extra debug printing meassages displayed, by default it is off.
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
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 scale_color_discrete
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 annotation_raster
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 coord_cartesian
#' @importFrom Seurat SingleRasterMap
#' @importFrom Seurat Idents
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
#' @importFrom scales hue_pal
#' @importFrom ggplot2 geom_text
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat %||%
#'
#' @keywords Heatmap
#' @examples getGoiHeatmap(rds, geneNames/goiFname, expCondCheck='sample/idents/expCond*')
#' @export
#' @return the Heatmap of provided GOI(gene of interest) provided by 'geneNames' directly or presented in 'goiFname'.
#'
## ---------------------------------------------------------------------------------------
getGoiHeatmap <- function(heatmap.view = 'idents', resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, slot = "scale.data",
                          goiFname = NULL, geneNames = NULL,
                          expCondCheck='sample', expCondCheckFname = NULL,
                          cellcluster = NULL , expCond = NULL, expCondReorderLevels = NULL,
                          plotFnamePrefix='goiHeatmap', plotWidth=25, plotHeight = 20,
                          res.output = as.logical(F), outputFname = NULL,
                          minVal = -2.5, maxVal = NULL, fontsize.legend = 20,
                          features.process = as.logical('F'), features.process.lists = NULL, features.process.fns = '+',
                          fontangle.top = 0, fontsize.top = 5, fontsize.y = 20, fontangle.y = 0, hjust.top = 0, vjust.top = 0, barHeight.top = 0.02,
                          display = 'cell', draw.lines = as.logical(T), debug = as.logical(F)) {
  ## ----
  maxVal <- maxVal %||% ifelse(test = slot == "scale.data", yes = 2.5, no = 6)
  sel.expConds <- expCond ## to avoid 'expCond' otpion with metadata 'expCond' rename this paratmer into sel.expConds
  if (is.null(geneNames) & is.null(goiFname)) {
    stop("Please provide genes to make heatmap with either 'goiFname' or 'geneNames' otpion.")
  } else if (!is.null(geneNames) & !is.null(goiFname)) {
    stop("Please provide genes to make heatmap with only 'goiFname' or 'geneNames' otpion.")
  } else if (!is.null(geneNames) & is.null(goiFname)) {
    features = geneNames
    if (length(unique(features)) < length(geneNames)) {
      print(sprintf("%s genes provided in the option 'geneNames' for heatmap, out of them, %s (%s) are duplicated genes.",
                    length(geneNames), (length(geneNames)-length(unique(features))), paste(unlist(geneNames[duplicated(geneNames)]), collapse = ', ') ))
    } else {
      print(sprintf("%s genes provided in the option 'geneNames' for heatmap.", length(geneNames)))
    }
  } else if (is.null(geneNames) & !is.null(goiFname)) {
    if (file_ext(basename(goiFname)) == 'xlsx') {
      markerGenesPrep         <- read.xlsx(file = as.character(goiFname), sheetIndex = 1, header = T)
    } else if (file_ext(basename(goiFname)) == 'txt') {
      markerGenesPrep         <- read.delim(file = as.character(goiFname), header = T, sep = '\t')
    }
    colnames(markerGenesPrep) <- tolower(colnames(markerGenesPrep))
    if ('gene' %in% colnames(markerGenesPrep)) {
      features = markerGenesPrep$gene
    } else {
      features = markerGenesPrep[,1]
    }
  }
  ## ----
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
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
  rdsDir                <- sprintf('%s/RDS_Dir', resDir)
  if (!dir.exists(rdsDir)) dir.create(rdsDir)
  ## -------------------------------------------------------------------------------------
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir                <- paste(resDir, 'results_wNewAnnotation', sep = '/')
  } else {
    resDir                <- paste(resDir, 'results_wOrgClusterAnnotation', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
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
  ## create corresponding 'plotResDir' for heatmap to save
  plotResDir            <- paste(resDir, sprintf('heatmap_%s', expCondCheckFname), sep = '/')
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
  print(sprintf('Heatmap plots with %s experimental condition will be saved at %s', expCondCheckFname, plotResDir))
  ## -------------------------------------------------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
    if (debug) {
      print('-=-=-=-')
      print('updated Idents are as below:')
      print(table(Idents(seuratObjFinal)))
      print('-=-=-=-')
    }
  }
  ##--------------------------------------------------------------------------------------##
  ## Assumption: original 'expCond' column represents samples information
  seuratObjFinal@meta.data$orgSample <- seuratObjFinal$expCond
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
  ##--------------------------------------------------------------------------------------##
  ## if provided, subset on 'cellcluster'
  clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% clusterLevels ) ) stop('Please provide the corresponding cell clusters in identfied idents.')
    print(sprintf('Subsetting %s specific cell clusters: %s', length(cellcluster), paste(cellcluster, collapse = ',')))
    seuratObjFinal        <- subset(seuratObjFinal, idents = cellcluster )
  }
  ## in provided, subset on 'expCond'
  expCondLevels           <- levels(factor(seuratObjFinal@meta.data$expCond))
  if (!is.null(sel.expConds)) {
    if (any(!sel.expConds %in% expCondLevels ) ) stop("Please provide the corresponding experimental condition levesl specified in 'expCondCheck' option.")
    print(sprintf('Subsetting %s specific expCond: %s', length(sel.expConds), paste(sel.expConds, collapse = ',')))

    if (length(sel.expConds)==1) {
      seuratObjFinal          <- subset(seuratObjFinal, expCond == sel.expConds)
    } else {
      for (i in 1:length(sel.expConds)) {
        seuratObjFinalPrep    <- subset(seuratObjFinal, subset = expCond == as.character(sel.expConds[i]))
        if (i ==1) {
          seuratObjFinalPrep2 = seuratObjFinalPrep
        } else {
          seuratObjFinalPrep2 <- merge(seuratObjFinalPrep2, seuratObjFinalPrep)
        }
      }
      seuratObjFinal          <- seuratObjFinalPrep2
    }
  }
  ##--------------------------------------------------------------------------------------##
  Seurat::DefaultAssay(seuratObjFinal) <- "RNA"

  plotTheme <- ggplot2::theme(axis.text.y = element_text(color = "black", size = fontsize.y, angle = fontangle.y),
                              legend.text = element_text(color = "black", size = fontsize.legend),
                              legend.title = element_blank() )
  if (slot == "scale.data") {
    seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T, features = features) ## by default, both do.scale/center are on, will scale the expression level for each gene feature by dividing the centered feature expression levels by their standard deviations
  } else {
    print("Heatmap is based on normalized data only.")
  }
  print('Start: Making feature genes heatmap plot')
  if (heatmap.view == 'idents') {
    ## ---
    print('Start to make heatmap seperated on Idents.')
    seuratObjFinal$expCond = Idents(seuratObjFinal)
    if (!is.null(expCondReorderLevels)) {
      if ( !all(expCondReorderLevels %in% levels(Idents(seuratObjFinal))) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort heatmap top legend bars.")
      seuratObjFinal$expCond <- factor(seuratObjFinal$expCond, levels = expCondReorderLevels)
      if(debug) print((sprintf("expCondReorderLevels: %s", paste(expCondReorderLevels, collapse = ','))))
    }
    seuratObjFinal@meta.data$orgSample <- paste(seuratObjFinal$expCond, seuratObjFinal@meta.data$orgSample, sep = '_')
    features.heatmap1 <- DoHeatmap2(object = seuratObjFinal, features = features, slot = slot,
                                    res.output = res.output, outputFname = outputFname,
                                    group.by = 'expCond', angle = fontangle.top, size = fontsize.top, hjust = hjust.top, vjust =  vjust.top, group.bar.height = barHeight.top,
                                    features.process = as.logical(features.process), features.process.lists = features.process.lists, features.process.fns = features.process.fns,
                                    disp.min = minVal, disp.max = maxVal, display = display, draw.lines = draw.lines, expCondReorderLevels = expCondReorderLevels, debug = debug)
    pdf(file = file.path(plotResDir, sprintf('%s_heatmap.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    print(features.heatmap1+ plotTheme + guides(colour = 'none'))
    dev.off()
    ## ---
  } else  if (heatmap.view == 'expCond') {
    ## ---
    print('Start to make heatmap seperated on expCond')
    if (!is.null(expCondReorderLevels)) {
      # print(levels(factor(seuratObjFinal$expCond)))
      if ( !all(expCondReorderLevels %in% levels(factor(seuratObjFinal$expCond))) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort heatmap top legend bars.")
      seuratObjFinal@meta.data$expCond <- factor(seuratObjFinal@meta.data$expCond, levels = expCondReorderLevels)
      if(debug) print((sprintf("expCondReorderLevels: %s", paste(expCondReorderLevels, collapse = ','))))
    }
    if (debug) {
      print(table(seuratObjFinal@meta.data$expCond))
    }
    features.heatmap1 <- DoHeatmap2(object = seuratObjFinal, features = features, assay = 'RNA', slot = slot,
                                    res.output = res.output, outputFname = outputFname,
                                    group.by = 'expCond', angle = fontangle.top, size = fontsize.top, hjust = hjust.top, vjust =  vjust.top, group.bar.height = barHeight.top,
                                    features.process = as.logical(features.process), features.process.lists = features.process.lists, features.process.fns = features.process.fns,
                                    disp.min = minVal, disp.max = maxVal, display = display, draw.lines = draw.lines, expCondReorderLevels = expCondReorderLevels, debug = debug)
    pdf(file = file.path(plotResDir, sprintf('%s_heatmap.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    print(features.heatmap1+ plotTheme + guides(colour = 'none'))
    # print(features.heatmap1+ plotTheme )
    dev.off()
    ## ---
  } else if (heatmap.view == 'expCond.idents') {
    ## ---
    print('Start to make heatmap seperated on expCond and idents combinations.')
    seuratObjFinal <- AddMetaData(object = seuratObjFinal, metadata = paste(seuratObjFinal$expCond, Idents(seuratObjFinal), sep = ' '), col.name = 'expIdentComb')
    seuratObjFinal$expCond <- seuratObjFinal$expIdentComb
    if (debug) {
      print(table(seuratObjFinal@meta.data$expCond))
    }
    ## -
    if (!is.null(expCondReorderLevels)) {
      if ( !all(expCondReorderLevels %in% levels(factor(seuratObjFinal$expCond))) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort heatmap top legend bars.")
      seuratObjFinal@meta.data$expCond <- factor(seuratObjFinal@meta.data$expCond, levels = expCondReorderLevels)
      if(debug) print((sprintf("expCondReorderLevels: %s", paste(expCondReorderLevels, collapse = ','))))
    }
    ## -
    seuratObjFinal@meta.data$orgSample <- paste(seuratObjFinal$expCond, seuratObjFinal@meta.data$orgSample, sep = '_')
    ## -
    features.heatmap1 <- DoHeatmap2(object = seuratObjFinal, features = features, assay = 'RNA', slot = slot,
                                    res.output = res.output, outputFname = outputFname,
                                    group.by = 'expCond', angle = fontangle.top, size = fontsize.top, hjust = hjust.top, vjust =  vjust.top, group.bar.height = barHeight.top,
                                    features.process = as.logical(features.process), features.process.lists = features.process.lists, features.process.fns = features.process.fns,
                                    disp.min = minVal, disp.max = maxVal, display = display, draw.lines = draw.lines, expCondReorderLevels = expCondReorderLevels, debug = debug)
    pdf(file = file.path(plotResDir, sprintf('%s_heatmap.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    print(features.heatmap1+ plotTheme + guides(colour = 'none'))
    dev.off()
    ## ---
  }
  print('END: Making selected genes heatmap plot')
  print('********************')
  ##--------------------------------------------------------------------------------------##
}

##----------------------------------------------------------------------------------------- ##
## DoHeatmap2 adding 'vjust' to adjust top color legend bar text position from Seurat::DoHeatmap
DoHeatmap2 <- function (object, features = NULL, cells = NULL, group.by = "expCond",
          group.bar = TRUE, group.colors = NULL, disp.min = -2.5, disp.max = NULL,
          slot = NULL, assay = 'RNA', label = TRUE, size = 5.5,
          hjust = 0, vjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE,
          features.process = features.process, features.process.lists = features.process.lists, features.process.fns = features.process.fns,
          lines.width = NULL, group.bar.height = 0.02, combine = TRUE, res.output = as.logical(F), outputFname = NULL,
          display = 'cell', expCondReorderLevels = expCondReorderLevels, debug = as.logical(F))
{
  if(debug) {
    print(sprintf("DoHeatmap2 use expCondReorderLevels: %s", paste(expCondReorderLevels, collapse = ',')))
  }
  if (features.process & is.null(features.process.lists)) stop("ERROR: 'features.process' is on, please provide corresponding 'features.process.lists'.")
  if(is.null(slot)){slot = "scale.data"}
  print(sprintf("Note: heatmap is using '%s'.", slot))
  if(display == 'expCond.merge') draw.lines = as.logical(F)
  if(debug) {
    if (display == 'expCond.merge') print("display == 'expCond.merge', 'draw.lines' is not working.")
  }
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% Seurat::DefaultAssay(object = object)
  Seurat::DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", yes = 2.5, no = 6)
  possible.features <- rownames(x = SeuratObject::GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot,
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ",
            slot, " slot for the ", assay, " assay: ", paste(bad.features,
                                                             collapse = ", "))
  }
  if (display == 'sample.merge') {
    data0 <- as.data.frame(x = as.matrix(x = t(x = as.matrix(SeuratObject::GetAssayData(object = object,slot = slot)[features, cells, drop = FALSE]))))
    object <- suppressMessages(expr = SeuratObject::StashIdent(object = object, save.name = "ident"))
    groups.use        <- object[[group.by]][cells, , drop = FALSE]
    groups.use.sample <- object[['orgSample']][cells, , drop = FALSE]
    if (sum(rownames(groups.use) == rownames(data0))!=dim(data0)[1]) stop("Error: heatmap plot input data and group levels dimension not match")
    # data0$cell.groups <- groups.use[[group.by]][match(rownames(data0), rownames(groups.use))]
    data0$cell.groups <- groups.use.sample[['orgSample']][match(rownames(data0), rownames(groups.use.sample))]
    if (features.process) {
      data1 = data0
      data2 = data0 %>% dplyr::select('cell.groups') %>% as.data.frame()
      if (length(features.process.fns)==1) features.process.fns = rep(features.process.fns, length(features.process.lists))
      for(i in 1:length(features.process.lists)){
        if (features.process.fns[i]=='/') {
          index.nominator   = grep(paste(sprintf('^%s$', features.process.lists[[i]][1]), collapse = '|'), colnames(data1))
          index.denominator = grep(paste(sprintf('^%s$', features.process.lists[[i]][2]), collapse = '|'), colnames(data1))
          data2[,i+1] = (data1[,index.nominator])/(data1[,index.denominator]+1e10)
        } else if (features.process.fns[i]=='+') {
          data2[,i+1] = rowSums(data1[,grep(paste(sprintf('^%s$', features.process.lists[[i]]), collapse = '|'), colnames(data1))])
        } else{
          print("only '+' or '/' (features.process.fns) functions are supported right now")
        }
      }
      if (!is.null(names(features.process.lists))) {
        features = rev(names(features.process.lists))
        colnames(data2) <- c('cell.groups', names(features.process.lists))
      } else {
        features = rev(paste('featuresCheck', 1:length(features.process.lists), sep = ''))
        colnames(data2) <- c('cell.groups', paste('featuresCheck', 1:length(features.process.lists), sep = ''))
      }
      data0 = data2
      if(debug) print(head(data0))
      if(debug) print(sprintf("features.process features is %s", paste(features, collapse = ', ')))
      if(debug) print(features.process.lists)
    }
    data0 <- data0%>% dplyr::group_by(cell.groups) %>% dplyr::summarise_all(mean) %>% as.data.frame()
    data  <- data0 %>% dplyr::select(-cell.groups)
    rownames(data) <- data0$cell.groups
    plots <- vector(mode = "list", length = ncol(x = groups.use))
  } else if (display == 'expCond.merge') {
    data0 <- as.data.frame(x = as.matrix(x = t(x = as.matrix(SeuratObject::GetAssayData(object = object,slot = slot)[features, cells, drop = FALSE]))))
    object <- suppressMessages(expr = SeuratObject::StashIdent(object = object, save.name = "ident"))
    groups.use        <- object[[group.by]][cells, , drop = FALSE]
    if (sum(rownames(groups.use) == rownames(data0))!=dim(data0)[1]) stop("Error: heatmap plot input data and group levels dimension not match")
    data0$cell.groups <- groups.use[['expCond']][match(rownames(data0), rownames(groups.use))]
    data0 <- data0 %>% dplyr::group_by(cell.groups) %>% dplyr::summarise_all(mean) %>% as.data.frame()
    data  <- data0 %>% dplyr::select(-cell.groups)
    rownames(data) <- data0$cell.groups
    plots <- vector(mode = "list", length = ncol(x = groups.use))
  } else if (display == 'cell'){
    data <- as.data.frame(x = as.matrix(x = t(x = as.matrix(SeuratObject::GetAssayData(object = object,slot = slot)[features, cells, drop = FALSE]))))
    object <- suppressMessages(expr = SeuratObject::StashIdent(object = object, save.name = "ident"))
    if (features.process) {
      data1 = data
      data2 = data.frame(matrix(ncol = length(features.process.lists), nrow = dim(data)[1]))
      if (length(features.process.fns)==1) features.process.fns = rep(features.process.fns, length(features.process.lists))
      for(i in 1:length(features.process.lists)){
        if (features.process.fns[i]=='/') {
          index.nominator   = grep(paste(sprintf('^%s$', features.process.lists[[i]][1]), collapse = '|'), colnames(data1))
          index.denominator = grep(paste(sprintf('^%s$', features.process.lists[[i]][2]), collapse = '|'), colnames(data1))
          data2[,i] = (data1[,index.nominator])/(data1[,index.denominator]+1e10)
        } else if (features.process.fns[i]=='+') {
          data2[,i] = rowSums(data1[,grep(paste(sprintf('^%s$', features.process.lists[[i]]), collapse = '|'), colnames(data1))])
        } else{
          print("only '+' or '/' (features.process.fns) functions are supported right now")
        }
      }
      rownames(data2) = rownames(data)
      if (!is.null(names(features.process.lists))) {
        features = rev(names(features.process.lists))
        colnames(data2) <- names(features.process.lists)
      } else {
        features = rev(paste('featuresCheck', 1:length(features.process.lists), sep = ''))
        colnames(data2) <- paste('featuresCheck', 1:length(features.process.lists), sep = '')
      }
      data = data2
    }
    group.by <- group.by %||% "ident"
    groups.use <- object[[group.by]][cells, , drop = FALSE]
    plots <- vector(mode = "list", length = ncol(x = groups.use))
  }

  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    if (display == 'sample.merge'){
      group.use0  <- merge(groups.use.sample, groups.use, by = 0)
      group.use0  <- group.use0[!duplicated(group.use0$orgSample),] %>% dplyr::select(c('orgSample', 'expCond'))
      group.use  <- group.use0$expCond
    } else if (display == 'expCond.merge') {
      group.use0  <- groups.use %>% dplyr::distinct(expCond, .keep_all = TRUE)
      group.use   <- group.use0$expCond
    } else if (display == 'cell') {
      group.use <- groups.use[, i, drop = TRUE]
    }
    ## ---
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    ## ---
    if (display == 'sample.merge') {
      names(x = group.use) <- group.use0$orgSample
      if (debug) {
        print('77777777')
        print(sprintf("%s group.use levels is %s", length(levels(group.use)), paste(levels(group.use), collapse = ', ')))
      }
      if (!is.null(expCondReorderLevels)) {
        group.use <- factor(x = group.use, levels = expCondReorderLevels)
      }
      if (debug) {
        print('-=-=-=')
        print(sprintf(" %s updated group.use levels is %s", length(levels(group.use)), paste(levels(group.use), collapse = ', ')))
        print('77777777')
      }
      if (draw.lines) {
        lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                  0.0025)
        placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) *
                                             lines.width), FUN = function(x) {
                                               return(SeuratObject::RandomName(length = 20))
                                             })
        placeholder.groups <- rep(x = levels(x = group.use),
                                  times = lines.width)
        group.levels <- levels(x = group.use)
        names(x = placeholder.groups) <- placeholder.cells
        group.use <- as.vector(x = group.use)
        names(x = group.use) <- group.use0$orgSample
        group.use <- factor(x = c(group.use, placeholder.groups),
                            levels = group.levels)
        na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                                ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                             colnames(x = data.group)))
        data.group <- rbind(data.group, na.data.group)
      }
    } else if (display == 'cell') {
      names(x = group.use) <- cells
      if (debug) {
        print('2222222222')
        print(sprintf("%s group.use levels is %s", length(levels(group.use)), paste(levels(group.use), collapse = ', ')))
      }

      if (!is.null(expCondReorderLevels)) {
        group.use <- factor(x = group.use, levels = expCondReorderLevels)
      }

      if(debug) {
        print('-=-=-=')
        print(sprintf(" %s updated group.use levels is %s", length(levels(group.use)), paste(levels(group.use), collapse = ', ')))
        print('2222222222')
      }

      if (draw.lines) {
        lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                  0.0025)
        placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) *
                                             lines.width), FUN = function(x) {
                                               return(SeuratObject::RandomName(length = 20))
                                             })
        placeholder.groups <- rep(x = levels(x = group.use),
                                  times = lines.width)
        group.levels <- levels(x = group.use)
        # print('333333333')
        # print(sprintf("group.use levels is %s", levels(group.use)))
        # if (!is.null(expCondReorderLevels)) {
        #   group.use <- factor(x = group.use, levels = expCondReorderLevels)
        # }
        # print('-=-=-=')
        # print(sprintf("updated group.use levels is %s", levels(group.use)))
        # print('3333333333')
        names(x = placeholder.groups) <- placeholder.cells
        group.use <- as.vector(x = group.use)
        names(x = group.use) <- cells
        group.use <- factor(x = c(group.use, placeholder.groups),
                            levels = group.levels)
        na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                                ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                             colnames(x = data.group)))
        data.group <- rbind(data.group, na.data.group)
      }
    } else if (display == 'expCond.merge') {
      names(x = group.use) <- group.use0$expCond
      if (debug) {
        print('99999999999')
        print(sprintf("%s group.use levels is %s", length(levels(group.use)), paste(levels(group.use), collapse = ', ')))
      }
      if (!is.null(expCondReorderLevels)) {
        group.use <- factor(x = group.use, levels = expCondReorderLevels)
      }
      if (debug) {
        print('-=-=-=')
        print(sprintf(" %s updated group.use levels is %s", length(levels(group.use)), paste(levels(group.use), collapse = ', ')))
        print('99999999999')
      }
      if (draw.lines) {
        lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                  0.0025)
        placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) *
                                             lines.width), FUN = function(x) {
                                               return(SeuratObject::RandomName(length = 20))
                                             })
        placeholder.groups <- rep(x = levels(x = group.use),
                                  times = lines.width)
        group.levels <- levels(x = group.use)
        names(x = placeholder.groups) <- placeholder.cells
        group.use <- as.vector(x = group.use)
        names(x = group.use) <- group.use0$orgSample
        group.use <- factor(x = c(group.use, placeholder.groups),
                            levels = group.levels)
        na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                                ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                             colnames(x = data.group)))
        data.group <- rbind(data.group, na.data.group)
      }

    }

    lgroup <- length(levels(group.use))
    if(res.output) {
      if(is.null(outputFname)) outputFname = file.path(getwd(), 'heatmap_countInput')
      print(sprintf("heatmap count input is also saved in %s", outputFname))
      save(data.group, group.use, file = sprintf('%s.Rdata', outputFname))
    }
    if(debug) print(sprintf("Final plot features is %s", paste(features, collapse = ', ')))
    plot <- SingleRasterMap(data = data.group, raster = raster,
                            disp.min = disp.min, disp.max = disp.max, feature.order = features,
                            cell.order = names(x = sort(x = group.use)), group.by = group.use)
    if (group.bar) {
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      if (!is.null(x = names(x = group.colors))) {
        cols <- unname(obj = group.colors[levels(x = group.use)])
      } else {
        cols <- group.colors[1:length(x = levels(x = group.use))] %||%
          default.colors
      }
      if (any(is.na(x = cols))) {
        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
        cols <- Col2Hex(cols)
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols,
                                                                        start = 1, stop = 7)))))
        through <- length(x = default.colors)
        while (length(x = col.dups) > 0) {
          pal.max <- length(x = col.dups) + through
          cols.extra <- hue_pal()(pal.max)[(through +
                                              1):pal.max]
          cols[col.dups] <- cols.extra
          col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols,
                                                                          start = 1, stop = 7)))))
        }
      }
      group.use2 <- sort(x = group.use)
      if (draw.lines) {
        na.group <- SeuratObject::RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2),
                                    na.group)
        group.use2[placeholder.cells] <- na.group
        cols <- c(cols, "#FFFFFF")
      }
      pbuild <- ggplot_build(plot = plot)
      names(x = cols) <- levels(x = group.use2)
      y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) +
        y.range * 0.015
      y.max <- y.pos + group.bar.height * y.range
      x.min <- min(pbuild$layout$panel_params[[1]]$x.range) +
        0.1
      x.max <- max(pbuild$layout$panel_params[[1]]$x.range) -
        0.1
      plot <- plot + annotation_raster(raster = t(x = cols[group.use2]),
                                       xmin = x.min, xmax = x.max, ymin = y.pos, ymax = y.max) +
        coord_cartesian(ylim = c(0, y.max), clip = "off") +
        scale_color_discrete(name = "Identity", na.translate = FALSE)
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        x.divs <- pbuild$layout$panel_params[[1]]$x.major %||%
          attr(x = pbuild$layout$panel_params[[1]]$x$get_breaks(),
               which = "pos")
        x <- data.frame(group = sort(x = group.use),
                        x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group,
                              FUN = function(y) {
                                if (isTRUE(x = draw.lines)) {
                                  mean(x = y[-length(x = y)])
                                }
                                else {
                                  mean(x = y)
                                }
                              })
        label.x.pos <- data.frame(group = names(x = label.x.pos),
                                  label.x.pos)
        plot <- plot + geom_text(stat = "identity", data = label.x.pos,
                                 aes_string(label = "group", x = "label.x.pos"),
                                 y = y.max + y.max * 0.03 * 0.5, angle = angle,
                                 hjust = hjust, vjust = vjust, size = size)
        plot <- suppressMessages(plot + coord_cartesian(ylim = c(0,
                                                                 y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) *
                                                                   size), clip = "off"))
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots)
  }
  return(plots)
}
##----------------------------------------------------------------------------------------- ##
