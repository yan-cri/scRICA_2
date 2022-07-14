#' getGoiHeatmap() Function
#' @details
#' This function is used to make heat-map of marker/features genes in provided 'goiFname' or geneNames.
#'
#' @param heatmap.view 3 options: 'ident', 'expCond', and 'expCond.ident' to make heatmap, by default 'idents'.
#' @param resDir full path of integration results analysis are saved, where RDS file is saved inside the 'RDS_Dir'. This path is also returned by getClusterMarkers() execution.
#' @param rds User also can provide the full path of RDS file instead of 'resDir' where RDS file is saved in. If this option is used, please also provide 'resDir' to specify where the analysis results will be saved.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param scaled whether input rds is scaled, by default False, will scale data based on selected 'cellcluster' and 'expCond'.
#' @param resSaveFname whether to save re-scaled seurat object on the defined 'cellcluster' and 'expCond', by default it is 'NULL' without saving re-scaled RDS object.
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
#' @param debug by default False, turn on to print extra messages to debug the analysis process
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
#'
#' @keywords Heatmap
#' @examples getGoiHeatmap(rds, geneNames/goiFname, expCondCheck='sample/idents/expCond*')
#' @export
#' @return the Heatmap of provided GOI(gene of interest) provided by 'geneNames' directly or presented in 'goiFname'.
#'
## ---------------------------------------------------------------------------------------
getGoiHeatmap <- function(heatmap.view = 'ident', resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL,
                          scaled = F, resSaveFname = NULL,
                          goiFname = NULL, geneNames = NULL,
                          expCondCheck='sample', expCondCheckFname = NULL,
                          cellcluster = NULL , expCond = NULL, expCondReorderLevels = NULL,
                          plotFnamePrefix='goiHeatmap', plotWidth=25, plotHeight = 20,
                          minVal = -2.5, maxVal = 2.5, fontsize.legend = 20,
                          fontangle.top = 0, fontsize.top = 5, fontsize.y = 20, fontangle.y = 0, hjust.top = 0, vjust.top = 0, barHeight.top = 0.02, debug = F) {
  ## ----
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
    print('-=-=-=-')
    print('updated Idents are as below:')
    print(table(Idents(seuratObjFinal)))
    print('-=-=-=-')
  }
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
  if (debug) print(expCondLevels)
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
  if (!scaled) {
    seuratObjFinal <- Seurat::ScaleData(object = seuratObjFinal, do.scale = T, do.center = T) ## by default, both do.scale/center are on, will scale the expression level for each feature by dividing the centered feature expression levels by their standard deviations
  } else {
    print("Input RDS is already scaled, no further scaling is implemented here.")
  }
  if (!is.null(resSaveFname)) {
    saveRDS(object = seuratObjFinal, file = file.path(rdsDir, sprintf("%s.rds", resSaveFname)) )
  }
  print('Start: Making feature genes heatmap plot')
  if (heatmap.view == 'ident') {
    ## ---
    print('Start to make heatmap seperated on Idents.')
    if (!is.null(expCondReorderLevels)) {
      # print(levels(Idents(seuratObjFinal)))
      if ( !all(expCondReorderLevels %in% levels(Idents(seuratObjFinal))) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort heatmap top legend bars.")
      Seurat::Idents(seuratObjFinal) <- factor(Seurat::Idents(seuratObjFinal), levels = expCondReorderLevels)
    }
    features.heatmap1 <- DoHeatmap2(object = seuratObjFinal, features = features, assay = 'RNA', slot = 'scale.data',
                                    group.by = 'ident', angle = fontangle.top, size = fontsize.top, hjust = hjust.top, vjust =  vjust.top, group.bar.height = barHeight.top,
                                    disp.min = minVal, disp.max = maxVal)
    pdf(file = file.path(plotResDir, sprintf('%s_heatmap.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    print(features.heatmap1+ plotTheme + guides(colour = 'none'))
    dev.off()
    # features.heatmap2 <- Seurat::DoHeatmap(object = seuratObjFinal, features = features, assay = 'RNA', slot = 'data', group.by = 'ident')
    # pdf(file = file.path(plotResDir, sprintf('%s_heatmap_norm.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    # print(features.heatmap2)
    # dev.off()
    ## ---
  } else  if (heatmap.view == 'expCond') {
    ## ---
    print('Start to make heatmap seperated on expCond')
    if (!is.null(expCondReorderLevels)) {
      # print(levels(factor(seuratObjFinal$expCond)))
      if ( !all(expCondReorderLevels %in% levels(factor(seuratObjFinal$expCond))) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort heatmap top legend bars.")
      seuratObjFinal@meta.data$expCond <- factor(seuratObjFinal@meta.data$expCond, levels = expCondReorderLevels)
    }
    features.heatmap1 <- DoHeatmap2(object = seuratObjFinal, features = features, assay = 'RNA', slot = 'scale.data',
                                    group.by = 'expCond', angle = fontangle.top, size = fontsize.top, hjust = hjust.top, vjust =  vjust.top, group.bar.height = barHeight.top,
                                    disp.min = minVal, disp.max = maxVal)
    pdf(file = file.path(plotResDir, sprintf('%s_heatmap.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    print(features.heatmap1+ plotTheme + guides(colour = 'none'))
    dev.off()
    # features.heatmap2 <- Seurat::DoHeatmap(object = seuratObjFinal, features = features, assay = 'RNA', slot = 'data', group.by = 'expCond')
    # pdf(file = file.path(plotResDir, sprintf('%s_heatmap_norm.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    # print(features.heatmap2)
    # dev.off()
    ## ---
  } else if (heatmap.view == 'expCond.ident') {
    ## ---
    print('Start to make heatmap seperated on expCond and idents combinations.')
    seuratObjFinal <- AddMetaData(object = seuratObjFinal, metadata = paste(seuratObjFinal$expCond, Idents(seuratObjFinal), sep = ' '), col.name = 'expIdentComb')
    ## -
    if (!is.null(expCondReorderLevels)) {
      if ( !all(expCondReorderLevels %in% levels(factor(seuratObjFinal$expIdentComb))) ) stop("Please provide correct corresponding 'expCondReorderLevels' to sort heatmap top legend bars.")
      seuratObjFinal@meta.data$expIdentComb <- factor(seuratObjFinal@meta.data$expIdentComb, levels = expCondReorderLevels)
    }
    ## -
    features.heatmap1 <- DoHeatmap2(object = seuratObjFinal, features = features, assay = 'RNA', slot = 'scale.data',
                                    group.by = 'expIdentComb', angle = fontangle.top, size = fontsize.top, hjust = hjust.top, vjust =  vjust.top, group.bar.height = barHeight.top,
                                    disp.min = minVal, disp.max = maxVal)
    pdf(file = file.path(plotResDir, sprintf('%s_heatmap.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    print(features.heatmap1+ plotTheme + guides(colour = 'none'))
    dev.off()
    # features.heatmap2 <- Seurat::DoHeatmap(object = seuratObjFinal, features = features, assay = 'RNA', slot = 'data', group.by = 'expIdentComb')
    # pdf(file = file.path(plotResDir, sprintf('%s_heatmap_norm.pdf', plotFnamePrefix )), width = plotWidth, height = plotHeight)
    # print(features.heatmap2)
    # dev.off()
    ## ---
  }
  print('END: Making selected genes heatmap plot')
  print('********************')
  ##--------------------------------------------------------------------------------------##
}

##----------------------------------------------------------------------------------------- ##
## DoHeatmap2 adding 'vjust' to adjust top color legend bar text position from Seurat::DoHeatmap
DoHeatmap2 <- function (object, features = NULL, cells = NULL, group.by = "ident",
          group.bar = TRUE, group.colors = NULL, disp.min = -2.5, disp.max = NULL,
          slot = "scale.data", assay = NULL, label = TRUE, size = 5.5,
          hjust = 0, vjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE,
          lines.width = NULL, group.bar.height = 0.02, combine = TRUE)
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object,
                                                 slot = slot))
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
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object,
                                                             slot = slot)[features, cells, drop = FALSE])))
  object <- suppressMessages(expr = StashIdent(object = object,
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  plots <- vector(mode = "list", length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) *
                                           lines.width), FUN = function(x) {
                                             return(RandomName(length = 20))
                                           })
      placeholder.groups <- rep(x = levels(x = group.use),
                                times = lines.width)
      group.levels <- levels(x = group.use)
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
    lgroup <- length(levels(group.use))
    plot <- SingleRasterMap(data = data.group, raster = raster,
                            disp.min = disp.min, disp.max = disp.max, feature.order = features,
                            cell.order = names(x = sort(x = group.use)), group.by = group.use)
    if (group.bar) {
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      if (!is.null(x = names(x = group.colors))) {
        cols <- unname(obj = group.colors[levels(x = group.use)])
      }
      else {
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
        na.group <- RandomName(length = 20)
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
                                 hjust = hjust, size = size, vjust = vjust)
        plot <- suppressMessages(plot + coord_cartesian(ylim = c(0,
                                                                 y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) *
                                                                   size), clip = "off"))
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}
##----------------------------------------------------------------------------------------- ##
