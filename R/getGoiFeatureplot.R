#' getGoiFeatureplot() Function
#' @details
#' This function is used to make dotplot of marker/features genes in provided 'goiFname'.
#'
#' @param resDir specify an exiting full path of directory, where results will be saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved with above rdsDir option.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param expCondCheckFname suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param cellcluster optional, if needed, this option can be used to specify cell clusters to be displayed on the feature plot
#' @param goiFname path to file, where a list of marker/features genes are provided in column 'Gene', if column 'Cell Type' is also provided, option 'geneCellTypeOrder' can be used to adjust orders.
#' @param featurePlotMinExpCutoff feature plot minimum expression value threshold, by default = 0.3.
#' @param featurePlotReductionMethod feature plot clustering reduction methods 'tsne' or 'umap', by default 'umap'.
#' @param featurePlotLabel whether to label feature plot cell types, by default True.
#' @param featurePlotFnamePrefix feature plot file name prefix, if not defined, by default 'goiFeaturePlot'.
#'
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 element_rect
#' @importFrom cowplot theme_cowplot
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat NoLegend
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat FeaturePlot
#' @importFrom Seurat FetchData
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom utils write.table
#' @importFrom tools file_ext
#' @importFrom utils read.delim
#' @importFrom xlsx read.xlsx
#'
#' @keywords GoiFeatureplot
#' @examples getGoiFeatureplot(rds, goiFname, expCondCheck='sample/idents/expCond*')
#' @export
#' @return the dotplots of provided GOI(gene of interest) saved in '' inside the provided 'resDir'
## ---------------------------------------------------------------------------------------
getGoiFeatureplot <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL,  expCondCheck='sample', expCondCheckFname = NULL, expCondReorderLevels = NULL, cellcluster = NULL, goiFname, featurePlotMinExpCutoff=0.3, featurePlotMaxExpCutoff=NULL, featurePlotReductionMethod='umap', featurePlotLabel = T, featurePlotFnamePrefix='goiFeaturePlot' ){
  ##--------------------------------------------------------------------------------------##
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) print("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  ##--------------------------------------------------------------------------------------##
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
  ##--------------------------------------------------------------------------------------##
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
    source(as.character(newAnnotationRscriptName))
  }
  ##--------------------------------------------------------------------------------------##
  # setup custom theme for plotting
  theme1noLegend <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                          # legend.key.size = unit(0.7, "cm"),
                          axis.title = element_text(size = 20),
                          axis.text = element_text(size = 25),
                          legend.position="bottom",
                          legend.text = element_text(size = 14) ) + NoLegend()
  theme1wLegend <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                         # legend.key.size = unit(0.7, "cm"),
                         axis.title = element_text(size = 15),
                         axis.text = element_text(size = 20),
                         legend.position="bottom",
                         legend.text = element_text(size = 15) )
  theme1wLegendRight <- theme(plot.title = element_text(size = 16, hjust = 0.5),
                              # legend.key.size = unit(0.7, "cm"),
                              axis.title = element_text(size = 20),
                              axis.text = element_text(size = 25),
                              legend.position="right",
                              legend.text = element_text(size = 14) )
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
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding 'plotResDir' for feature-plot to save
  plotResDir            <- paste(resDir, sprintf('featurePlot_selected_markers_%s', expCondCheckFname), sep = '/')
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
  ## ---
  plotResDir            <- paste(plotResDir, featurePlotFnamePrefix, sep = '/')
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
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  ##--------------------------------------------------------------------------------------##
  # print(table(Seurat::Idents(seuratObjFinal)))
  clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% clusterLevels ) ) stop('Please provide the corresponding cell clusters ')
    print(sprintf('Subsetting %s specific cell clusters: %s', length(cellcluster), paste(cellcluster, collapse = ',')))
    seuratObjFinal                          <- subset(seuratObjFinal, idents = cellcluster )
  }
  ##--------------------------------------------------------------------------------------##
  print('=====================================================')
  print(sprintf('START feature plot, the plots will be save at %s', plotResDir))
  print(table(Idents(seuratObjFinal)))
  if (file_ext(basename(goiFname)) == 'xlsx') {
    markerGenesPrep       <- read.xlsx(file = as.character(goiFname), sheetIndex = 1, header = T)
  } else if (file_ext(basename(goiFname)) == 'txt') {
    markerGenesPrep       <- read.delim(file = as.character(goiFname), header = T, sep = '\t')
  }
  print("----------------")
  print(sprintf("A total of %s genes will be plotted", length(unique(markerGenesPrep$Gene))))
  if (sum(duplicated(markerGenesPrep$Gene))>0) print(sprintf("%s genes are duplicated genes, they are: %s", sum(duplicated(markerGenesPrep$Gene)), paste(markerGenesPrep$Gene[duplicated(markerGenesPrep$Gene)], collapse = ', ' ) ))
  print("----------------")
  ## ---
  if ("Gene" %in% colnames(markerGenesPrep)) {
    markerGenes           <- as.character(unique(markerGenesPrep$Gene))
  } else {
    markerGenes           <- as.character(unique(markerGenesPrep[,1]))
  }
  ## -
  markerGenes = markerGenes[!is.na(markerGenes)]
  ## -
  if ('geneType' %in% colnames(markerGenesPrep) & "Gene" %in% colnames(markerGenesPrep)) {
    markerGenesPrep2      <- markerGenesPrep %>% dplyr::select(c('Gene', 'geneType'))
    geneTypes             <- markerGenesPrep2[match(markerGenes, markerGenesPrep2$Gene),]
    colnames(geneTypes) <- c('Gene', 'geneType')
  } else if (!'geneType' %in% colnames(markerGenesPrep) & "Gene" %in% colnames(markerGenesPrep)) {
    markerGenesPrep2      <- markerGenesPrep %>% dplyr::select(grep('Gene', colnames(markerGenesPrep)))
    geneTypes             <- as.data.frame(markerGenesPrep2[match(markerGenes, markerGenesPrep2$Gene),])
    colnames(geneTypes) <- c('Gene')
  } else if ('geneType' %in% colnames(markerGenesPrep) & !"Gene" %in% colnames(markerGenesPrep)) {
    markerGenesPrep2      <- markerGenesPrep %>% dplyr::select(c(1, grep('geneType', colnames(markerGenesPrep))))
    geneTypes             <- markerGenesPrep2[match(markerGenes, markerGenesPrep2[,1]),]
    colnames(geneTypes) <- c('Gene', 'geneType')
  } else {
    markerGenesPrep2      <- markerGenesPrep[,1]
    geneTypes             <- markerGenesPrep2[match(markerGenes, markerGenesPrep2[,1]), 1]
    colnames(geneTypes) <- c('Gene')
  }
  ##--------------------------------------------------------------------------------------##
  ## 2. make featurePlot for each above selected marker/features
  Seurat::DefaultAssay(seuratObjFinal)   <- "RNA"
  # Seurat::DefaultAssay(seuratObjFinal) <- "integrated"
  # featurePlotMinExpCutoff        = 0.3
  # featurePlotReductionMethod     = 'umap'
  # # featurePlotReductionMethod     = 'tsne'
  if (expCondCheck == 'idents') {
    expValSummary       <- matrix(NA, ncol = 6, nrow = length(markerGenes))
  } else {
    maxExpValSummary    <- matrix(NA, ncol = length(levels(as.factor(seuratObjFinal@meta.data$expCond))), nrow = length(markerGenes))
  }
  ##--------------------------------------------------------------------------------------##
  for (i in 1:length(markerGenes)) {
    print(sprintf('%s. START feature plot for gene marker %s', i, markerGenes[i]))
    # if (seuratObjFinal@active.assay == 'RNA' & markerGenes[i] %in% rownames(seuratObjFinal@assays$RNA@data)) {
    if ( markerGenes[i] %in% rownames(seuratObjFinal@assays$RNA@data) ) {
      ## if max gene expression value >
      if ( max(Seurat::FetchData(seuratObjFinal, markerGenes[i]))<featurePlotMinExpCutoff ) {
        print(sprintf("maximum %s expression values is less than overall defined 'featurePlotMinExpCutoff'=%s, no feature plot is generated", markerGenes[i], featurePlotMinExpCutoff))
      } else {
        if (expCondCheck == 'idents') {
          # -
          maxExpVals      <- summary(Seurat::FetchData(seuratObjFinal, markerGenes[i]))
          print(sprintf('summary expression values for gene marker %s at combined expCond are as below: ', markerGenes[i] ))
          print(maxExpVals)
          print('-=-=-=')
          ## -
          expValSummary[i,] <- sapply(strsplit(summary(Seurat::FetchData(seuratObjFinal, markerGenes[i])), split = ':'), '[[', 2)
          ## -
          if (is.null(featurePlotMaxExpCutoff)) {
            if (featurePlotLabel) {
              featurePlot     <- FeaturePlot2(object = seuratObjFinal, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                              order = T, label = T, pt.size = 0.5, ncol = 1, label.size = 5,
                                              # split.by = 'expCond',
                                              cols = c('#D3D3D3', '#CC0000'),
                                              min.cutoff = featurePlotMinExpCutoff )
            } else {
              featurePlot     <- FeaturePlot2(object = seuratObjFinal, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                              order = T, label = F, pt.size = 0.5, ncol = 1, label.size = 5,
                                              # split.by = 'expCond',
                                              cols = c('#D3D3D3', '#CC0000'),
                                              min.cutoff = featurePlotMinExpCutoff )
            }

          } else {
            if (featurePlotLabel) {
              featurePlot     <- FeaturePlot2(object = seuratObjFinal, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                              order = T, label = T, pt.size = 0.5, ncol = 1, label.size = 5,
                                              # split.by = 'expCond',
                                              cols = c('#D3D3D3', '#CC0000'),
                                              min.cutoff = featurePlotMinExpCutoff,
                                              max.cutoff = featurePlotMaxExpCutoff)
            } else {
              featurePlot     <- FeaturePlot2(object = seuratObjFinal, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                              order = T, label = F, pt.size = 0.5, ncol = 1, label.size = 5,
                                              # split.by = 'expCond',
                                              cols = c('#D3D3D3', '#CC0000'),
                                              min.cutoff = featurePlotMinExpCutoff,
                                              max.cutoff = featurePlotMaxExpCutoff)
            }

          }

          if ('geneType' %in% colnames(geneTypes)) {
            if (featurePlotLabel) {
              pdf(file = file.path(sprintf('%s/%s_featurePlot_%s_%s_expCondComb.pdf', plotResDir, featurePlotReductionMethod, gsub('/', '_', geneTypes$geneType[match(markerGenes[i], geneTypes$Gene)]), markerGenes[i])), width = 5.5, height = 6)
            } else {
              pdf(file = file.path(sprintf('%s/%s_featurePlot_woLable_%s_%s_expCondComb.pdf', plotResDir, featurePlotReductionMethod, gsub('/', '_', geneTypes$geneType[match(markerGenes[i], geneTypes$Gene)]), markerGenes[i])), width = 5.5, height = 6)
            }

          } else {
            if (featurePlotLabel) {
              pdf(file = file.path(sprintf('%s/%s_featurePlot_%s_expCondComb.pdf', plotResDir, featurePlotReductionMethod, markerGenes[i])), width = 5.5, height = 6)
            } else {
              pdf(file = file.path(sprintf('%s/%s_featurePlot_woLable_%s_expCondComb.pdf', plotResDir, featurePlotReductionMethod, markerGenes[i])), width = 5.5, height = 6)
            }

          }
          print(featurePlot+theme1wLegend)
          dev.off()
          ## -
        } else {
          ## ---
          ## separate feature plots by updated expCond except for expCondCheck == 'idents'
          if (!is.null(expCondReorderLevels)){
            seuratObjFinal@meta.data$expCond <- factor(seuratObjFinal@meta.data$expCond, levels = expCondReorderLevels)
          }
          expCondLevels <- levels(as.factor(seuratObjFinal@meta.data$expCond))
          seuratSubset         <- list()
          for (l in 1:length(expCondLevels)) seuratSubset[[l]] <- Seurat::FetchData(subset(seuratObjFinal, subset = expCond == as.character(expCondLevels[l])), markerGenes[i])
          names(seuratSubset)  <- expCondLevels
          maxExpVals           <- sapply(seuratSubset, function(x) round(max(x), digits = 2))
          print(sprintf('max expression values for gene marker %s at %s levels expCond are as below: ', markerGenes[i], length(expCondLevels) ))
          print(maxExpVals)
          print('-=-=-=')
          maxExpValSummary[i,] <- maxExpVals
          if (any(maxExpVals < featurePlotMinExpCutoff)) {
            if (length(which(maxExpVals < featurePlotMinExpCutoff)) != 1) {
              ## -
              for (r in 1:length(which(maxExpVals < featurePlotMinExpCutoff))){
                if (r ==1) {
                  seuratObjFinalUpdate <- subset(seuratObjFinal, subset = expCond != expCondLevels[which(maxExpVals < featurePlotMinExpCutoff)][r] )
                } else {
                  seuratObjFinalUpdate <- subset(seuratObjFinalUpdate, subset = expCond != expCondLevels[which(maxExpVals < featurePlotMinExpCutoff)][r] )
                }
              }
              ## -
            } else {
              seuratObjFinalUpdate <- subset(seuratObjFinal, subset = expCond != as.character(expCondLevels[which(maxExpVals < featurePlotMinExpCutoff)]))
            }

          } else {
            seuratObjFinalUpdate <- seuratObjFinal
          }
          expCondLevelsUpdate  <- levels(as.factor(seuratObjFinalUpdate@meta.data$expCond))
          if (length(expCondLevelsUpdate) < length(expCondLevels)) print(sprintf('%s expCond (%s) is/are removed due to low expression value for gene marker %s', length(which(maxExpVals < featurePlotMinExpCutoff)), paste(as.character(expCondLevels[which(maxExpVals < featurePlotMinExpCutoff)]), collapse = ', '), markerGenes[i]))
          ## -
          if (is.null(featurePlotMaxExpCutoff)) {
            if (featurePlotLabel) {
              featurePlotExpSplit  <- FeaturePlot(object = seuratObjFinalUpdate, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                                  order = T, label = T, pt.size = 0.5, ncol = 1, label.size = 5,
                                                  split.by = 'expCond',
                                                  cols = c('#D3D3D3', '#CC0000'),
                                                  min.cutoff = featurePlotMinExpCutoff )
            } else {
              featurePlotExpSplit  <- FeaturePlot(object = seuratObjFinalUpdate, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                                  order = T, label = F, pt.size = 0.5, ncol = 1, label.size = 5,
                                                  split.by = 'expCond',
                                                  cols = c('#D3D3D3', '#CC0000'),
                                                  min.cutoff = featurePlotMinExpCutoff )
            }

          } else {
            if (featurePlotLabel){
              featurePlotExpSplit  <- FeaturePlot(object = seuratObjFinalUpdate, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                                  order = T, label = T, pt.size = 0.5, ncol = 1, label.size = 5,
                                                  split.by = 'expCond',
                                                  cols = c('#D3D3D3', '#CC0000'),
                                                  min.cutoff = featurePlotMinExpCutoff,
                                                  max.cutoff = featurePlotMaxExpCutoff)
            } else {
              featurePlotExpSplit  <- FeaturePlot(object = seuratObjFinalUpdate, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                                  order = T, label = F, pt.size = 0.5, ncol = 1, label.size = 5,
                                                  split.by = 'expCond',
                                                  cols = c('#D3D3D3', '#CC0000'),
                                                  min.cutoff = featurePlotMinExpCutoff,
                                                  max.cutoff = featurePlotMaxExpCutoff)
            }

          }
          if ('geneType' %in% colnames(geneTypes)) {
            if (featurePlotLabel){
              pdf(file = file.path(sprintf('%s/%s_featurePlot_%s_%s_%s.pdf', plotResDir, featurePlotReductionMethod, gsub('/', '_', geneTypes$geneType[match(markerGenes[i], geneTypes$Gene)]), markerGenes[i], expCondCheckFname)), width = 5.5*length(expCondLevelsUpdate), height = 6)
            } else {
              pdf(file = file.path(sprintf('%s/%s_featurePlot_woLabel_%s_%s_%s.pdf', plotResDir, featurePlotReductionMethod, gsub('/', '_', geneTypes$geneType[match(markerGenes[i], geneTypes$Gene)]), markerGenes[i], expCondCheckFname)), width = 5.5*length(expCondLevelsUpdate), height = 6)
            }

          } else {
            if (featurePlotLabel) {
              pdf(file = file.path(sprintf('%s/%s_featurePlot_%s_%s.pdf', plotResDir, featurePlotReductionMethod, markerGenes[i], expCondCheckFname)), width = 5.5*length(expCondLevelsUpdate), height = 6)
            } else {
              pdf(file = file.path(sprintf('%s/%s_featurePlot_woLabel_%s_%s.pdf', plotResDir, featurePlotReductionMethod, markerGenes[i], expCondCheckFname)), width = 5.5*length(expCondLevelsUpdate), height = 6)
            }

          }
          print(featurePlotExpSplit+theme1wLegend)
          dev.off()
        }
        ## ---
      }
      ## -
    } else {
      print(sprintf('%s. No featureplot for gene marker %s, which cannot be found in this data', i, markerGenes[i]))
      print('======')
    }
    print(sprintf('%s. Complete featureplot for gene marker %s', i, markerGenes[i]))
    print('======')
  }
  ##--------------------------------------------------------------------------------------##
  if (expCondCheck == 'idents') {
    colnames(expValSummary) <- sapply(strsplit(summary(Seurat::FetchData(seuratObjFinal, markerGenes[i])), split = ':'), '[[', 1)
    rownames(expValSummary) <- markerGenes
    print(head(expValSummary))
    write.table(x = expValSummary, file = file.path(sprintf('%s/%s_expValSummary_%s.txt', plotResDir, featurePlotReductionMethod, expCondCheckFname)), quote = F, col.names = NA, row.names = T, sep = '\t' )
  } else {
    print(head(maxExpValSummary))
    colnames(maxExpValSummary) <- names(maxExpVals)
    rownames(maxExpValSummary) <- markerGenes
    write.table(x = maxExpValSummary, file = file.path(sprintf('%s/%s_expMaxValSummary_%s.txt', plotResDir, featurePlotReductionMethod, expCondCheckFname)), quote = F, col.names = NA, row.names = T, sep = '\t' )
  }
  print(sprintf('END feature plot, the plots were saved at %s', plotResDir))
  print('=========')
  ##--------------------------------------------------------------------------------------##
}

##---
## update Seurat::FeaturePlot into FeaturePlot2 on low ('min.cutoff') and high 'max.cutoff'
FeaturePlot2 <- function (object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {c("lightgrey", "#ff0000", "#00ff00")} else { c("lightgrey", "blue")}, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA, reduction = NULL, split.by = NULL, keep.scale = "feature", shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5, label = FALSE, label.size = 4, label.color = "black", repel = FALSE, ncol = NULL, coord.fixed = FALSE, by.col = TRUE, sort.cell = NULL, interactive = FALSE, combine = TRUE, raster = NULL, raster.dpi = c(512, 512)) {
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ",
            "parameter instead for equivalent functionality.",
            call. = FALSE, immediate. = TRUE)
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(IFeaturePlot(object = object, feature = features[1],
                        dims = dims, reduction = reduction, slot = slot))
  }
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature",
                                                        "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",size = 14, margin = margin(r = 7)))
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
                   `0` = {
                     warning("No colors provided, using default colors",
                             call. = FALSE, immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors",
                             call. = FALSE, immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '",
                             default.colors[1], "' for double-negatives",
                             call. = FALSE, immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three",
                             call. = FALSE, immediate. = TRUE)
                     cols[1:3]
                   })
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(dims, "ident",
                                              features), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ", paste(features,
                                                              collapse = ", "), " in slot ", slot, call. = FALSE)
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[,
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[,
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
  ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data),
                                     FUN = function(index) {
                                       data.feature <- as.vector(x = data[, index])
                                       min.use <- SetQuantile(cutoff = min.cutoff[index -
                                                                                    3], data.feature)
                                       max.use <- SetQuantile(cutoff = max.cutoff[index -
                                                                                    3], data.feature)
                                       data.feature[data.feature < min.use] <- min.use
                                       data.feature[data.feature > max.use] <- max.use
                                       if (brewer.gran == 2) {
                                         return(data.feature)
                                       }
                                       data.cut <- if (all(data.feature == 0)) {
                                         0
                                       }
                                       else {
                                         as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                                                                          breaks = brewer.gran)))
                                       }
                                       return(data.cut)
                                     })
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  }  else {
    switch(EXPR = split.by, ident = Idents(object = object)[cells,
                                                            drop = TRUE], object[[split.by, drop = TRUE]][cells,
                                                                                                          drop = TRUE])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend,
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[,
                                                                   dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[,
                                                               dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold,
                                negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ],
                   as.vector(x = color.matrix))
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident,
                      , drop = FALSE]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[,
                                                       features]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ",
             paste(no.expression, collapse = ", "), call. = FALSE)
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")],
                         BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[,
                                                              feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident", feature, shape.by)]
      plot <- SingleDimPlot(data = data.single, dims = dims,
                            col.by = feature, order = order, pt.size = pt.size,
                            cols = cols.use, shape.by = shape.by, label = FALSE,
                            raster = raster, raster.dpi = raster.dpi) + scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) + theme_cowplot() +
        CenterTitle()
      if (label) {
        plot <- LabelClusters(plot = plot, id = "ident",
                              repel = repel, size = label.size, color = label.color)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA,
                                                         colour = "black"))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        }
        else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident),
                                                                    limits = ylims) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(),
                               axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(),
                               axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                               axis.title.x = element_blank())
        }
      } else {
        plot <- plot + labs(title = feature)
      }

      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (",
                    unique.feature.exp, ") of ", feature, ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, guide = "colorbar", limits = c(min.cutoff, max.cutoff) ))
        }
      }

      if (!(is.null(x = keep.scale)) && keep.scale == "feature" &&
          !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        print(sprintf("max.cutoff is %s, min.cutoff is %s", max.cutoff, min.cutoff))
        if (is.na(min.cutoff) & is.na(max.cutoff)) {
          plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
        } else {
          plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, limits = c(min.cutoff, max.cutoff)))
        }

      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(x = plots,
                                              values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                                                                                                                                  1, yes = levels(x = data$split)[ii], no = "")),
                                                                                              expand = c(0, 0)) + labs(x = features[1], y = features[2],
                                                                                                                       title = if (ii == 1) {
                                                                                                                         paste("Color threshold:", blend.threshold)
                                                                                                                       } else {
                                                                                                                         NULL
                                                                                                                       }) + no.right), after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol,
                 no = length(x = features))
  legend <- if (blend) {
    "none"
  }
  else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() +
                                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""),
                                                                   limits = ylims) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) +
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]),
                                                            limits = ylims) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots)%%length(x = features) ==
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      }
      else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots),
                                                          f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    }
    else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
                            length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" &&
        !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      # print(sprintf("22222. max.cutoff is %s, min.cutoff is %s", max.cutoff, min.cutoff)) ## for debug
      if (is.na(min.cutoff) & is.na(max.cutoff)) {
        plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, limits = c(min.feature.value, max.feature.value)))
      } else {
        plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, limits = c(min.cutoff, max.cutoff)))
      }
    }
  }
  return(plots)
}
## ---------------------------------------------------------------------------------------
