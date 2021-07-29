#' getGoiFeatureplot() Function
#' @details
#' This function is used to make dotplot of marker/features genes in provided 'goiFname'.
#'
#' @param resDir full path of integration results analysis returned in getClusterMarkers()
#' @param newAnnotation logical value, whether to provide manual annotation
#' @param newAnnotationRscriptName if newAnnotation == T, this script is used to redefine the old clusters
#' @param expCondSepName character string, user defined name either to be 'org' or any character string
#' @param expCondName2change if above 'expCondSepName' is defined not as 'org', provide the name to be changed
#' @param goiFname full path of a file name, where a list of marker/features genes provided
#' @param featurePlotMinExpCutoff featurePlot miniumn expression threshold
#' @param featurePlotReductionMethod featurePlot reduction methods: options are 'tsne' or 'umap'
#' @param featurePlotFnamePrefix featurePlot file name prefix
#'
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 labs
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
#' @examples getGoiFeatureplot()
#' @export
#' @return the dotplots of provided GOI(gene of interest) saved in '' inside the provided 'resDir'
## ---------------------------------------------------------------------------------------
getGoiFeatureplot <- function(resDir=NULL, rdsFname=NULL, newAnnotation=F, newAnnotationRscriptName=NULL,  expCondCheck='sample', expCondName2change=NULL, goiFname, featurePlotMinExpCutoff=0.3, featurePlotReductionMethod='umap', featurePlotFnamePrefix='goiFeaturePlot' ){
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
  print('Done for RDS readin')
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
  # setup cutome theme for plotting
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
  ## -------------------------------------------------------------------------------------
  if (expCondCheck == 'sample') {
    expCondSepName        <- 'expCond_sample'
  } else {
    expCondSepName        <- as.character(expCondCheck)
  }
  ## -------------------------------------------------------------------------------------
  ## update 'seuratObjFinal@meta.data$expCond' and create corresponding 'plotResDir' for feature-plot to save
  plotResDir            <- paste(resDir, sprintf('featurePlot_selected_markers_%s', expCondSepName), sep = '/')
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
  ## ---
  plotResDir            <- paste(plotResDir, featurePlotFnamePrefix, sep = '/')
  if (!dir.exists(plotResDir)) dir.create(plotResDir)
  ## ------
  if (expCondCheck == 'sample') {
    seuratObjFinal        <- seuratObjFinal
  } else if (expCondCheck == 'comb') {
    seuratObjFinal@meta.data$expCond   <- Idents(seuratObjFinal)
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
  print('=====================================================')
  print(sprintf('START feature plot, the plots will be save at %s', plotResDir))
  print(table(Idents(seuratObjFinal)))
  if (file_ext(basename(goiFname)) == 'xlsx') {
    markerGenesPrep       <- read.xlsx(file = as.character(goiFname), sheetIndex = 1, header = T)
  } else if (file_ext(basename(goiFname)) == 'txt') {
    markerGenesPrep       <- read.delim(file = as.character(goiFname), header = T, sep = '\t')
  }
  print("----------------")
  print(sprintf("A total of %s genes will be ploted", length(unique(markerGenesPrep$Gene))))
  if (sum(duplicated(markerGenesPrep$Gene))>0) print(sprintf("%s genes are duplicated genes, they are: %s", sum(duplicated(markerGenesPrep$Gene)), paste(markerGenesPrep$Gene[duplicated(markerGenesPrep$Gene)], collapse = ', ' ) ))
  print("----------------")
  ## ---
  if ("Gene" %in% colnames(markerGenesPrep)) {
    markerGenes           <- as.character(unique(markerGenesPrep$Gene))
  } else {
    markerGenes           <- as.character(unique(markerGenesPrep[,1]))
  }
  markerGenes = markerGenes[!is.na(markerGenes)]
  ## ---
  ## 2. make featurePlot for each above selected marker/features
  Seurat::DefaultAssay(seuratObjFinal)   <- "RNA"
  # Seurat::DefaultAssay(seuratObjFinal) <- "integrated"
  # featurePlotMinExpCutoff        = 0.3
  # featurePlotReductionMethod     = 'umap'
  # # featurePlotReductionMethod     = 'tsne'
  if (expCondSepName == 'comb') {
    expValSummary       <- matrix(NA, ncol = 6, nrow = length(markerGenes))
  } else {
    maxExpValSummary    <- matrix(NA, ncol = length(levels(as.factor(seuratObjFinal@meta.data$expCond))), nrow = length(markerGenes))
  }
  ## ---
  for (i in 1:length(markerGenes)) {
    print(sprintf('%s. START featureplot for gene marker %s', i, markerGenes[i]))
    # if (seuratObjFinal@active.assay == 'RNA' & markerGenes[i] %in% rownames(seuratObjFinal@assays$RNA@data)) {
    if ( markerGenes[i] %in% rownames(seuratObjFinal@assays$RNA@data) ) {
      ## if max gene expression value >
      if ( max(Seurat::FetchData(seuratObjFinal, markerGenes[i]))<featurePlotMinExpCutoff ) {
        print(sprintf("maximun %s expression values is less than overall defined 'featurePlotMinExpCutoff'=%s, no feature plot is generated", markerGenes[i], featurePlotMinExpCutoff))
      } else {
        if (expCondSepName == 'comb') {
          # -
          maxExpVals      <- summary(Seurat::FetchData(seuratObjFinal, markerGenes[i]))
          print(sprintf('summary expression values for gene marker %s at combined expCond are as below: ', markerGenes[i] ))
          print(maxExpVals)
          print('-=-=-=')
          ## -
          expValSummary[i,] <- sapply(strsplit(summary(Seurat::FetchData(seuratObjFinal, markerGenes[i])), split = ':'), '[[', 2)
          ## -
          featurePlot     <- Seurat::FeaturePlot(object = seuratObjFinal, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                                 order = T, label = T, pt.size = 0.5, ncol = 1, label.size = 5,
                                                 # split.by = 'expCond',
                                                 cols = c('#D3D3D3', '#CC0000'),
                                                 min.cutoff = featurePlotMinExpCutoff )
          pdf(file = file.path(sprintf('%s/%s_featurePlot_%s_expCondComb.pdf', plotResDir, featurePlotReductionMethod, markerGenes[i])), width = 5.5, height = 6)
          print(featurePlot+theme1wLegend)
          dev.off()
          ## -
        } else {
          ## ---
          ## separate feature plots by updated expCond except for expCondSepName == 'comb'
          expCondLevels        <- levels(as.factor(seuratObjFinal@meta.data$expCond))
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
          featurePlotExpSplit  <- Seurat::FeaturePlot(object = seuratObjFinalUpdate, features = markerGenes[i], reduction = featurePlotReductionMethod,
                                                      order = T, label = T, pt.size = 0.5, ncol = 1, label.size = 5,
                                                      split.by = 'expCond',
                                                      cols = c('#D3D3D3', '#CC0000'),
                                                      min.cutoff = featurePlotMinExpCutoff )
          pdf(file = file.path(sprintf('%s/%s_featurePlot_%s_expCond%s.pdf', plotResDir, featurePlotReductionMethod, markerGenes[i], expCondSepName)), width = 5.5*length(expCondLevelsUpdate), height = 6)
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
  ## ---
  if (expCondSepName == 'comb') {
    colnames(expValSummary) <- sapply(strsplit(summary(Seurat::FetchData(seuratObjFinal, markerGenes[i])), split = ':'), '[[', 1)
    rownames(expValSummary) <- markerGenes
    print(head(expValSummary))
    write.table(x = expValSummary, file = file.path(sprintf('%s/%s_expValSummary_expCond%s.txt', plotResDir, featurePlotReductionMethod, expCondSepName)), quote = F, col.names = NA, row.names = T, sep = '\t' )
  } else {
    print(head(maxExpValSummary))
    colnames(maxExpValSummary) <- names(maxExpVals)
    rownames(maxExpValSummary) <- markerGenes
    write.table(x = maxExpValSummary, file = file.path(sprintf('%s/%s_expMaxValSummary_expCond%s.txt', plotResDir, featurePlotReductionMethod, expCondSepName)), quote = F, col.names = NA, row.names = T, sep = '\t' )
  }
  print(sprintf('END feature plot, the plots were saved at %s', plotResDir))
  print('=========')
  ## ---------
}
## ---------------------------------------------------------------------------------------
