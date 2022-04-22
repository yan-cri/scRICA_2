## getExpCondClusterPseudotime(): perform functional pseudotime analysis on cell clusters specified from different experimental condition types##
## Developed by Yan Li, July, 2021
##--------------------------------------------------------------------------------------##
#' getExpCondClusterPseudotime() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param resDir full path of integration results analysis are saved, where RDS file is saved inside the 'RDS_Dir'. This path is also returned by getClusterMarkers() execution.
#' @param rdsFname User also can provide the full path of RDS file instead of 'resDir' where RDS file is saved in. If this option is used, please also provide 'resDir' to specify where the analysis results will be saved.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param expCondCheckFname suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param expCond specify the specific experimental condition for pseudo time trajectory analysis, if not specified, all experimental conditions will be performed .
#' @param cellcluster specify the specific cell cluster name for pseudo time trajectory analysis, if not specified, all cell clusters will be performed .
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @importFrom Seurat RenameIdents
#'
#' @keywords getExpCondClusterPseudotime
#' @examples getExpCondClusterPseudotime()
#' @export
#'
#' @return
#' a SingleCellExperiment where 3 methods functional pseudotime analysis results are saved in colData
## ------------------------------------------------------------------------------------ ##
getExpCondClusterPseudotime <- function(resDir=NULL, rdsFname=NULL, newAnnotation = 'F', newAnnotationRscriptName = NULL, expCondCheck='sample', expCondSepName = NULL, expCond = NULL, cellcluster = NULL){
  ##--------------------------------------------------------------------------------------##
  newAnnotation             <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Please provide corresponding 'newAnnotationRscriptName', becasue 'newAnnotation' == True.")
  ##--------------------------------------------------------------------------------------##
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
  seuratObjFinal            <<- readRDS(file = as.character(rdsFname))
  print('Done for RDS readin')
  ##--------------------------------------------------------------------------------------##
  ## update results directory if new annotation is used
  if (newAnnotation) {
    resDir <- paste(resDir, 'results_wNewAnnotation_pseudoTime', sep = '/')
  } else {
    resDir <- paste(resDir, 'results_wOrgClusterAnnotation_pseudoTime', sep = '/')
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
  resDir                <- paste(sprintf('%s/%s', resDir, expCondSepName ))
  if (!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('Pseudotime analysis results will be saved at %s', resDir))
  ##--------------------------------------------------------------------------------------##
  if (newAnnotation) {
    ## Assign cell type identity to originally identified cell clusters
    source(as.character(newAnnotationRscriptName))
  }
  print(sprintf("Current idents in 'seuratObjFinal' are as below:"))
  print(table(Seurat::Idents(seuratObjFinal)))
  print('*******************')
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
  print(sprintf('Start preparing to conduct pseudotime analysis on all identfied unsupervise/annotated cell clusters on experimental condition %s', expCondSepName))
  print(sprintf("Current experimental condition levels cell clusters partition in 'seuratObjFinal' are as below:"))
  print(table(seuratObjFinal@meta.data$expCond))
  print('*******************')
  print('=====================================================')
  ## --------------------------------------------
  ## updated expCond factor after above updating
  expCondLevels           <- levels(factor(seuratObjFinal@meta.data$expCond))
  ptGmmRes                <- list()
  if (is.null(expCond) & is.null(cellcluster)) {
    for (i in 1:length(expCondLevels)) {
      expCondLevel          <- expCondLevels[i]
      ptGmmClusterRes       <- list()
      ptGmmClusterResNames  <- c()
      ## ------
      ## subsetting different experimental conditions
      seuratObjFinalexpCond <- subset(seuratObjFinal, expCond == expCondLevel)
      ## ------
      print(sprintf('START %s: pseudotime analysis will be performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', i, expCondSepName, expCondLevel ))
      cellclusterLevels     <- levels(Idents(seuratObjFinalexpCond))
      print("pseudotime analysis will be performed on below unsupervised/annotated cell clusters:")
      print(table(Idents(seuratObjFinalexpCond)))
      print('===================')
      ## loop over each cellclusterName in 'cellclusterLevels' for pseudotime analysis
      for (c in 1:length(cellclusterLevels)) {
        tryCatch({
          cellclusterName              <- cellclusterLevels[c]
          print(sprintf("START %s.%s: pseudotime analysis with GMM clustering for cell cluster = %s in experimental condition %s: %s", i, c, cellclusterName, expCondSepName, expCondLevel ))
          seuratObjFinalexpCondCluster <- subset(seuratObjFinalexpCond, idents = as.character(cellclusterName) )
          pseudoRes                    <- NULL
          pseudoRes                    <- calcPCApseudo(obj = seuratObjFinalexpCondCluster, slingshotclusterLabels = 'GMM',
                                                        resSave = 'T', resFnamePrefix = paste(sprintf('%s/expCondLevel_%s_cluster_%s', resDir, expCondLevel, gsub('/|-|[.]','', cellclusterName))))

          if (!is.null(pseudoRes)) {
            print(sprintf("Complete %s.%s calcPCApseudo() analysis", i, c))
            print('Start to make correspdoning calcPCApseudo() plots')
            makeGMMptPlots(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, expCondLevel, gsub('/|-|[.]','', cellclusterName)))
            plotPseudotimeLineages(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, expCondLevel, gsub('/|-|[.]','', cellclusterName)))
            plotPseudotimeHeatmap(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, expCondLevel, gsub('/|-|[.]','', cellclusterName)))
            print('END to make correspdoning calcPCApseudo() plots')
            print(sprintf("FINISH %s.%s: pseudotime analysis for cell cluster = %s in experimental condition %s: %s", i, c, cellclusterName, expCondSepName, expCondLevel ))
            ptGmmClusterRes[[c]]       <- pseudoRes
            ptGmmClusterResNames       <- c(ptGmmClusterResNames, cellclusterName)
          }
          print('***********************************')
        }, error=function(e){cat(paste(sprintf("ERROR NOTE: NOT complete pseudotime analysis in cluster levels %s:", cellclusterName)),conditionMessage(e), "\n")})
      }
      names(ptGmmClusterRes)           <- ptGmmClusterResNames
      print(sprintf('END %s: FINISH pseudotime analysis with GMM clustering were performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', i, expCondSepName, expCondLevel ))
      ptGmmRes[[i]]                    <- ptGmmClusterRes
      print(sprintf("COMPLETE pseudotime analysis on experimental condtion %s", expCondLevels[i]))
      print('==========*********========*********==========*********=========')
    }
    names(ptGmmRes)                  <- expCondLevels
    ## ------------------------------------------------------------------------------------
  } else if (!is.null(expCond) & is.null(cellcluster)) {
    ## ------------------------------------------------------------------------------------
    expCondLevel          <- expCond
    ptGmmClusterRes       <- list()
    ptGmmClusterResNames  <- c()
    ## ------
    ## subsetting different experimental conditions
    seuratObjFinalexpCond <- subset(seuratObjFinal, expCond == expCondLevel)
    ## ------
    print(sprintf('START: pseudotime analysis will be performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', expCondSepName, expCondLevel ))
    cellclusterLevels     <- levels(Idents(seuratObjFinalexpCond))
    print("pseudotime analysis will be performed on below unsupervised/annotated cell clusters:")
    print(table(Idents(seuratObjFinalexpCond)))
    print('===================')
    ## loop over each cellclusterName in 'cellclusterLevels' for pseudotime analysis
    for (c in 1:length(cellclusterLevels)) {
      tryCatch({
        cellclusterName              <- cellclusterLevels[c]
        print(sprintf("START %s: pseudotime analysis with GMM clustering for cell cluster = %s in experimental condition %s: %s", c, cellclusterName, expCondSepName, expCondLevel ))
        seuratObjFinalexpCondCluster <- subset(seuratObjFinalexpCond, idents = as.character(cellclusterName) )
        pseudoRes                    <- NULL
        pseudoRes                    <- calcPCApseudo(obj = seuratObjFinalexpCondCluster, slingshotclusterLabels = 'GMM',
                                                      resSave = 'T', resFnamePrefix = paste(sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName))))

        if (!is.null(pseudoRes)) {
          print(sprintf("Complete %s calcPCApseudo() analysis", c))
          print('Start to make correspdoning calcPCApseudo() plots')
          makeGMMptPlots(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          plotPseudotimeLineages(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          plotPseudotimeHeatmap(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          print('END to make correspdoning calcPCApseudo() plots')
          print(sprintf("FINISH %s: pseudotime analysis for cell cluster = %s in experimental condition %s: %s", c, cellclusterName, expCondSepName, expCondLevel ))
          ptGmmClusterRes[[c]]       <- pseudoRes
          ptGmmClusterResNames       <- c(ptGmmClusterResNames, cellclusterName)
        }
        print('***********************************')
      }, error=function(e){cat(paste(sprintf("ERROR NOTE: NOT complete pseudotime analysis in cluster levels %s:", cellclusterName)),conditionMessage(e), "\n")})
    }
    names(ptGmmClusterRes)           <- ptGmmClusterResNames
    print(sprintf('END: FINISH pseudotime analysis with GMM clustering were performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', expCondSepName, expCondLevel ))
    ptGmmRes[[1]]                    <- ptGmmClusterRes
    print(sprintf("COMPLETE pseudotime analysis on experimental condtion %s", expCond))
    print('==========*********========*********==========*********=========')
    names(ptGmmRes)                  <- expCondLevel
    ## ------------------------------------------------------------------------------------
  } else if (is.null(expCond) & !is.null(cellcluster)) {
    ## ------------------------------------------------------------------------------------
    for (i in 1:length(expCondLevels)) {
      tryCatch({
        expCondLevel          <- expCondLevels[i]
        ptGmmClusterRes       <- list()
        ptGmmClusterResNames  <- c()
        ## ------
        ## subsetting different experimental conditions
        seuratObjFinalexpCond <- subset(seuratObjFinal, expCond == expCondLevel)
        ## ------
        print(sprintf('START %s: pseudotime analysis will be performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', i, expCondSepName, expCondLevel ))
        cellclusterLevels     <- levels(Idents(seuratObjFinalexpCond))
        print("pseudotime analysis will be performed on below unsupervised/annotated cell clusters:")
        print(table(Idents(seuratObjFinalexpCond)))
        print('===================')
        ## conduct analysis on the specified cellcluster
        cellclusterName              <- cellcluster
        print(sprintf("START %s: pseudotime analysis with GMM clustering for cell cluster = %s in experimental condition %s: %s", i, cellclusterName, expCondSepName, expCondLevel ))
        seuratObjFinalexpCondCluster <- subset(seuratObjFinalexpCond, idents = as.character(cellclusterName) )
        pseudoRes                    <- NULL
        pseudoRes                    <- calcPCApseudo(obj = seuratObjFinalexpCondCluster, slingshotclusterLabels = 'GMM',
                                                      resSave = 'T', resFnamePrefix = paste(sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName) )))

        if (!is.null(pseudoRes)) {
          print(sprintf("Complete %s calcPCApseudo() analysis", i))
          print('Start to make correspdoning calcPCApseudo() plots')
          makeGMMptPlots(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          plotPseudotimeLineages(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          plotPseudotimeHeatmap(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          print('END to make correspdoning calcPCApseudo() plots')
          print(sprintf("FINISH %s: pseudotime analysis for cell cluster = %s in experimental condition %s: %s", i, cellclusterName, expCondSepName, expCondLevel ))
          ptGmmClusterRes[[1]]       <- pseudoRes
          print('***********************************')
          names(ptGmmClusterRes)     <- cellclusterName
          print(sprintf('END %s: FINISH pseudotime analysis with GMM clustering were performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', i, expCondSepName, expCondLevel ))
        }
        ptGmmRes[[i]]                <- ptGmmClusterRes
        print(sprintf("COMPLETE pseudotime analysis on experimental condtion %s", expCondLevel))
        print('==========*********========*********==========*********=========')
      }, error=function(e){cat(paste(sprintf("ERROR NOTE: NOT complete pseudotime analysis in expCond %s cluster levels %s :", expCondLevel, cellclusterName)),conditionMessage(e), "\n")})
    }
    names(ptGmmRes)                  <- expCondLevels
    ## ------------------------------------------------------------------------------------
  } else if (!is.null(expCond) & !is.null(cellcluster)) {
    ## ------------------------------------------------------------------------------------
    expCondLevel          <- expCond
    ptGmmClusterRes       <- list()
    ptGmmClusterResNames  <- c()
    ## ------
    ## subsetting different experimental conditions
    seuratObjFinalexpCond <- subset(seuratObjFinal, expCond == expCondLevel)
    ## ------
    print(sprintf('START: pseudotime analysis will be performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', expCondSepName, expCondLevel ))
    cellclusterLevels     <- levels(Idents(seuratObjFinalexpCond))
    print("pseudotime analysis will be performed on below unsupervised/annotated cell clusters:")
    print(table(Idents(seuratObjFinalexpCond)))
    print('===================')
    ## conduct analysis on the specified cellcluster
    cellclusterName              <- cellcluster
    print(sprintf("START pseudotime analysis with GMM clustering for cell cluster = %s in experimental condition %s: %s", cellclusterName, expCondSepName, expCondLevel ))
    seuratObjFinalexpCondCluster <- subset(seuratObjFinalexpCond, idents = as.character(cellclusterName) )
    pseudoRes                    <- NULL
    tryCatch({
      pseudoRes                  <- calcPCApseudo(obj = seuratObjFinalexpCondCluster, slingshotclusterLabels = 'GMM',
                                                  resSave = 'T', resFnamePrefix = paste(sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName))))
    }, error=function(e){cat(paste(sprintf("ERROR NOTE: NOT complete pseudotime analysis in expCond %s cluster levels %s :", expCondLevel, cellclusterName)),conditionMessage(e), "\n")})

    if (!is.null(pseudoRes)) {
      print(sprintf("Complete calcPCApseudo() analysis"))
      print('Start to make correspdoning calcPCApseudo() plots')
      makeGMMptPlots(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
      plotPseudotimeLineages(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
      plotPseudotimeHeatmap(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
      print('END to make correspdoning calcPCApseudo() plots')
      print(sprintf("FINISH pseudotime analysis for cell cluster = %s in experimental condition %s: %s", cellclusterName, expCondSepName, expCondLevel ))
      ptGmmClusterRes[[1]]       <- pseudoRes
      names(ptGmmClusterRes)     <- cellclusterName
      print(sprintf('END: FINISH pseudotime analysis with GMM clustering were performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', expCondSepName, expCondLevel ))
      ptGmmRes[[1]]              <- ptGmmClusterRes
      print(sprintf("COMPLETE pseudotime analysis on experimental condtion %s", expCondLevel))
      print('==========*********========*********==========*********=========')
      names(ptGmmRes)            <- expCondLevel
    }
    ## ------------------------------------------------------------------------------------
  }
  return(ptGmmRes)
}
## ------------------------------------------------------------------------------------ ##
