## getExpCondClusterPseudotime(): perform functional pseudotime analysis on cell clusters specified from different experimental condition types##
## Developed by Yan Li, July, 2021
##--------------------------------------------------------------------------------------##
#' getExpCondClusterPseudotime() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param resDir specify an exiting full path of directory, where results will be saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved with above rdsDir option.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param expCondCheckFname suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param expCond specify a specific experimental condition for pseudo time trajectory analysis, if not specified, all experimental conditions will be performed .
#' @param cellcluster specify the specific cell cluster name for pseudo time trajectory analysis, if not specified, all cell clusters will be performed.
#' @param slingshotclusterLabels specify slighshot cluster label options, 3 options are available here, they are 'NULL', 'seurat_clusters', or 'GMM'. By default is 'GMM'.
#' @param topFeatureNo specify the number of features used in PCA pseudo-time trajectory calculation.
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
getExpCondClusterPseudotime <- function(resDir=NULL, rds=NULL, newAnnotation = 'F',
                                        newAnnotationRscriptName = NULL,
                                        expCondCheck='sample', expCondCheckFname = NULL,
                                        expCond = NULL, cellcluster = NULL,
                                        slingshotclusterLabels = 'GMM', topFeatureNo = 2000){
  ##--------------------------------------------------------------------------------------##
  newAnnotation             <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Please provide corresponding 'newAnnotationRscriptName', becasue 'newAnnotation' == True.")
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
  resDir                <- paste(sprintf('%s/%s', resDir, expCondCheckFname ))
  if (!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('Pseudotime analysis results will be saved at %s', resDir))
  ##--------------------------------------------------------------------------------------##
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    print("Before adding new annotation")
    print(table(Idents(seuratObjFinal)))
    source(newAnnotationRscriptName)
  }
  print(sprintf("Current idents in 'seuratObjFinal' are as below:"))
  print(table(Seurat::Idents(seuratObjFinal)))
  print('*******************')
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond'
  if (expCondCheck == 'sample') {
    seuratObjFinal                     <- seuratObjFinal
  } else {
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  ##--------------------------------------------------------------------------------------##
  print(sprintf('Start preparing to conduct pseudotime analysis on all identfied unsupervise/annotated cell clusters on experimental condition %s', expCondCheckFname))
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
      print(sprintf('START %s: pseudotime analysis will be performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', i, expCondCheck, expCondLevel ))
      cellclusterLevels     <- levels(Idents(seuratObjFinalexpCond))
      print("pseudotime analysis will be performed on below unsupervised/annotated cell clusters:")
      print(table(Idents(seuratObjFinalexpCond)))
      print('===================')
      ## loop over each cellclusterName in 'cellclusterLevels' for pseudotime analysis
      for (c in 1:length(cellclusterLevels)) {
        tryCatch({
          cellclusterName              <- cellclusterLevels[c]
          print(sprintf("START %s.%s: pseudotime analysis with GMM clustering for cell cluster = %s in experimental condition %s: %s", i, c, cellclusterName, expCondCheck, expCondLevel ))
          seuratObjFinalexpCondCluster <- subset(seuratObjFinalexpCond, idents = as.character(cellclusterName) )
          pseudoRes                    <- NULL
          pseudoRes                    <- calcPCApseudo(obj = seuratObjFinalexpCondCluster, slingshotclusterLabels = slingshotclusterLabels, topFeatureNo = topFeatureNo,
                                                        resSave = 'T', resFnamePrefix = paste(sprintf('%s/expCondLevel_%s_cluster_%s', resDir, expCondLevel, gsub('/|-|[.]','', cellclusterName))))

          if (!is.null(pseudoRes)) {
            print(sprintf("Complete %s.%s calcPCApseudo() analysis", i, c))
            print('Start to make correspdoning calcPCApseudo() plots')
            makeGMMptPlots(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, expCondLevel, gsub('/|-|[.]','', cellclusterName)))
            plotPseudotimeLineages(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, expCondLevel, gsub('/|-|[.]','', cellclusterName)))
            plotPseudotimeHeatmap(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, expCondLevel, gsub('/|-|[.]','', cellclusterName)))
            print('END to make correspdoning calcPCApseudo() plots')
            print(sprintf("FINISH %s.%s: pseudotime analysis for cell cluster = %s in experimental condition %s: %s", i, c, cellclusterName, expCondCheck, expCondLevel ))
            ptGmmClusterRes[[c]]       <- pseudoRes
            ptGmmClusterResNames       <- c(ptGmmClusterResNames, cellclusterName)
          }
          print('***********************************')
        }, error=function(e){cat(paste(sprintf("ERROR NOTE: NOT complete pseudotime analysis in cluster levels %s:", cellclusterName)),conditionMessage(e), "\n")})
      }
      names(ptGmmClusterRes)           <- ptGmmClusterResNames
      print(sprintf('END %s: FINISH pseudotime analysis with GMM clustering were performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', i, expCondCheck, expCondLevel ))
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
    print(sprintf('START: pseudotime analysis will be performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', expCondCheck, expCondLevel ))
    cellclusterLevels     <- levels(Idents(seuratObjFinalexpCond))
    print("pseudotime analysis will be performed on below unsupervised/annotated cell clusters:")
    print(table(Idents(seuratObjFinalexpCond)))
    print('===================')
    ## loop over each cellclusterName in 'cellclusterLevels' for pseudotime analysis
    for (c in 1:length(cellclusterLevels)) {
      tryCatch({
        cellclusterName              <- cellclusterLevels[c]
        print(sprintf("START %s: pseudotime analysis with GMM clustering for cell cluster = %s in experimental condition %s: %s", c, cellclusterName, expCondCheck, expCondLevel ))
        seuratObjFinalexpCondCluster <- subset(seuratObjFinalexpCond, idents = as.character(cellclusterName) )
        pseudoRes                    <- NULL
        pseudoRes                    <- calcPCApseudo(obj = seuratObjFinalexpCondCluster, slingshotclusterLabels = slingshotclusterLabels, topFeatureNo = topFeatureNo,
                                                      resSave = 'T', resFnamePrefix = paste(sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName))))

        if (!is.null(pseudoRes)) {
          print(sprintf("Complete %s calcPCApseudo() analysis", c))
          print('Start to make correspdoning calcPCApseudo() plots')
          makeGMMptPlots(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          plotPseudotimeLineages(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          plotPseudotimeHeatmap(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          print('END to make correspdoning calcPCApseudo() plots')
          print(sprintf("FINISH %s: pseudotime analysis for cell cluster = %s in experimental condition %s: %s", c, cellclusterName, expCondCheck, expCondLevel ))
          ptGmmClusterRes[[c]]       <- pseudoRes
          ptGmmClusterResNames       <- c(ptGmmClusterResNames, cellclusterName)
        }
        print('***********************************')
      }, error=function(e){cat(paste(sprintf("ERROR NOTE: NOT complete pseudotime analysis in cluster levels %s:", cellclusterName)),conditionMessage(e), "\n")})
    }
    names(ptGmmClusterRes)           <- ptGmmClusterResNames
    print(sprintf('END: FINISH pseudotime analysis with GMM clustering were performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', expCondCheck, expCondLevel ))
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

        # print("pseudotime analysis will be performed on below unsupervised/annotated cell clusters:")
        # print(table(Idents(seuratObjFinalexpCond)))
        # print('===================')
        ## conduct analysis on the specified cellcluster
        cellclusterName              <- cellcluster
        print(sprintf("START %s: pseudotime analysis will be performed on experimental condition %s level: %s for '%s' cell clusters", i, expCondCheck, expCondLevel, cellcluster ))
        cellclusterLevels     <- levels(Idents(seuratObjFinalexpCond))
        if (any(!cellcluster %in% cellclusterLevels ) ) stop('Please provide the corresponding cell clusters in identfied idents for pseudotime analysis.')
        # print(sprintf("START %s: pseudotime analysis with GMM clustering for cell cluster = %s in experimental condition %s: %s", i, cellclusterName, expCondCheck, expCondLevel ))
        seuratObjFinalexpCondCluster <- subset(seuratObjFinalexpCond, idents = as.character(cellclusterName) )
        pseudoRes                    <- calcPCApseudo(obj = seuratObjFinalexpCondCluster, slingshotclusterLabels = slingshotclusterLabels, topFeatureNo = topFeatureNo,
                                                      resSave = 'T', resFnamePrefix = paste(sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName) )) )

        if (!is.null(pseudoRes)) {
          print(sprintf("Complete %s calcPCApseudo() analysis", i))
          print('Start to make correspdoning calcPCApseudo() plots')
          makeGMMptPlots(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          plotPseudotimeLineages(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          plotPseudotimeHeatmap(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
          print('END to make correspdoning calcPCApseudo() plots')
          print(sprintf("FINISH %s: pseudotime analysis for cell cluster = %s in experimental condition %s: %s", i, cellclusterName, expCondCheck, expCondLevel ))
          ptGmmClusterRes[[1]]       <- pseudoRes
          print('***********************************')
          names(ptGmmClusterRes)     <- cellclusterName
          print(sprintf('END %s: FINISH pseudotime analysis with GMM clustering were performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', i, expCondCheck, expCondLevel ))
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
    print(sprintf('START: pseudotime analysis will be performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', expCondCheck, expCondLevel ))
    cellclusterLevels     <- levels(Idents(seuratObjFinalexpCond))
    print("pseudotime analysis will be performed on below unsupervised/annotated cell clusters:")
    print(table(Idents(seuratObjFinalexpCond)))
    print('===================')
    ## conduct analysis on the specified cellcluster
    cellclusterName              <- cellcluster
    print(sprintf("START pseudotime analysis with GMM clustering for cell cluster = %s in experimental condition %s: %s", cellclusterName, expCondCheck, expCondLevel ))
    seuratObjFinalexpCondCluster <- subset(seuratObjFinalexpCond, idents = as.character(cellclusterName) )
    pseudoRes                    <- NULL
    tryCatch({
      pseudoRes                  <- calcPCApseudo(obj = seuratObjFinalexpCondCluster, slingshotclusterLabels = slingshotclusterLabels, topFeatureNo = topFeatureNo,
                                                  resSave = 'T', resFnamePrefix = paste(sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName))))
    }, error=function(e){cat(paste(sprintf("ERROR NOTE: NOT complete pseudotime analysis in expCond %s cluster levels %s :", expCondLevel, cellclusterName)),conditionMessage(e), "\n")})

    if (!is.null(pseudoRes)) {
      print(sprintf("Complete calcPCApseudo() analysis"))
      print('Start to make correspdoning calcPCApseudo() plots')
      makeGMMptPlots(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
      plotPseudotimeLineages(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
      plotPseudotimeHeatmap(pseudoRes = pseudoRes, plotname = sprintf('%s/expCondLevel_%s_cluster_%s', resDir, gsub('/|-|[.]','', expCondLevel), gsub('/|-|[.]','', cellclusterName)))
      print('END to make correspdoning calcPCApseudo() plots')
      print(sprintf("FINISH pseudotime analysis for cell cluster = %s in experimental condition %s: %s", cellclusterName, expCondCheck, expCondLevel ))
      ptGmmClusterRes[[1]]       <- pseudoRes
      names(ptGmmClusterRes)     <- cellclusterName
      print(sprintf('END: FINISH pseudotime analysis with GMM clustering were performed on experimental condition %s level: %s for all seurat identified unsupervised/annotated cell clusters', expCondCheck, expCondLevel ))
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
