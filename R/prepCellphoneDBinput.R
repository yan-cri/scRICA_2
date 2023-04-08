## prepCellphoneDBinput(): prepare the cellphoneDB input files from the normalized count ##
## Developed by Yan Li
##-------------------------------------------------------------------------------------- ##
#' getExpCondClusterPseudotime() Function
#' @details
#' This function is used to perform functional pseudotime analysis via PCA, Diffusion Map, and slingshot
#' @param rds User also can provide the full path of RDS file instead of 'resDir' where RDS file is saved in. If this option is used, please also provide 'resDir' to specify where the analysis results will be saved.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param resDir full path of folder directory where results to be saved.
#' @param resFnamePrefix result file name prefix.
#' @param expCondCheck specify which experimental conditions to be explored, including sample or expCond1/2/... presented in the input RDS metadata table column names.
#' @param expCond specify the specific experimental condition for pseudo time trajectory analysis, if not specified, all experimental conditions will be performed.
#' @param cellcluster specify the specific cell cluster name for pseudo time trajectory analysis, if not specified, all cell clusters will be performed.
#' @param debug whether to turn on for debug check, by default FALSE (turned off).
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
prepCellphoneDBinput <- function(rds, newAnnotation = 'F', newAnnotationRscriptName = NULL,
                                 resDir = NULL, resFnamePrefix = 'testInput', expCondCheck='sample', expCond = NULL, cellcluster = NULL, debug = as.logical(F)){
  ##--------------------------------------------------------------------------------------##
  newAnnotation             <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Please provide corresponding 'newAnnotationRscriptName', becasue 'newAnnotation' == True.")
  ##--------------------------------------------------------------------------------------##
  ## stop program if input options do not match
  if (is.list(rds)) {
    if(newAnnotation) {
      if(length(newAnnotationRscriptName) != length(rds)) stop("Error: please provide corresponding 'newAnnotationRscriptName' files correctly.")
    }
    if (!is.null(expCond)) {
      if(length(expCond) != length(rds)) stop("Error: 'expCond' has different length from input 'rds' list, please provide corresponding 'expCond' list option correctly.")
      if (!is.list(expCond)) stop("Error: 'expCond' is not list item, please provide corresponding 'expCond' list option correctly.")
    }
    if (!is.null(cellcluster)) {
      if(length(cellcluster) != length(rds)) stop("Error: 'cellcluster' has different length from input 'rds' list, please provide corresponding 'cellcluster' list option correctly.")
      if (!is.list(cellcluster)) stop("Error: 'cellcluster' is not list item, please provide corresponding 'cellcluster' list option correctly.")
    }
    if(length(expCondCheck) != length(rds)) stop("Error: 'expCondCheck' has different length from input 'rds' list, please provide corresponding 'expCondCheck' list option correctly.")
    if (!is.list(expCondCheck)) stop("Error: 'expCondCheck' is not list item, please provide corresponding 'expCondCheck' list option correctly.")
  }
  ##--------------------------------------------------------------------------------------##
  ## creat resDir for results to save
  if (is.null(resDir)) {
    resDir                <- getwd()
  } else {
    if (!dir.exists(resDir)) dir.create(resDir, mode = '0700')
  }
  # ##--------------------------------------------------------------------------------------##
  # if (is.list(rds)) {
  #   if (length(expCondCheck)==1) {
  #     if (expCondCheck == 'sample') {
  #       expCondCheckFname = 'expCond_sample'
  #     } else {
  #       expCondCheckFname   = paste('cond', expCondCheck, sep = '_')
  #     }
  #   } else {
  #     expCondCheckFname = paste(paste(paste('cond',1:2, sep = ''),unlist(expCondCheck), sep = ''), collapse = '_')
  #   }
  # }
  # ##----------------------
  # ## update resDir based on 'expCondCheck'
  resDir                <- paste(sprintf('%s/%s_cellphoneDB_inputFiles', resDir, resFnamePrefix ))
  if (!dir.exists(resDir)) dir.create(resDir, mode = '0700')
  print(sprintf('cellphoneDB analysis input files will be saved at %s', resDir))
  ##--------------------------------------------------------------------------------------##
  ##--------------------------------------------------------------------------------------##
  seuratObjFinalList <- list()
  if (is.list(rds)) {
    for (r in 1:length(rds)) {
      if (class(rds[[r]]) == 'Seurat') {
        print(sprintf("Read in rds %s", r))
        seuratObjFinal      <<- rds[[r]]
        print('RDS is provided with rds option')
        if (newAnnotation) source(newAnnotationRscriptName[[r]])
        print('-=-=-=-=-')
      } else {
        print(sprintf("Read in rds %s", r))
        rdsFname            <- rds[[r]]
        if (!file.exists(rdsFname)) stop("Please provide correct rds file name in 'rds'.")
        seuratObjFinal      <<- readRDS(file = as.character(rdsFname))
        print('Done for RDS read in')
        print('-=-=-=-=-')
      }
      seuratObjFinalList[[r]] <- seuratObjFinal
    }
  } else {
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
    if (newAnnotation) source(newAnnotationRscriptName)
    seuratObjFinalList[[1]] <- seuratObjFinal
  }
  if (debug) print(length(seuratObjFinalList))
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond'
  for (i in 1:length(seuratObjFinalList)) {
    seuratObjFinal <- seuratObjFinalList[[i]]
    if (length(expCondCheck)==1) {
      expCondCheck = expCondCheck
    } else {
      expCondCheck = expCondCheck[[i]]
    }
    ## -
    if (expCondCheck == 'sample' | expCondCheck == 'ALL'| expCondCheck == 'all' | expCondCheck == 'All') {
      seuratObjFinal                     <- seuratObjFinal
    } else {
      if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
        stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
      } else {
        seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
      }
    }
    ## -
    if (expCondCheck == 'sample' | expCondCheck == 'ALL'| expCondCheck == 'all' | expCondCheck == 'All') {
      print(sprintf("Current experimental condition cell clusters partition in rds.%s are on idents only as below:", i))
      print(table(Seurat::Idents(seuratObjFinal)))
      print('*******************')
    } else {
      print(sprintf("Current experimental condition cell clusters partition in rds.%s are as below:", i))
      print(table(seuratObjFinal@meta.data$expCond))
      print('*******************')
    }
    seuratObjFinalList[[i]] <- seuratObjFinal
  }
  ##--------------------------------------------------------------------------------------##
  ## subset cell clusters if cellcluster is not NULL or all|All|ALL
  if (!is.null(cellcluster)) {
    if (is.list(rds)) {
      for (r in 1:length(seuratObjFinalList)) {
        if (length(cellcluster[[r]])==1) {
          if (grepl('all|All|ALL', cellcluster[[r]])) {
            seuratObjFinalList[[r]] <-  seuratObjFinalList[[r]]
          } else {
            seuratObjFinalList[[r]] <-  subset(seuratObjFinalList[[r]], idents = as.character(cellcluster[[r]]) )
          }
        } else {
          if (any(!cellcluster[[r]] %in% levels(Idents(seuratObjFinalList[[r]])) ) ) stop("Error: subsetting 'cellcluster' does not exist in the rds input")
          seuratObjFinalList[[r]]   <-  subset(seuratObjFinalList[[r]], idents = as.character(cellcluster[[r]]) )
        }
      }
    } else {
      if (length(cellcluster) ==1) {
        if (grepl('all|All|ALL', cellcluster)) {
          seuratObjFinalList[[1]]   <-  seuratObjFinalList[[1]]
        } else {
          if (any(!cellcluster %in% levels(Idents(seuratObjFinalList[[1]])) ) ) stop("Error: subsetting 'cellcluster' does not exist in the rds input")
          seuratObjFinalList[[1]]   <-  subset(seuratObjFinalList[[1]], idents = as.character(cellcluster) )
        }
      } else {
        if (any(!cellcluster %in% levels(Idents(seuratObjFinalList[[1]])) ) ) stop("Error: subsetting 'cellcluster' does not exist in the rds input")
        seuratObjFinalList[[1]]     <-  subset(seuratObjFinalList[[1]], idents = as.character(cellcluster) )
      }
    }
  }
  ##--------------------------------------------------------------------------------------##
  ##subsetting expCond to prepare cellphoneDB input
  for ( e in 1:length(expCond) ) {
    expConds1             <- expCond[[e]]
    if (length(expConds1)==1) {
      if (expConds1 == 'org' | expConds1 == 'All' | expConds1 == 'all' | expConds1 =='ALL') {
        seuratObj4cellphonedb1 = seuratObjFinalList[[e]]
        # seuratObj4cellphonedb1@meta.data$expCond <- gsub(pattern = as.character(expCondName2change[[e]]), replacement = '', x = seuratObj4cellphonedb1@meta.data$expCond)
      } else{
        seuratObj4cellphonedb1 <- subset(seuratObjFinalList[[e]],  subset = expCond == as.character(expConds1))
      }
      countMatrix1Gname        <- as.data.frame(seuratObj4cellphonedb1@assays$RNA@data)
      countMatrix1Gname$Gene   <- rownames(countMatrix1Gname)
      countMatrix1Gname2       <- countMatrix1Gname %>% select(Gene, everything())
      ## -
      cellMeta1                <- seuratObj4cellphonedb1@meta.data %>% select(expCond) ##%>% rename(expCond = expCond)
      cellMeta1$idents         <- as.factor(Idents(seuratObj4cellphonedb1))
      if (expCondCheck=='sample'|expCondCheck=='all'|expCondCheck=='ALL'|expCondCheck=='All') {
        cellMeta1$cellType       <- cellMeta1$idents
      } else {
        cellMeta1$cellType       <- paste(cellMeta1$expCond, cellMeta1$idents, sep = '_')
      }
      cellMeta12               <- cellMeta1 %>% select(cellType)
      cellMeta12               <- tibble::rownames_to_column(cellMeta12, "Cell")
      print("cellphoneDB interactions analysis on:")
      print(table(cellMeta12$cellType))
      print('***********************')
      ## -
    } else {
      ## ---
      seuratObj4cellphonedb1   = seuratObjFinalList[[e]]
      for(i in 1: length(expConds1) ) {
        seuratObjFinalPrep    <- subset(seuratObj4cellphonedb1,  subset = expCond == as.character(expConds1[i]))
        ## -
        countMatrix1PrepGname  <- as.data.frame(seuratObjFinalPrep@assays$RNA@data)
        countMatrix1PrepGname$Gene <- rownames(countMatrix1PrepGname)
        countMatrix1PrepGname2 <- countMatrix1PrepGname %>% select(Gene, everything())
        ## -
        cellMeta1Prep          <- seuratObjFinalPrep@meta.data %>% select(expCond) %>% rename(expCond = expCond)
        cellMeta1Prep$idents   <- as.factor(Idents(seuratObjFinalPrep))
        cellMeta1Prep$cellType <- paste(cellMeta1Prep$expCond, cellMeta1Prep$idents, sep = '_')
        cellMeta1Prep2         <- cellMeta1Prep %>% select(cellType)
        cellMeta1Prep2         <- tibble::rownames_to_column(cellMeta1Prep2, "Cell")
        ## -
        ## -
        if (i ==1) {
          seuratObjFinalPrep2  = seuratObjFinalPrep
          countMatrix1Gname2   = countMatrix1PrepGname2
          cellMeta12           = cellMeta1Prep2
        } else {
          seuratObjFinalPrep2 <- merge(seuratObjFinalPrep2, seuratObjFinalPrep)
          if (all(countMatrix1Gname2$Gene == countMatrix1PrepGname2$Gene)){
            countMatrix1Gname2   = cbind(countMatrix1Gname2, (countMatrix1PrepGname2%>%select(-Gene)))
          } else {
            print('Gene name from scale count slot do not match, join are used here')
            countMatrix1Gname2   = data.table::merge.data.table(countMatrix1Gname2, countMatrix1PrepGname2, by = 'Gene')
          }
          ## -
          cellMeta12           = rbind(cellMeta12, cellMeta1Prep2)
        }
      }
      seuratObj4cellphonedb1  = seuratObjFinalPrep2
      # if (!all(colnames(countMatrix1Gname2%>%select(-Gene)) == cellMeta12$Cell)) stop('extract cellMeta12 & countMatrix1Gname2 for cellphoneDB do not match')
      print(sprintf("cellphoneDB interactions analysis on rds%s:", e))
      print(table(cellMeta12$cellType))
      print('***********************')
      ## ---
    }
    # cellMetaResList[[e]]     <- cellMeta12
    # countMatrixResList[[e]]  <- countMatrix1Gname2
    if (e == 1) {
      cellMetaRes    <- cellMeta12
      countMatrixRes <- countMatrix1Gname2
    } else {
      cellMetaRes    <- rbind(cellMetaRes, cellMeta12)
      countMatrixRes <- data.table::merge.data.table(countMatrixRes, countMatrix1Gname2, by = 'Gene')
    }
  }
  ## update resDir to save input files
  print('writting cellphonDB input tables')
  countFname        <- file.path(resDir, sprintf('cellphoneDb_%s_Count.txt', resFnamePrefix) )
  metaFname         <- file.path(resDir, sprintf('cellphoneDb_%s_Meta.txt', resFnamePrefix) )
  ## ---
  if (sum(colnames(countMatrixRes)[-1] == cellMetaRes$Cell)!=length(cellMetaRes$Cell)) stop('countmatrix does not mathch count meta table')
  print(sprintf("count table size, %s, %s", dim(countMatrixRes)[1], dim(countMatrixRes)[2]))
  print(sprintf("cell table size, %s, %s", dim(cellMetaRes)[1], dim(cellMetaRes)[2]))
  print('writting count table')
  write.table(x = countMatrixRes, file = countFname, quote = F, sep = '\t', row.names = F, col.names = T)
  print('writting meta table')
  write.table(x = cellMetaRes, file = metaFname, quote = F, sep = '\t', row.names = F, col.names = T)
  print(sprintf("%s cells will be processed", dim(cellMeta12)[1]))
  print("END-=-=-=-END-=-=-=-END")
  ## ------------------------------------------------------------------------------------ ##
  return(list(resDir = resDir, countFname = countFname, metaFname = metaFname))
}
## ------------------------------------------------------------------------------------ ##
