## processQC(): readin cellranger results into seurat object w/wo doublets removal##
## Developed by Yan Li, Jan, 2021
##--------------------------------------------------------------------------------------##
#' processQC() Function
#' @details
#' This function is used to readin cellranger results into seurat object by defining w/wo doublets removal
#' @param metadata txt file with 4 columns,
#'                 where sample specifiy the sample names,
#'                 path specify full path to cellranger analysis results for that sample,
#'                 doubletsRmMethods specify doublet removal methods,
#'                 doubletsResDir specify the full path of saved doublet removal results.
#' @param resDirName optional, define folder/directory name where integration analysis results will be saved
#' @param genomeSpecies To be added
#' @param minCells optional, default = 3, used in creating seurat object to indicate the minumn detected number of cells to be included.
#' @param minFeatures optional, default = 200, used in creating seurat object to indicate the minimum detected number of gene features (genes) of included cells.
#' @param mtFiltering default = F, logical, whehter to filter mitochondrial content
#' @param mtPerCutoff optional, if mtFiltering = T, default = 20, indicating the percentage cut-off of mitochondrial content
#'
#' @importFrom Seurat Read10X
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat AddMetaData
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom Seurat FeatureScatter
#' @importFrom Seurat LabelPoints
#' @importFrom Seurat VlnPlot
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat VariableFeaturePlot
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom stats filter
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom utils write.table
#'
#' @keywords processQC
#' @examples processQC()
#' @export
#'
#' @return
#' a list item including 3 elements:
#' 1. 'countReadInOjb': a list of orignial counts seurat objects
#' 2. 'qcProcessObj': a list of QC passed seurat objects, such as mitochondrial content removal
#' 3. 'resDir': full path of directory/folder name where the entire QC analysis results are saved, it includes
#' TO BE ADDED.
##----------------------------------------------------------------------------------------##
# library(Seurat)
# library(ggplot2)
# library(dplyr)
##----------------------------------------------------------------------------------------
## A total of 9 options for this function if 'rmDoublets' is on, otherwise, only first 5 options need to be defined.
## This function operation will automately create 'resDir' under current working directory (getwd()) and corresponding 'rdsDir' inside 'resDir'.
## Due to time consuming data readin/out to garder server, temperarely 'rdsDir' is redifined locally.
## returning 4 objects: created 'resDir', 'rdsDir', executed 'seuratObjList' (original readin from 10x w/wo doublet removal), and 'seuratQcProcessObjList' (after QC process)
## 1. 'Org_cellNoSummary.txt': summarize No. of cells from each item of provided 'cellrangerResList'
## 2. 'MT_percentage_summary_' + names(cellrangerResList)[x].txt: for each items in 'cellrangerResList' to summarize MT percentage.
## 3. if 'mtFiltering' is on: 'No_mtFiltered_cells_summary_'+ names(cellrangerResList)[x].txt: summarize the number of cells after MT and low RNA feature cells removal
## 4. 'featureScatter_'+ names(seuratObjList)[x].pdf: featureScatter plot for each items in 'cellrangerResList'
##    if 'mtFiltering' is on: 'featureScatter_'+ names(seuratObjList)[x] +'_afterFiltering.pdf'
## 5. 'featureViolin_', names(seuratObjList)[x].pdf: featureViolin plot for each items in 'cellrangerResList'
##    if 'mtFiltering' is on: 'featureViolin_', names(seuratObjList)[x] +'_afterFiltering.pdf'
## 6. 'topVariableFeature_', names(seuratObjList)[x], '.pdf': topVariableFeature plot for each items in 'cellrangerResList'
processQC <- function(metadata, resDirName=NULL, genomeSpecies=NULL, minCells=3, minFeatures=200, mtFiltering=F, mtPerCutoff=NULL) {
  ## ---
  if (!all(c("sample", "path", 'expCond1',"doubletsRmMethod" ) %in% colnames(metadata))) stop('Please provide metadata table with at least 4 columns: sample, path, expCond1, and doubletsRmMethod')
  ## ---
  if (is.null(resDirName)) resDirName <- as.character('scRICA_results')
  # if (!is.list(cellrangerResList)) stop("Please provide list item of 'cellrangerResList', which is required.")
  minCells                       <- as.numeric(minCells)
  minFeatures                    <- as.numeric(minFeatures)
  mtFiltering                    <- as.logical(mtFiltering)
  if (is.null(genomeSpecies)) genomeSpecies <- as.character('human')
  if (mtFiltering & is.null(mtPerCutoff)) stop("Mitochondrial content filtering option ('mtFiltering') is on, please provide corresponding Mitochondrial content percentage filtering option in 'mtPerCutoff'.  ")
  mtPerCutoff                    <- as.numeric(mtPerCutoff)
  ## ---
  if ( 'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata)) {
    print(sprintf('Only 1 experimental condition are provided with %s experimental factor levels for comparisions.', length(levels(factor(metadata$expCond1))) ))
  } else if ( 'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
    print(sprintf('2 experimental conditions are provided, experimental condition 1 has %s experimental factor levels, and experimental condition 2 has %s experimental factor levels for comparisions', length(levels(factor(metadata$expCond1))), length(levels(factor(metadata$expCond2))) ))
  } else if ( !'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata) ) {
    print(sprintf('No experimental condition factors are provided, the analysis will be conducted only based samples integration'))
  } else if ( !'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata)) {
    print(sprintf('Only 1 experimental condition are provided with %s experimental factor levels for comparisions.', length(levels(factor(metadata$expCond2))) ))
  }
  ## ---
  ## updating metadata to include doublets detection results
  if ( !'doubletsResDir' %in% colnames(metadata) ) {
    metadataOrg                  <- metadata
    metadataOrg$doubletsRmMethod <- gsub('[]|[ ]', 'none', metadataOrg$doubletsRmMethod)
    metadataOrg$doubletsRmMethod[is.na(metadataOrg$doubletsRmMethod)] <- 'none'
    if(all(metadataOrg$doubletsRmMethod=='None'|metadataOrg$doubletsRmMethod=='none'|metadataOrg$doubletsRmMethod=='NONE')) {
      metadata                   <- metadataOrg
      metadata$doubletsResDir    <- as.character('NA')
    } else {
      md4doublet                 <- metadataOrg %>% dplyr::filter(doubletsRmMethod != 'None' ) %>% dplyr::filter(doubletsRmMethod != 'none' ) %>% dplyr::filter(doubletsRmMethod != 'NONE' )
      print('===================================================================')
      print('START: Doublets identfication analysis before processing to the next step.')
      doubletsRes                <- findDoublets(metadata = md4doublet, genomeSpecies = genomeSpecies, resFilename = paste(resDirName, 'doublet_results', sep = '/') )
      print('END: Doublets identfication, process to next step QC.')
      print('===================================================================')
      md4doublet$doubletsResDir  <- rep(as.character(doubletsRes), length(md4doublet$sample))
      doubletsMd                 <- md4doublet %>% dplyr::select(c('sample', 'doubletsResDir'))
      metadata                   <- dplyr::left_join(x = metadataOrg, y = doubletsMd, by = 'sample', )
    }
  } else {
    print('===================================================================')
    print('NO doublets identfication analysis is needed, because doublets identification either conducted or not needed.')
    metadataOrg                  <- metadata
    metadataOrg$doubletsRmMethod <- gsub('[]|[ ]', 'none', metadataOrg$doubletsRmMethod)
    metadataOrg$doubletsRmMethod[is.na(metadataOrg$doubletsRmMethod)] <- 'none'
    metadataOrg$doubletsRmMethod[metadataOrg$doubletsRmMethod=='NONE'] <- 'none'
    metadataOrg$doubletsRmMethod[metadataOrg$doubletsRmMethod=='None'] <- 'none'
    metadataOrgDbNotNa           <- metadataOrg %>% dplyr::filter(doubletsRmMethod != 'None') %>% dplyr::filter(doubletsRmMethod != 'NONE') %>% dplyr::filter(doubletsRmMethod != 'none')
    if(dim(metadataOrgDbNotNa)[1]>0) {
      if (!all(dir.exists(metadataOrgDbNotNa$doubletsResDir))) stop("please provide correct corresponding 'doubletsResDir' in 'metadata' for column 'doubletsRmMethod' specified NOT as 'none'.")
    }
    metadata                     <- metadataOrg
    print('Process to QC.')
    print('===================================================================')
  }
  # print(metadata)
  # print('*****************************')
  ## ---
  cellrangerResList              <- meata2list(metadata = metadata) ## by default names(cellrangerResList) = metadata$sample
  ## ---
  doubletsRmMethods              <- as.character(metadata$doubletsRmMethod)
  doubletsResDirs                <- as.character(metadata$doubletsResDir)
  ## 0. create 'resDir' based on provided 'resDirName' under current workDir
  resDir                         <- paste(getwd(), resDirName, sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  # ## intermediate RDS result dir
  # rdsDir               <- paste(resDir, 'RDS_Dir', sep = '/')
  # if (!dir.exists(rdsDir)) dir.create(rdsDir)
  ## -------
  ## 1. setup the original seurat object into a list for each item in the 'cellrangerResList'
  # Sys.time()
  print(sprintf('Step 1: read in 10X data into Seurat object at %s', Sys.time()))
  seuratObjList                  <- list()
  for (x in 1:length(cellrangerResList)) {
    print('---===------------')
    print(sprintf("Processing sample %s: '%s'.", x, as.character(cellrangerResList[[x]])))
    cellrangerCounts           <- Read10X(data.dir = cellrangerResList[[x]])
    print(sprintf('Orignially it has %s cells and %s features originated from cellranger to import into seurat object', length(cellrangerCounts@Dimnames[[2]]), length(cellrangerCounts@Dimnames[[1]]) ))
    ## include feature detected in at least 'min.cells = 3', and include cells where at least 'min.features = 200' detected
    seuratObjOrg               <- CreateSeuratObject(counts = cellrangerCounts,  project = names(cellrangerResList)[x], min.cells = as.numeric(minCells), min.features = as.numeric(minFeatures))
    ## add metadata feature2 into object, here is 'expCond' for metadata$sample
    seuratObjOrg               <- AddMetaData(object = seuratObjOrg,  col.name = 'expCond', metadata = as.factor(metadata$sample[x]))
    ## extra metadata columns in metadata columns 'expCond1' and 'expCond2' will be added respectively
    if ( 'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
      seuratObjOrg             <- AddMetaData(object = seuratObjOrg,  col.name = 'expCond1', metadata = as.factor(metadata$expCond1[x]))
      seuratObjOrg             <- AddMetaData(object = seuratObjOrg,  col.name = 'expCond2', metadata = as.factor(metadata$expCond2[x]))
    } else if ( 'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata) ) {
      seuratObjOrg             <- AddMetaData(object = seuratObjOrg,  col.name = 'expCond1', metadata = as.factor(metadata$expCond1[x]))
    } else if ( !'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
      seuratObjOrg             <- AddMetaData(object = seuratObjOrg,  col.name = 'expCond2', metadata = as.factor(metadata$expCond2[x]))
    }
    ## ---
    orgCellNoSummary           <- data.frame('cellrangeRcellNo' = dim(seuratObjOrg)[2], 'cellrangeRfeatureNo' = dim(seuratObjOrg)[1] )
    if (x == 1) {
      orgCellNoSummarySampComb <- orgCellNoSummary
    } else {
      orgCellNoSummarySampComb <- rbind(orgCellNoSummarySampComb, orgCellNoSummary)
    }
    ## ---
    doubletsMethod             <- doubletsRmMethods[x]
    if (doubletsMethod != 'none') {
      if (is.na(doubletsResDirs[x])) stop('doubletsMethod is on, please provide corresponding full path to the saved doublets removal results')
      print(sprintf('%s doublet removal methods were implemented.', doubletsMethod))
      if (doubletsMethod == 'centroids') {
        doubletFname           <- paste(doubletsResDirs[x], '/', names(cellrangerResList), '_centroids_doublet_cells_name.Rdata', sep = '')
        load(doubletFname[x])
        print(sprintf('%s doublets were estimated from %s sample', length(centroidsDoubletCells), names(cellrangerResList)[x] ))
        doubletCellsUpdate <- gsub(pattern = '[.]', replacement = '-', centroidsDoubletCells)
      } else if (doubletsMethod == 'medoids') {
        doubletFname           <- paste(doubletsResDirs[x], '/', names(cellrangerResList), '_medoids_doublet_cells_name.Rdata', sep = '')
        load(doubletFname[x])
        print(sprintf('%s doublets were estimated from %s sample', length(medoidsDoubletCells), names(cellrangerResList)[x] ))
        doubletCellsUpdate <- gsub(pattern = '[.]', replacement = '-', medoidsDoubletCells)
      } else if (doubletsMethod == 'OL') {
        doubletFname           <- paste(doubletsResDirs[x], '/', names(cellrangerResList), '_OL_doublet_cells_name.Rdata', sep = '')
        load(doubletFname[x])
        print(sprintf('%s doublets were estimated from %s sample', length(olDoubletCells), names(cellrangerResList)[x] ))
        doubletCellsUpdate <- gsub(pattern = '[.]', replacement = '-', olDoubletCells)
      }
      doubletDf                <- data.frame(cells = rownames(seuratObjOrg@meta.data))
      doubletDf                <- doubletDf %>% dplyr::mutate(doublet = ifelse(rownames((seuratObjOrg@meta.data)) %in% doubletCellsUpdate, 'TRUE', 'FASLE') )
      seuratObjOrg             <- AddMetaData(object = seuratObjOrg,  col.name = 'doublet', metadata = as.factor(doubletDf$doublet) )
      seuratObj                <- subset(seuratObjOrg, doublet == 'FASLE')
      print('Orignal seurat object without doublet removal:')
      print(seuratObjOrg)
      print('---')
      print(sprintf('%s Doublets were removal, %s cells were left', length(doubletCellsUpdate), (dim(seuratObjOrg@meta.data)[1] - length(doubletCellsUpdate)) ))
      print('doublet removed object:')
      print(seuratObj)
      print('---')
    } else {
      seuratObj                <- seuratObjOrg
      print('No doublet removal is executed here')
      print(seuratObjOrg)
      print('---')
    }
    seuratObjList[[x]]         <- seuratObj
  }
  print(sprintf('Step 1: END read in 10X data into Seurat object at %s', Sys.time()))
  print('---===---')
  names(seuratObjList)         <- names(cellrangerResList)
  cellNoSummary                <- data.frame('afterDoubletsRmCellNo' = unlist( lapply(seuratObjList, function(x) dim(x)[2])), 'afterDoubletsRmfeatureNo' = unlist(lapply(seuratObjList, function(x) dim(x)[1])) )
  print('---===---')
  print('Original cellrangR processed cell No. and features')
  print(orgCellNoSummarySampComb)
  print('---')
  print('To be processed cell No. and features')
  print(cellNoSummary)
  print('---===---')
  cellNoSummaryComb            <- cbind(orgCellNoSummarySampComb, cellNoSummary)
  rownames(cellNoSummaryComb)  <- names(cellrangerResList)
  write.table(x = cellNoSummaryComb, file = file.path(resDir, 'org_doubletsRemoval_cellNoSummary.txt'), quote = F, row.names = T, col.names = NA, sep = '\t')
  ##--------------------------------------------------------------------------------------
  ## 2. pre-processing: 1) check/filter the mitochodrial content, 2) normalization, 3) find variable features, and scaling
  print('Step 2: check mitochodrial content & normalization')
  ## 2.1 create 'seuratObjListPlotsDir' for each item in processed 'seuratObjList' based on input 'cellrangerResList'
  qcPlotsDir                  <- paste(resDir, 'QC_plots', sep = '/')
  if (!dir.exists(qcPlotsDir)) dir.create(qcPlotsDir)
  ## ---
  seuratQcProcessObjList      <- list()
  for (x in 1:length(seuratObjList)) {
    ## ---------
    seuratObj                 <- seuratObjList[[x]]
    ## 1).1 calculating mitochodrial content
    # print(sprintf('%s mitochodiral genes are processed', length(grep('^MT-', seuratObj@assays$RNA@data@Dimnames[[1]])) ))
    # seuratObj[['percent.mt']] <- PercentageFeatureSet(object = seuratObj, pattern = '^MT-')
    print(sprintf('%s mitochodiral genes are processed', length(grep(mtPatten(as.character(genomeSpecies)), seuratObj@assays$RNA@data@Dimnames[[1]])) ))
    seuratObj[['percent.mt']] <- PercentageFeatureSet(object = seuratObj, pattern = mtPatten(as.character(genomeSpecies)) )
    print(head(seuratObj@meta.data, 5))
    ## calculating no. of cells with certain mitochodrial percentage
    mtPer <- list()
    for(k in seq_along(1:4)){
      mtPer[[k]] <- sum(seuratObj@meta.data$percent.mt < 5*k) / length(seuratObj@meta.data$percent.mt)
    }
    mtPer                     <- do.call(rbind, mtPer)
    mtPerTab                  <- data.frame(`MT cutoff` = c("<5%", "<10%", "<15%", "<20%"), `Percentage of cells` =  scales::label_percent()(mtPer[,1]))
    print(mtPerTab)
    if (x == 1) {
      mtPerTabSampComb        <- mtPerTab
    } else {
      mtPerTabSampComb        <- dplyr::full_join(x = mtPerTabSampComb, y = mtPerTab, by = 'MT.cutoff')
    }
    # print(mtPerTabSampComb)
    # write.table(x = mtPerTab, file = paste(resDir, '/ToDelete_MT_percentage_summary_', names(seuratObjList)[x], '.txt', sep = ''), quote = F, sep = '\t', row.names = F, col.names = T)
    ## 1).2 plot mitochodrial content
    pdf(file = paste(qcPlotsDir, '/featureScatter_', names(seuratObjList)[x], '.pdf', sep = ''), width = 10, height = 6)
    # VlnPlot(object = seuratObj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
    plot1                     <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2                     <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + plot2)
    dev.off()
    ## 1). 3 visualize QC metrics as a violin plot
    pdf(file = paste(qcPlotsDir, '/featureViolin_', names(seuratObjList)[x], '.pdf', sep = ''), width = 10, height = 6)
    vlnFeaturePlot            <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(vlnFeaturePlot)
    dev.off()
    ## --
    if (mtFiltering) {
      ## Filtering is on, update 'seuratObj' with subsetted obj
      ## 1). 3 remove unwanted features and cells
      seuratObjBeforeFilter   <- seuratObj
      # seuratObj              <- subset(seuratObjBeforeFilter, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10) ##resDir: *_wFiltering_nFeatures, note 10% setup might be too low
      seuratObj               <- subset(seuratObjBeforeFilter, subset = nFeature_RNA > 200 & percent.mt < mtPerCutoff) ##resDir: *_wFiltering_mtContentOnly20_standardWorkflow
      ## 1).2 plot mitochodrial content
      pdf(file = paste(qcPlotsDir, '/featureScatter_', names(seuratObjList)[x], '_afterFiltering.pdf', sep = ''), width = 10, height = 6)
      # VlnPlot(object = seuratObj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
      plot1 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      print(plot1 + plot2)
      dev.off()
      pdf(file = paste(qcPlotsDir, '/featureViolin_', names(seuratObjList)[x], '_afterFiltering.pdf', sep = ''), width = 10, height = 6)
      vlnFeaturePlot <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      print(vlnFeaturePlot)
      dev.off()
      print(sprintf('Processing sample %s (%s), After filtering out high percentage contained mitochondrial genes, low and high gene feature expressed cells, cell number reduced from %s to %s.', x, names(seuratObjList)[x], dim(seuratObjBeforeFilter@meta.data)[1], dim(seuratObj@meta.data)[1] ))
      noFilteredCells           <- data.frame(`filtering` = c('before', 'after'), `No cells` = c(dim(seuratObjBeforeFilter@meta.data)[1], dim(seuratObj@meta.data)[1]) )

      # write.table(x = noFilteredCells, file = paste(resDir, '/ToDelete_No_mtFiltered_cells_summary_', names(seuratObjList)[x], '.txt', sep = ''), quote = F, sep = '\t', row.names = F, col.names = T)
      ## Complete filtering
    } else {
      noFilteredCells           <- data.frame(`filtering` = c('before'), `No cells` = c(dim(seuratObj@meta.data)[1]) )
    }
    if (x == 1) {
      noFilteredCellsSampComb <- noFilteredCells
    } else {
      noFilteredCellsSampComb <- dplyr::full_join(x = noFilteredCellsSampComb, y = noFilteredCells, by = 'filtering')
    }
    ## -
    ## 2). Normalization: 'normalization.method' & 'scale.factor' are default options
    seuratObj                   <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
    ## 3). Find variable features, by default top 2000
    ## by default, select top 2000 features, vst is also default method, other options are 'mean.var.plot(mvp)' and 'dispersion (disp)'
    seuratObj                   <- FindVariableFeatures(seuratObj, selection.method = 'vst', nfeatures = 5000)
    topFeatures                 <- head(VariableFeatures(seuratObj), n = 10)
    ## 3).2 plot top 100 variable features
    pdf(file = paste(qcPlotsDir, '/topVariableFeature_', names(seuratObjList)[x], '.pdf', sep = ''), width = 8, height = 6)
    plot1 <- VariableFeaturePlot(seuratObj)
    plot2 <- LabelPoints(plot = plot1, points = topFeatures, repel = TRUE)
    print(plot2)
    dev.off()
    seuratQcProcessObjList[[x]] <- seuratObj
  }
  names(seuratQcProcessObjList) <- names(seuratObjList)
  colnames(mtPerTabSampComb)    <- c(colnames(mtPerTabSampComb)[1], names(seuratObjList))
  write.table(x = mtPerTabSampComb, file = file.path(resDir, 'MT_percentage_summary.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
  colnames(noFilteredCellsSampComb)   <- c(colnames(noFilteredCellsSampComb)[1], names(seuratObjList))
  noFilteredCellsSampComb$Total       <- rowSums(noFilteredCellsSampComb[,-1])
  write.table(x = noFilteredCellsSampComb, file = file.path(resDir, 'No_filtered_cells_summary.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
  print('Step 2: END check mitochodrial content & normalization')
  print('---===---')
  ## ------
  return(list('countReadInOjb' = seuratObjList, 'qcProcessObj' = seuratQcProcessObjList, 'resDir' = resDir))
}
## ---------
## Minor fns to return mitochodrial content search pattern
mtPatten <- function(genomeSpecies) {
  if (genomeSpecies == 'human') return('^MT-')
  if (genomeSpecies == 'mouse') return('^mt-')
  if (genomeSpecies == 'human_mouse') return('^MT-|mt-')
}
meata2list <- function(metadata) {
  cellrangerResList    <- list()
  for (i in 1:length(metadata$sample)) {
    cellrangerResList[[i]] <- as.character(metadata$path[i])
  }
  names(cellrangerResList) <- metadata$sample
  return(cellrangerResList)
}
##----------------------------------------------------------------------------------------
