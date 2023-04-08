## processQC(): read in counts in MEX (cellranger results out format), hd5/hdf5, txt, or RDS formats into R
##              to establish the seurat object with inherited expCond* included
##--------------------------------------------------------------------------------------##
#' processQC() Function
#' @details
#' This function is used to read-in cellranger results into Seurat object by defining w/wo doublets removal.
#' @param metadata metadata table with at least 4 required columns: 'sample', 'path', 'doubletsRmMethod', and 'expCond1'.
#' @param multiomics logical option to specify whether multiomics sequencing data is used, by default = F.
#' @param extraFilter logical option to specify whether extra cells in excel file specified in the column 'filterFname' of metadata table need to be removded from the analysis, by default = F.
#' @param resDirName define the folder/directory name where integration analysis results will be saved, if not defined, by default results will be saved at the current working directory in a folder named as 'scRICA_results'.
#' @param genomeSpecies specify sample's genome species, by default 'human', currently supporting human, mouse, and rate.
#' @param minCells specify the minimum detected number of cells with an gene expression to be included in the analysis, this parameter is used by Seurat, bu default = 3.
#' @param minFeatures specify the minimum number of expressed gene features for cells per cell included in the analysis, this parameter is used by Seurat, by default = 200.
#' @param mtFiltering logical option to indicate whether to filter high mitochondrial content specified in option 'mtPerCutoff', by default 'FALSE'.
#' @param mtPerCutoff if 'mtFiltering = T', this option is required to indicate the percentage cut-off for mitochondrial content filtering.
#'
#' @importFrom Seurat Read10X
#' @importFrom tools file_ext
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
#' @examples processQC(metadata = metadata, resDirName = 'my_integration_results')
#' @export
#'
#' @return
#' the QC analysis results saved in defined 'resDir' under current working directory.
#' Meanwhile, this will also return a list item including 3 elements:
#' 1. 'countReadInOjb': a list of Seurat objects of original counts matrix without QC filtering.
#' 2. 'qcProcessObj': a list of Seurat objects with QC passed, such as doublets removal and mitochondrial content removal.
#' 3. 'resDir': full path of the directory/folder name where the QC results are saved.
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
processQC <- function(metadata, multiomics = F, extraFilter=F, resDirName=NULL, genomeSpecies=NULL, minCells=3, minFeatures=200, mtFiltering=F, mtPerCutoff=NULL, nfeatures = 5000) {
  ##--------------------------------------------------------------------------------------##
  if (!all(c("sample", "path", 'expCond1',"doubletsRmMethod" ) %in% colnames(metadata))) stop('Please provide metadata table with at least 4 columns: sample, path, expCond1, and doubletsRmMethod')
  ## ---
  if (is.null(resDirName)) resDirName <- as.character('scRICA_results')
  # if (!is.list(cellrangerResList)) stop("Please provide list item of 'cellrangerResList', which is required.")
  minCells                       <- as.numeric(minCells)
  minFeatures                    <- as.numeric(minFeatures)
  mtFiltering                    <- as.logical(mtFiltering)
  multiomics                     <- as.logical(multiomics)
  extraFilter                    <- as.logical(extraFilter)
  if (is.null(genomeSpecies)) genomeSpecies <- as.character('human')
  if (mtFiltering & is.null(mtPerCutoff)) stop("Mitochondrial content filtering option ('mtFiltering') is on, please provide corresponding Mitochondrial content percentage filtering option in 'mtPerCutoff'.  ")
  mtPerCutoff                    <- as.numeric(mtPerCutoff)
  ##--------------------------------------------------------------------------------------##
  # if ( 'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata)) {
  #   print(sprintf('Only 1 experimental condition are provided with %s experimental factor levels for comparisons.', length(levels(factor(metadata$expCond1))) ))
  # } else if ( 'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
  #   print(sprintf('2 experimental conditions are provided, experimental condition 1 has %s experimental factor levels, and experimental condition 2 has %s experimental factor levels for comparisons', length(levels(factor(metadata$expCond1))), length(levels(factor(metadata$expCond2))) ))
  # } else if ( !'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata) ) {
  #   print(sprintf('No experimental condition factors are provided, the analysis will be conducted only based samples integration'))
  # } else if ( !'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata)) {
  #   print(sprintf('Only 1 experimental condition are provided with %s experimental factor levels for comparisons.', length(levels(factor(metadata$expCond2))) ))
  # }
  print("-=-=-=-=-=-=-")
  expCond.no    <- length(grep('expCond', colnames(metadata)))
  expCond.index <- grep('expCond', colnames(metadata))
  print(sprintf("%s experimental condition are provided in the metadata columns: %s.", expCond.no, paste(colnames(metadata)[expCond.index], collapse = ', ')) )
  for (i in 1:expCond.no) {
    print(sprintf("%s: %s has %s experimental factor levels", i, colnames(metadata)[expCond.index], length(levels(factor(metadata[,expCond.index[i]]))) ))
  }
  print("-=-=-=-=-=-=-")
  ##--------------------------------------------------------------------------------------##
  ## 0. create main directory for results to be save in
  resDir                         <- paste(getwd(), resDirName, sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  ## run findDoublet and update metadata to include doublets detection results
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
      print('START: Doublets identification analysis before processing to the next step.')
      doubletsRes                <- findDoublets(metadata = md4doublet, multiomics = multiomics, genomeSpecies = genomeSpecies, resFilename = paste(resDirName, 'doublet_results', sep = '/') )
      print('END: Doublets identification, process to next step QC.')
      print('===================================================================')
      md4doublet$doubletsResDir  <- rep(as.character(doubletsRes), length(md4doublet$sample))
      doubletsMd                 <- md4doublet %>% dplyr::select(c('sample', 'doubletsResDir'))
      metadata                   <- dplyr::left_join(x = metadataOrg, y = doubletsMd, by = 'sample', )
    }
  } else {
    print('===================================================================')
    print('NO doublets identification analysis is needed, because doublets identification either conducted or not needed.')
    metadataOrg                  <- metadata
    metadataOrg$doubletsRmMethod <- gsub('[]|[ ]', 'none', metadataOrg$doubletsRmMethod)
    metadataOrg$doubletsRmMethod[is.na(metadataOrg$doubletsRmMethod)] <- 'none'
    metadataOrg$doubletsRmMethod[metadataOrg$doubletsRmMethod=='NONE'] <- 'none'
    metadataOrg$doubletsRmMethod[metadataOrg$doubletsRmMethod=='None'] <- 'none'
    metadataOrgDbNotNa           <- metadataOrg %>% dplyr::filter(doubletsRmMethod != 'None') %>% dplyr::filter(doubletsRmMethod != 'NONE') %>% dplyr::filter(doubletsRmMethod != 'none')
    if(dim(metadataOrgDbNotNa)[1]>0) {
      if (!all(dir.exists(metadataOrgDbNotNa$doubletsResDir))){
        print(sprintf("Stop at sample %s:", metadata$sample[which(!dir.exists(metadataOrgDbNotNa$doubletsResDir))]))
        stop("please provide correct corresponding 'doubletsResDir' in 'metadata' for column 'doubletsRmMethod' specified NOT as 'none'.")
      }
    }
    metadata                     <- metadataOrg
    print('Process to QC.')
    print('===================================================================')
  }
  # print(metadata)
  # print('*****************************')
  ##--------------------------------------------------------------------------------------##
  cellrangerResList              <- meata2list(metadata = metadata) ## by default names(cellrangerResList) = metadata$sample
  ## ---
  doubletsRmMethods              <- as.character(metadata$doubletsRmMethod)
  doubletsResDirs                <- as.character(metadata$doubletsResDir)
  ## 1. setup the original Seurat object into a list for each item in the 'cellrangerResList'
  # Sys.time()
  print(sprintf('Step 1: read in 10X data into Seurat object at %s', Sys.time()))
  seuratObjList                  <- list()
  for (x in 1:length(cellrangerResList)) {
    print('---===------------')
    print(sprintf("Processing sample %s: '%s'.", x, as.character(cellrangerResList[[x]])))
    if (dir.exists(cellrangerResList[[x]])) {
      print("read in count data in MEX format")
      cellrangerCountsOrg          <- Seurat::Read10X(data.dir = cellrangerResList[[x]])
    } else {
      if (tools::file_ext(cellrangerResList[[x]])=='h5' | tools::file_ext(cellrangerResList[[x]]) == 'hdf5') {
        print("read in count data in h5/hdf5 format")
        cellrangerCountsOrg        <- Seurat::Read10X_h5(filename = as.character(cellrangerResList[[x]]), use.names = TRUE, unique.features = TRUE)
      } else if (tools::file_ext(cellrangerResList[[x]])=='txt' | tools::file_ext(gsub('.gz', '', cellrangerResList[[x]])) == 'txt') {
        print("read in count data in txt format")
        cellrangerCountsOrg        <- read.delim2(file = as.character(cellrangerResList[[x]]))
      } else if (tools::file_ext(cellrangerResList[[x]])=='rds') {
        print("read in count data in rds format")
        rdsOrg                     <- readRDS(file = as.character(cellrangerResList[[x]]))
        cellrangerCountsOrg        <- rdsOrg@assays$RNA@counts
      } else {
        stop("input file in metadata table cannot be read in, it should be any of these 3 formats: txt, hdf5, or MEX in a directory.")
      }
    }
    ##--------------------------------------------------------------------------------------##
    if (multiomics) {
      cellrangerCounts           <- cellrangerCountsOrg$`Gene Expression`
    } else {
      cellrangerCounts           <- cellrangerCountsOrg
    }
    # print(colnames(cellrangerCounts)[1:5])
    ## -
    if (dir.exists(cellrangerResList[[x]]) | tools::file_ext(cellrangerResList[[x]])=='h5' | tools::file_ext(cellrangerResList[[x]]) == 'hdf5') {
      print(sprintf('Originally it has %s cells and %s features originated from cellranger (MEX or hdf5 format) to import into Seurat object', length(cellrangerCounts@Dimnames[[2]]), length(cellrangerCounts@Dimnames[[1]]) ))
    } else if (tools::file_ext(cellrangerResList[[x]])=='txt' | tools::file_ext(gsub('.gz', '', cellrangerResList[[x]])) == 'txt') {
      print(sprintf('Originally it has %s cells and %s features originated from a txt format file to import into Seurat object', dim(cellrangerCounts)[2], dim(cellrangerCounts)[1] ))
    }
    ##--------------------------------------------------------------------------------------##
    if(extraFilter) {
      if (!'filterFname' %in% colnames(metadata)) stop("Option extraFilter is on, but no filter files is provided in the metadata table column 'filterFname'.")
      if (is.na(metadata$filterFname[x])) {
        print(sprintf("'extraFilter' is on for other samples, but not here; %s cells will be used for next steps analysis.", length(cellrangerCounts@Dimnames[[2]])))
        cellrangerCounts <- cellrangerCounts
      } else {
        if (tools::file_ext(metadata$filterFname[x])=='csv') {
          filterRes   <- read.csv(file = metadata$filterFname[x], header = T)
        } else  if (tools::file_ext(metadata$filterFname[x])=='txt') {
          filterRes   <- read.delim(file = metadata$filterFname[x], header = T, sep = '\t')
        }
        if (sum(grepl('barcode|filter',colnames(filterRes)))!=2) stop("Please make sure columns 'barcode' & 'filter' are inside provided ''.")
        if (!all(filterRes$barcode %in% colnames(cellrangerCounts) )) stop("provided filter files does not include all corresponding cells information.")
        filter.index  <- grep('filter', colnames(filterRes))
        print(sprintf("%s cells will be filtered based on input 'filterFname' in metadata table; %s cells will be used for next steps analysis.",
                      sum(filterRes[,filter.index]==as.logical(T)), sum(filterRes[,filter.index]==as.logical(F)) ))
        cellrangerCountsFilter <- cellrangerCounts[, match(filterRes$barcode[filterRes[,filter.index]==as.logical(F)], colnames(cellrangerCounts))]
        cellrangerCounts       <- cellrangerCountsFilter
      }
    }
    ##--------------------------------------------------------------------------------------##
    seuratObjOrg               <- Seurat::CreateSeuratObject(counts = cellrangerCounts,  project = names(cellrangerResList)[x], min.cells = as.numeric(minCells), min.features = as.numeric(minFeatures))
    ## add metadata feature2 into object, here is 'expCond' for metadata$sample
    seuratObjOrg               <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond', metadata = as.factor(metadata$sample[x]))
    # print(head(seuratObjOrg@meta.data)) ##for debug
    ## extra metadata columns in metadata columns 'expCond*' will be added respectively
    for (i in 1:expCond.no ) {
      seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = sprintf("expCond%s",i), metadata = as.factor(metadata[x, expCond.index[i]]) )
    }
    # print(head(seuratObjOrg@meta.data)) ##for debug
    # print(table(seuratObjOrg@meta.data$expCond3)) ##for debug
    ## -------------------------
    # if ( 'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
    #   seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond1', metadata = as.factor(metadata$expCond1[x]))
    #   seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond2', metadata = as.factor(metadata$expCond2[x]))
    # } else if ( 'expCond1' %in% colnames(metadata) & !'expCond2' %in% colnames(metadata) ) {
    #   seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond1', metadata = as.factor(metadata$expCond1[x]))
    # } else if ( !'expCond1' %in% colnames(metadata) & 'expCond2' %in% colnames(metadata) ) {
    #   seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'expCond2', metadata = as.factor(metadata$expCond2[x]))
    # }
    ## -------------------------
    ##--------------------------------------------------------------------------------------##
    orgCellNoSummary           <- data.frame('cellrangeRcellNo' = dim(seuratObjOrg)[2], 'cellrangeRfeatureNo' = dim(seuratObjOrg)[1] )
    if (x == 1) {
      orgCellNoSummarySampComb <- orgCellNoSummary
    } else {
      orgCellNoSummarySampComb <- rbind(orgCellNoSummarySampComb, orgCellNoSummary)
    }
    ##--------------------------------------------------------------------------------------##
    doubletsMethod             <- tolower(doubletsRmMethods[x])
    if (doubletsMethod == 'none') {
      seuratObj                <- seuratObjOrg
      print('No doublet removal is executed here')
      print(seuratObjOrg)
      print('---')
    } else {
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
      } else if (doubletsMethod == 'ol') {
        doubletFname           <- paste(doubletsResDirs[x], '/', names(cellrangerResList), '_OL_doublet_cells_name.Rdata', sep = '')
        load(doubletFname[x])
        print(sprintf('%s doublets were estimated from %s sample', length(olDoubletCells), names(cellrangerResList)[x] ))
        doubletCellsUpdate <- gsub(pattern = '[.]', replacement = '-', olDoubletCells)
      }
      doubletDf                <- data.frame(cells = rownames(seuratObjOrg@meta.data))
      doubletDf                <- doubletDf %>% dplyr::mutate(doublet = ifelse(rownames((seuratObjOrg@meta.data)) %in% doubletCellsUpdate, 'TRUE', 'FASLE') )
      seuratObjOrg             <- Seurat::AddMetaData(object = seuratObjOrg,  col.name = 'doublet', metadata = as.factor(doubletDf$doublet) )
      seuratObj                <- subset(seuratObjOrg, doublet == 'FASLE')
      print('Original Seurat object without doublet removal:')
      print(seuratObjOrg)
      print('---')
      print(sprintf('%s Doublets were removal, %s cells were left', length(doubletCellsUpdate), (dim(seuratObjOrg@meta.data)[1] - length(doubletCellsUpdate)) ))
      print('doublet removed object:')
      print(seuratObj)
      print('---')
    }
    seuratObjList[[x]]         <- seuratObj
  }
  print(sprintf('Step 1: END read in 10X data into Seurat object at %s', Sys.time()))
  print('---===---')
  ##--------------------------------------------------------------------------------------##
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
  ##--------------------------------------------------------------------------------------##
  ## 2. pre-processing: 1) check/filter the mitochondrial content, 2) normalization, 3) find variable features, and scaling
  print('Step 2: check mitochondrial content & normalization')
  ## 2.1 create 'seuratObjListPlotsDir' for each item in processed 'seuratObjList' based on input 'cellrangerResList'
  qcPlotsDir                  <- paste(resDir, 'QC_plots', sep = '/')
  if (!dir.exists(qcPlotsDir)) dir.create(qcPlotsDir)
  ## ---
  seuratQcProcessObjList      <- list()
  for (x in 1:length(seuratObjList)) {
    ## ---------
    seuratObj                 <- seuratObjList[[x]]
    ## 1).1 calculating mitochondrial content
    # print(sprintf('%s mitochondrial genes are processed', length(grep('^MT-', seuratObj@assays$RNA@data@Dimnames[[1]])) ))
    # seuratObj[['percent.mt']] <- PercentageFeatureSet(object = seuratObj, pattern = '^MT-')
    print(sprintf('%s mitochondrial genes are processed', length(grep(mtPatten(as.character(genomeSpecies)), seuratObj@assays$RNA@data@Dimnames[[1]])) ))
    seuratObj[['percent.mt']] <- Seurat::PercentageFeatureSet(object = seuratObj, pattern = mtPatten(as.character(genomeSpecies)) )
    print(head(seuratObj@meta.data, 5))
    ## calculating no. of cells with certain mitochondrial percentage
    mtPer <- list()
    rRNAper <- list()
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
    ## ---
    ## calculate rRNA content
    print(sprintf('%s rRNA genes are processed', length(grep(rRNAcontent(as.character(genomeSpecies)), seuratObj@assays$RNA@data@Dimnames[[1]])) ))
    seuratObj[['rRNA.content']] <- Seurat::PercentageFeatureSet(object = seuratObj, pattern = rRNAcontent(as.character(genomeSpecies)) )
    for(k in seq_along(1:20)){
      rRNAper[[k]] <- sum(seuratObj@meta.data$rRNA.content < 5*k) / length(seuratObj@meta.data$rRNA.content)
    }
    rRNAper                     <- do.call(rbind, rRNAper)
    rRNAperTab                  <- data.frame(`rRNA content` = c("<5%", "<10%", "<15%", "<20%",
                                                                 "<25%", "<30%", "<35%", "<40%",
                                                                 "<45%", "<50%", "<55%", "<60%",
                                                                 "<65%", "<70%", "<75%", "<80%",
                                                                 "<85%", "<90%", "<95%", "<100%"),
                                              `Percentage of cells` =  scales::label_percent()(rRNAper[,1]))
    if (x == 1) {
      rRNAperSampComb           <- rRNAperTab
    } else {
      rRNAperSampComb           <- dplyr::full_join(x = rRNAperSampComb, y = rRNAperTab, by = 'rRNA.content')
    }
    ## ------
    # print(rRNAperSampComb)
    # write.table(x = mtPerTab, file = paste(resDir, '/ToDelete_MT_percentage_summary_', names(seuratObjList)[x], '.txt', sep = ''), quote = F, sep = '\t', row.names = F, col.names = T)
    ## 1).2 plot mitochondrial content
    pdf(file = paste(qcPlotsDir, '/featureScatter_', names(seuratObjList)[x], '.pdf', sep = ''), width = 10, height = 6)
    # VlnPlot(object = seuratObj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
    plot1                     <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2                     <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + plot2)
    dev.off()
    ## 1). 3 visualize QC metrics as a violin plot
    pdf(file = paste(qcPlotsDir, '/featureViolin_', names(seuratObjList)[x], '.pdf', sep = ''), width = 10, height = 6)
    vlnFeaturePlot            <- Seurat::VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rRNA.content"), ncol = 4)
    print(vlnFeaturePlot)
    dev.off()
    ## --
    if (mtFiltering) {
      ## Filtering is on, update 'seuratObj' with subsetted obj
      ## 1). 3 remove unwanted features and cells
      seuratObjBeforeFilter   <- seuratObj
      seuratObj               <- subset(seuratObjBeforeFilter, subset = nFeature_RNA > 200 & percent.mt < mtPerCutoff)
      ## 1).2 plot mitochondrial content
      pdf(file = paste(qcPlotsDir, '/featureScatter_', names(seuratObjList)[x], '_afterFiltering.pdf', sep = ''), width = 10, height = 6)
      # VlnPlot(object = seuratObj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
      plot1 <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      print(plot1 + plot2)
      dev.off()
      pdf(file = paste(qcPlotsDir, '/featureViolin_', names(seuratObjList)[x], '_afterFiltering.pdf', sep = ''), width = 10, height = 6)
      vlnFeaturePlot <- Seurat::VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rRNA.content"), ncol = 4)
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
    seuratObj                   <- Seurat::NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
    ## 3). Find variable features, by default top 2000
    ## by default, select top 2000 features, vst is also default method, other options are 'mean.var.plot(mvp)' and 'dispersion (disp)'
    seuratObj                   <- Seurat::FindVariableFeatures(seuratObj, selection.method = 'vst', nfeatures = nfeatures)
    topFeatures                 <- head(Seurat::VariableFeatures(seuratObj), n = 10)
    ## 3).2 plot top 100 variable features
    pdf(file = paste(qcPlotsDir, '/topVariableFeature_', names(seuratObjList)[x], '.pdf', sep = ''), width = 8, height = 6)
    plot1 <- Seurat::VariableFeaturePlot(seuratObj)
    plot2 <- Seurat::LabelPoints(plot = plot1, points = topFeatures, repel = TRUE)
    print(plot2)
    dev.off()
    seuratQcProcessObjList[[x]] <- seuratObj
  }
  names(seuratQcProcessObjList) <- names(seuratObjList)
  colnames(mtPerTabSampComb)    <- c(colnames(mtPerTabSampComb)[1], names(seuratObjList))
  write.table(x = mtPerTabSampComb, file = file.path(resDir, 'MT_percentage_summary.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
  ## ---
  colnames(rRNAperSampComb)    <- c(colnames(rRNAperSampComb)[1], names(seuratObjList))
  write.table(x = rRNAperSampComb, file = file.path(resDir, 'rRNA_content_summary.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
  ## ---
  colnames(noFilteredCellsSampComb)   <- c(colnames(noFilteredCellsSampComb)[1], names(seuratObjList))
  if (length(seuratQcProcessObjList)>1) {
    noFilteredCellsSampComb$Total       <- rowSums(noFilteredCellsSampComb[,-1])
  }
  write.table(x = noFilteredCellsSampComb, file = file.path(resDir, 'No_filtered_cells_summary.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
  print('Step 2: END check mitochondrial content & normalization')
  print('---===---')
  ##--------------------------------------------------------------------------------------##
  return(list('countReadInOjb' = seuratObjList, 'qcProcessObj' = seuratQcProcessObjList, 'resDir' = resDir))
}
##--------------------------------------------------------------------------------------##
## Minor fns to return mitochondrial content search pattern
mtPatten <- function(genomeSpecies) {
  if (genomeSpecies == 'human') return('^MT-')
  if (genomeSpecies == 'mouse') return('^mt-')
  if (genomeSpecies == 'human_mouse') return('^MT-|mt-')
}

rRNAcontent <-  function(genomeSpecies) {
  if (genomeSpecies == 'human') return('^RP[SL]')
  if (genomeSpecies == 'mouse') return('^rp[sl]')
  if (genomeSpecies == 'human_mouse') return('^RP[SL]|^rp[sl]')
}
meata2list <- function(metadata) {
  cellrangerResList    <- list()
  for (i in 1:length(metadata$sample)) {
    cellrangerResList[[i]] <- as.character(metadata$path[i])
  }
  names(cellrangerResList) <- metadata$sample
  return(cellrangerResList)
}
##--------------------------------------------------------------------------------------##
