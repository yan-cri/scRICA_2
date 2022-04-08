## R functions used to identify doublets/multiplets with DoubletDecon                   ##
## Developed by Yan Li, Jan, 2021                                                       ##
##--------------------------------------------------------------------------------------##
#' findDoublets() Function
#'
#' This function allows you to identify doublets with DoubletDecon for the provide input in metadata
#' @param metadata required, metadata table with at least 2 required columns: 'sample' and 'path', where path shows the count table matrix is located.
#' @param multiomics optional, default = F
#' @param extraFilter optional, default = F, if 'TRUE', column 'filterFname' needs to be included in the provided metadata table.
#' @param genomeSpecies genome species, currently supporting human, mouse, and rate.
#' @param doubletDeconRhop option used for DoubletDecon detection, by default = 0.5.
#' @param doubletDeconPMF option used for DoubletDecon detection, by default 'FALSE'.
#' @param doubletDeconNoCore option used for DoubletDecon detection, by default '-1' no parallel.
#' @param resFilename define the folder/directory name where doublets detection results will be saved, if not defined, by default results will be saved at the current working directory in a folder named as 'doublets_results'.
#'
#' @importFrom SeuratObject UpdateSeuratObject
#' @importFrom tools file_ext
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom DoubletDecon Main_Doublet_Decon
#' @keywords findDoublets
#' @export
#' @examples findDoublets()
#' @return generate a result directory/folder specified with 'resFilename' to save
#' doublets identification analysis results for all samples in the metadata with a data structure shown as below:
#' 1. 'metadata$sample' + '_medoids_doublet_cells_name.Rdata': detected doublet cells with mediods, object named as 'centroidsDoubletCells'
#' 2. 'metadata$sample' + '_centroids_doublet_cells_name.Rdata': detected doublet cells with centroids, object named as 'centroidsDoubletCells'
#' 3. 'metadata$sample' + '_OL_doublet_cells_name.Rdata': detected overlapped doublet cells with centroids and mediods, object named as 'olDoubletCells'
#' 4. 'metadata$sample' + '_doubletDeconRes.Rdata': Main_Doublet_Decon() analysis results for Centroids (object named as 'doubletDeconResCentroids') and Medoids (object named as 'doubletDeconResMedoids')
#' 5. 'metadata$sample' + '_doubletDecon_mediod_centroid_vennComp.pdf': venn-diagram plot of detected doublets with centroids and mediods algorithms
#' 6. 'metadata$sample' + '_doublets_no_summary.txt': centroids and mediods algorithms detected doublets summary
##----------------------------------------------------------------------------------------
# library(Seurat)
# library(DoubletDecon)
# library(dplyr)
# library(gplots)
# Sys.setenv('R_MAX_VSIZE'=32000000000)
##----------------------------------------------------------------------------------------
findDoublets <- function(metadata, multiomics = F, extraFilter=F, genomeSpecies=NULL, doubletDeconRhop=0.5, doubletDeconPMF=F, doubletDeconNoCore=-1, resFilename=NULL) {
  ## ---
  if (is.null(genomeSpecies)) genomeSpecies <- 'human'
  if (is.null(resFilename)) resFilename <- 'doublets_results'
  doubletDeconRhop               <- as.numeric(doubletDeconRhop)
  doubletDeconPMF                <- as.logical(doubletDeconPMF)
  multiomics                     <- as.logical(multiomics)
  extraFilter                    <- as.logical(extraFilter)
  ## ---
  cellrangerResList              <- meata2list(metadata = metadata)
  ## ---
  ## prepare results saving directory, if not exist, create one
  resDir  <- sprintf('%s/%s', getwd(), resFilename)
  if(!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('Doublets identification results will be saved in %s', resDir))
  ## intermediate 'resProcessDir' under/inside provided 'resDir'
  resProcessDir                  <- paste(resDir, 'doubletDecon_preProcessed_results', sep = '/')
  if(!dir.exists(resProcessDir)) dir.create(resProcessDir)
  ## ---
  ## loop over provided 'cellrangerResList' to identify/estimate doublets in each list item of provided 'cellrangerResList'
  for ( x in 1:length(cellrangerResList)) {
    ## ---
    ## 1. creat seurat object as DoubletDecon suggested
    print('---')
    print(sprintf("Processing sample '%s'.", as.character(cellrangerResList[[x]])))
    cellrangerCountsOrg          <- Seurat::Read10X(data.dir = cellrangerResList[[x]])
    if (multiomics) {
      cellrangerCounts           <- cellrangerCountsOrg$`Gene Expression`
    } else {
      cellrangerCounts           <- cellrangerCountsOrg
    }
    print(sprintf('Step1: Orignially it has %s cells and %s features originated from cellranger to import into seurat object', length(cellrangerCounts@Dimnames[[2]]), length(cellrangerCounts@Dimnames[[1]]) ))
    ## include feature detected in at least 'min.cells = 3', and include cells where at least 'min.features = 200' detected
    ## ---
    if(extraFilter) {
      if (!'filterFname' %in% colnames(metadata)) stop("Option extraFilter is on, but no filter files is provided in the metadata table column 'filterFname'.")
      if (file_ext(metadata$filterFname[x])=='csv') {
        filterRes   <- read.csv(file = metadata$filterFname[x], header = T)
      } else  if (file_ext(metadata$filterFname[x])=='txt') {
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
    ## ---
    seuratObject                 <- Seurat::CreateSeuratObject(counts = cellrangerCounts,  project = names(cellrangerResList)[x], min.cells = 3, min.features = 200)
    ## add metadata feature into object, here is 'expCond'
    seuratObject                 <- Seurat::AddMetaData(object = seuratObject,  col.name = 'expCond', metadata = as.factor(names(cellrangerResList)[x]))
    ## -
    seuratObject[['percent.mt']] <- Seurat::PercentageFeatureSet(object = seuratObject, pattern = as.character(mtPatten(genomeSpecies)) )
    seuratObject                 <- Seurat::NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
    print(sprintf('Complete Log Normlization'))
    seuratObject                 <- Seurat::FindVariableFeatures(seuratObject, selection.method = 'vst', nfeatures = 2000)
    print(sprintf('Complete vst top variable feature identification'))
    if (dim(seuratObject@meta.data)[1] > 10000) {
      seuratObject               <- Seurat::ScaleData( object = seuratObject, features = VariableFeatures(seuratObject) )
    } else {
      seuratObject               <- Seurat::ScaleData(object = seuratObject, features = rownames(seuratObject))
    }
    print(sprintf('Complete data centering/scaling'))
    seuratObject                 <- Seurat::RunPCA(object = seuratObject)
    print(sprintf('Complete PCA clustering'))
    seuratObject                 <- Seurat::RunUMAP(seuratObject, dims = 1:10)
    print(sprintf('Complete UMAP clustering'))
    seuratObject                 <- Seurat::RunTSNE(seuratObject, dims = 1:10)
    print(sprintf('Complete tSNE clustering'))
    seuratObject                 <- Seurat::FindNeighbors(seuratObject, dims = 1:10)
    print(sprintf('Complete finding neighbors'))
    seuratObject                 <- Seurat::FindClusters(seuratObject, resolution = 0.8)
    print(sprintf('Complete finding clustering'))
    print('END step1 for orignal seurat processing')
    ## 2. Improved_Seurat_Pre_Process() from DoubletDecon on established seurat object, and output results to 'resDir' for next step usage
    print(sprintf("Start Step2: 'improve seurat pre process'."))
    filename                     <- names(cellrangerResList)[x]
    newFiles                     <- Improved_Seurat_Pre_Process(seuratObject = seuratObject, num_genes=50, write_files=FALSE)
    write.table(newFiles$newExpressionFile, paste0(resProcessDir, '/', filename, "_expression"), sep="\t")
    write.table(newFiles$newFullExpressionFile, paste0(resProcessDir, '/', filename, "_fullExpression"), sep="\t")
    write.table(newFiles$newGroupsFile, paste0(resProcessDir, '/', filename , "_groups"), sep="\t", col.names = F)
    print(sprintf("END Step2: 'improve seurat pre process'."))
    ## 3.1 DoubletDecon doublets detection with Main_Doublet_Decon() based on centroid method
    print('=========')
    print(sprintf("Start Step3: DoubletDecon (centroids & medoids) use rhop = %s based on %s genome with PMF = %s", doubletDeconRhop, as.character(doubletDeconSpecies(genomeSpecies)), doubletDeconPMF ))
    print('Start Step 3.1: DoubletDecon centroid detection')
    doubletDeconResCentroids     <- DoubletDecon::Main_Doublet_Decon(rawDataFile = paste0(resProcessDir, '/', filename, "_expression"),
                                                                     groupsFile = paste0(resProcessDir, '/', filename , "_groups"),
                                                                     filename = filename,
                                                                     location = paste(resProcessDir, 'Centroids_', sep = '/'),
                                                                     removeCC = F, ## default is FALSE
                                                                     species = as.character(doubletDeconSpecies(genomeSpecies)), ## default is 'mmu'
                                                                     rhop = doubletDeconRhop, ## Default is 1, x in mean+x*SD to determine upper cutoff for correlation in the blacklist.
                                                                     PMF = doubletDeconPMF, ## default = T, Use step 2 (unique gene expression) in doublet determination criteria
                                                                     useFull =  F, ## default = F, Use full gene list for PMF analysis
                                                                     heatmap = F,
                                                                     centroids = TRUE,
                                                                     nCores = doubletDeconNoCore)
    print('Doublet centroids detection results table:')
    print( table(doubletDeconResCentroids$DRS_doublet_table$isADoublet) )
    print('END Step 3.1: DoubletDecon centroid detection')
    ## 3.2 DoubletDecon doublets detection with Main_Doublet_Decon() based on medoids method
    print('Start Step 3.2: DoubletDecon medoids detection')
    print('=========')
    doubletDeconResMedoids       <- DoubletDecon::Main_Doublet_Decon(rawDataFile = paste0(resProcessDir, '/', filename, "_expression"),
                                                                     groupsFile = paste0(resProcessDir, '/', filename , "_groups"),
                                                                     filename = filename,
                                                                     location = paste(resProcessDir, 'Medoids_', sep = '/'),
                                                                     removeCC = F, ## default is FALSE
                                                                     species = as.character(doubletDeconSpecies(genomeSpecies)), ## default is 'mmu'
                                                                     rhop = doubletDeconRhop, ## x in mean+x*SD to determine upper cutoff for correlation in the blacklist. Default is 1
                                                                     PMF = doubletDeconPMF, ## default = T, Use step 2 (unique gene expression) in doublet determination criteria
                                                                     useFull =  F, ## default = F, Use full gene list for PMF analysis
                                                                     heatmap = F,
                                                                     centroids = FALSE,
                                                                     nCores = doubletDeconNoCore)
    print('Doublet medoids detection results table:')
    print( table(doubletDeconResMedoids$DRS_doublet_table$isADoublet) )
    print('=========')
    print('END Step 3.2: DoubletDecon medoids detection')
    print(sprintf("Doublets identification results ('doubletDeconResCentroids' & 'doubletDeconResMedoids') saved in '%s'.", as.character(file.path(resDir, sprintf('%s_doubletDeconRes.Rdata', filename)))))
    save(doubletDeconResCentroids, doubletDeconResMedoids, file = file.path(resDir, sprintf('%s_doubletDeconRes.Rdata', filename)) )
    ## -
    print('Start Step 4: Doublets identification results summary')
    ## 4. summarize doubletDecon detection results.
    ## 4.1 Medoids detection summary
    doubletNoSummaryMedoids             <- as.data.frame(table(doubletDeconResMedoids$DRS_doublet_table$isADoublet))
    doubletNoSummaryMedoids$Per         <- round(x = doubletNoSummaryMedoids[,2] / sum(doubletNoSummaryMedoids[,2]) * 100, digits = 1)
    colnames(doubletNoSummaryMedoids)   <- c('Doublet', 'Medioids detection No', 'Medioids detection Per')
    print('Doublet detection results table with medoids:')
    print( table(doubletDeconResMedoids$DRS_doublet_table$isADoublet) )
    ## Centroids detection
    doubletNoSummaryCentroids           <- as.data.frame(table(doubletDeconResCentroids$DRS_doublet_table$isADoublet))
    doubletNoSummaryCentroids$Per       <- round(x = doubletNoSummaryCentroids[,2] / sum(doubletNoSummaryCentroids[,2]) * 100, digits = 1)
    colnames(doubletNoSummaryCentroids) <- c('Doublet', 'Centroids detection No', 'Centroids detection Per')
    print('Doublet detection results table with centroids:')
    print( table(doubletDeconResCentroids$DRS_doublet_table$isADoublet) )
    if (sum(doubletNoSummaryMedoids[,2]) != sum(doubletNoSummaryCentroids[,2])) stop('Error: total number of cells different from medoids and centroids detection methods.')
    medoidsDoubletCells = rownames(doubletDeconResMedoids$DRS_doublet_table %>% dplyr:: filter(isADoublet == TRUE))
    save(medoidsDoubletCells, file = file.path(resDir, sprintf("%s_medoids_doublet_cells_name.Rdata", filename) ) )
    ## -
    centroidsDoubletCells = rownames(doubletDeconResCentroids$DRS_doublet_table %>% dplyr:: filter(isADoublet == TRUE))
    save(centroidsDoubletCells, file = file.path(resDir, sprintf("%s_centroids_doublet_cells_name.Rdata", filename) ) )
    ## -
    print('Start overlapping medoids and centroids doublet detection results')
    olRes <- gplots::venn(list('medoids'  = rownames(doubletDeconResMedoids$DRS_doublet_table %>% dplyr::filter(isADoublet == TRUE)),
                               'centroid' = rownames(doubletDeconResCentroids$DRS_doublet_table %>% dplyr::filter(isADoublet == TRUE)) ))
    plot(olRes)
    ## -
    pdf(file = file.path(resDir, sprintf("%s_doubletDecon_mediod_centroid_vennComp.pdf", filename) ), width = 3, height = 4 )
    plot(olRes)
    dev.off()
    ## -
    olDoubletCells <- attr(olRes, 'intersection')[[grep('medoids:centroid', names(attr(olRes, 'intersection')) )]]
    print('Finish overlapping medoids and centroids doublet detection results')
    ## -
    doubletNoSummaryComb                <- data.frame(Doublet = c('TRUE', 'FALSE'),
                                                      No = c(length(olDoubletCells), sum(doubletNoSummaryCentroids[,2])-length(olDoubletCells) ))
    doubletNoSummaryComb$Per            <- round(x = doubletNoSummaryComb[,2] / sum(doubletNoSummaryComb[,2]) * 100, digits = 1)
    colnames(doubletNoSummaryComb)      <- c('Doublet', 'Combined detection No', 'Combined detection Per')
    doubletNoSummary                    <- dplyr::left_join(doubletNoSummaryMedoids, doubletNoSummaryCentroids, by = 'Doublet') %>% dplyr::left_join(doubletNoSummaryComb, by = 'Doublet')
    write.table(x = doubletNoSummary, file = file.path(resDir, sprintf("%s_doublets_no_summary.txt", filename)), quote = F, sep = '\t', row.names = F, col.names = T)
    ## -
    save(olDoubletCells, file = file.path(resDir, sprintf("%s_OL_doublet_cells_name.Rdata", filename) ) )
    print('END Step 4: Doublets identification results summary')
    print('=========')
    print(sprintf("COMPLETE doublets identification for sample '%s'.", as.character(cellrangerResList[[x]])))
    print('---END---END---END---')
    ## ---
  }
  return(resDir)
  ## ---
}
##----------------------------------------------------------------------------------------
## Two minor fns to return mitochodrial content search pattern
## and Main_Doublet_Decon() species based on genome input used in findDoublets()
mtPatten            <- function(genomeSpecies) {
  if (genomeSpecies == 'human') return('^MT-')
  if (genomeSpecies == 'mouse') return('^mt-')
  if (genomeSpecies == 'human_mouse') return('^MT-|mt-')
}
## ---
doubletDeconSpecies <- function(genomeSpecies) {
  if (genomeSpecies == 'human') return('hsa')
  if (genomeSpecies == 'mouse') return('mmu')
  if (genomeSpecies == 'rat') return('rno')
}
##----------------------------------------------------------------------------------------
