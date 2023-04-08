## Developed by Yan Li, June, 2021                                                        ##
##----------------------------------------------------------------------------------------##
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(gridExtra)
## -------------------------------------------------------------------------------------- ##
#' getExpCondClusterMarkers() Function
#' @details
#' This function is used to identify positively expressed cluster marker genes for specified cell clusters with respect to different samples attributes levels specified in the metadata table column.
#'
#' @param resDir specify an exiting full path of directory, where results will be saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved with above rdsDir option.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, allExpConds, or expCond1/2/....
#' @param expCondCheckFname suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param cellcluster specify cell clusters to be extracted for the cluster markers identification.
#' @param pAdjValCutoff adjusted p-value cutoff for significant positively expressed cluster markers, by default = 0.05.
#' @param topNo specify the top number of significantly over expressed cluster markers in each identified/annotated clusters presented in heatmap , by default = 10.
#' @param deMethod DE test method with options: 'wilcox', 't', 'negbinom', 'poisson', 'MAST', 'DESeq2', and 'DESeq2.bulk', default = 'wilcox'.
#' @param min.pct only test genes that are detected in this specified minimum fraction of cells in either of these 2 comparison populations, default is 0.1.
#' @param logfc.threshold only test genes that are detected in this specified X-fold difference (log-scale) between these 2 comparison populations cells, default is 0.25 (around 1.19 FC).
#' @param min.cells.group Minimum number of cells in one of the comparison groups, by default 10.
#' @param deseq2bulk.metaCovariateInput when deMethod = 'DESeq2.bulk', use this option to provide the meta data table for covariate correction. This input should be a data frame object with all corresponding sample items in rows and corresponding group information on columns.
#' @param norm.method DESeq2 counts across samples normalization method, options are TMM, UQ or DEseq2.
#' @param run.dispersion whether to run DESeq2 counts across genes dispersion normalization, if norm.method='DEseq2', this option is automatically on.
#' @param debug whether to turn on for debug check, by default FALSE (turned off).
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @importFrom Seurat RenameIdents
#' @importFrom Seurat FindAllMarkers
#' @importFrom utils write.table
#' @importFrom Seurat DoHeatmap
#' @importFrom Seurat NoLegend
#' @importFrom dplyr %>%
#' @importFrom dplyr top_n
#' @importFrom dplyr group_by
#' @importFrom dplyr distinct
#' @importFrom dplyr summarise
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @keywords getExpCondClusterMarkers
#' @examples getExpCondClusterMarkers()
#' @export
#' @return
#' a list item including 2 elements: All positively expressed genes from each experimental condition cell clusters in 'expCondPosMarkers' and all positively significantly expressed genes at FDR corrected p-value of 0.05 (by default) from each experimental condition cell clusters in'expCondSigPosMarkers'
## -------------------------------------------------------------------------------------- ##
getExpCondClusterMarkers <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL,
                                     expCondCheck='allExpConds', expCondCheckFname = NULL,
                                     cellcluster = NULL,
                                     deMethod = 'wilcox',
                                     deseq2bulk.metaCovariateInput = NULL,
                                     min.pct = 0.1, logfc.threshold = 0.25, min.cells.group = 10, pAdjValCutoff = 0.05, topNo = 10,
                                     norm.method = 'TMM', run.dispersion = as.logical(T),
                                     debug = F) {
  options(java.parameters = "-Xmx32000m")
  ###--------------------------------------------------------------------------------------##
  pAdjValCutoff           <- as.numeric(pAdjValCutoff)
  topNo                   <- as.numeric(topNo)
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
    resDir <- paste(resDir, 'clusterMarkerGenes_results_wNewAnnotation', sep = '/')
  } else {
    resDir <- paste(resDir, 'clusterMarkerGenes_results_wOrgClusterAnnotation', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  ## -------------------------------------------------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  }
  # print(table(Idents(seuratObjFinal)))
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
  print(sprintf('Cluster marker genes identification on cell clusters with respect experimental condition %s will be saved at %s', expCondCheckFname, resDir))
  resDir                <- paste(sprintf('%s/%s', resDir, expCondCheckFname ))
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond'
  if (expCondCheck == 'sample') {
    seuratObjFinal                     <- seuratObjFinal
  } else if (expCondCheck == 'allExpConds') {
    seuratObjFinal@meta.data$expCond   <- 'allExpConds'
  } else {
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  ##--------------------------------------------------------------------------------------##
  ## if 'cellcluster' is provided, subset on 'cellcluster'
  clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% clusterLevels ) ) stop('Please provide the corresponding cell clusters in identfied idents.')
    print(sprintf('Subsetting %s specific cell clusters: %s', length(cellcluster), paste(cellcluster, collapse = ',')))
    seuratObjFinal        <- subset(seuratObjFinal, idents = cellcluster )
  }
  ##--------------------------------------------------------------------------------------##
  ## updated expCond factor after above updating
  expCondLevels           <- levels(factor(seuratObjFinal@meta.data$expCond))
  ## ---
  expCondPosMarkers       <- list()
  expCondSigPosMarkers    <- list()
  if (Seurat::DefaultAssay(seuratObjFinal)!='integrated') {
    stop("Please conduct data integration before using this function.")
  }
  print('-=-=-=--=-=-=-=-=-=-')
  print("Conducting cluster markers analysis on below experimental conditions respectively.")
  print(table(seuratObjFinal@meta.data$expCond))
  ## -------------------------------------------------------------------------------------
  for (l in 1:length(expCondLevels)) {
    ## ---
    systime1              <- Sys.time()
    ## update 'resDir' to ceate dir under 'result_wNewAnnotation' or 'results_wOrgClusterAnnotation'
    expCondLevel          <- expCondLevels[l]
    ## subsetting different experimental conditions
    seuratObjFinalexpCond <- subset(seuratObjFinal, expCond == expCondLevel)
    ## identify cluster positively expressed markers
    print(sprintf('Start %s: finding positive regulated cluster marker genes for experimental condition: %s', l, expCondLevel))
    ## identify significant positively expressed cluster markers
    if (deMethod == 'DESeq2.bulk') {
      Seurat::DefaultAssay(seuratObjFinalexpCond) <- 'RNA'
      seuratObjFinalexpCond@meta.data$deseq2bulk = seuratObjFinalexpCond@meta.data$orig.ident
      clusterLevels2      <- levels(Seurat::Idents(seuratObjFinalexpCond))
      seuratObjFinalexpCond$expCond  <- Seurat::Idents(seuratObjFinalexpCond)
      if (debug) print(sprintf("%s cell types used for DEseq2 bulk analysis are %s.", length(clusterLevels2), paste(as.character(clusterLevels2), collapse = ', ')))
      for (c in 1:length(clusterLevels2)) {
        if (debug) print(sprintf("%s.Processing %s cell clusters for DEseq2 bulk analysis", c, as.character(clusterLevels2[c])))
        deseq2bulkRes     <- runBulkDEseq2(object = seuratObjFinalexpCond, ident.1 = clusterLevels2[c], min.pct = min.pct, metaCovariateInput = deseq2bulk.metaCovariateInput, debug = debug, min.cells.group = min.cells.group, norm.method = norm.method, run.dispersion = run.dispersion)
        deseq2bulkRes$cluster = as.character(clusterLevels2[c])
        deseq2bulkResPos  <- deseq2bulkRes %>% dplyr::filter(log2FoldChange > 0)
        if (c ==1) {
          allPosMarkers   = deseq2bulkResPos
        } else {
          allPosMarkers   = rbind(allPosMarkers, deseq2bulkResPos)
        }
        if (debug) print(sprintf("%END %s cell clusters for DEseq2 bulk analysis", c, as.character(clusterLevels2[c])))
        if (debug) print(head(deseq2bulkRes))
        if (debug) print(sprintf("deseq2bulkRes dimension [%s, %s].", dim(deseq2bulkRes)[1], dim(deseq2bulkRes)[2]))
        if (debug) print(sprintf("allPosMarkers dimension [%s, %s].", dim(allPosMarkers)[1], dim(allPosMarkers)[2]))
        if (debug) print('0000000000000000000000000')
      }
    } else {
      # allPosMarkers       <- FindAllMarkers(seuratObjFinalexpCond, assay = 'integrated', slot = "scale.data", only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold, min.cells.group = min.cells.group, test.use = deMethod)
      allPosMarkers       <- FindAllMarkers(seuratObjFinalexpCond, assay = 'integrated', slot = "data", only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold, min.cells.group = min.cells.group, test.use = deMethod)
    }
    ## ---
    if (deMethod == 'DESeq2.bulk'){
      allPosMarkersAdjSig   <- allPosMarkers %>% tibble::rownames_to_column(var = "gene")%>% dplyr::filter(padj <= pAdjValCutoff) %>% dplyr::mutate(FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -2^(-log2FoldChange)) ) %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(FC), .by_group = TRUE) %>% as.data.frame()
      seuratObjFinalexpCond <- Seurat::ScaleData(object = seuratObjFinalexpCond, do.scale = T, do.center = T)
    } else {
      # allPosMarkersAdjSig   <- allPosMarkers %>% dplyr::filter(p_val_adj <= pAdjValCutoff) %>% dplyr::mutate(perDiff = pct.1-pct.2) %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(avg_diff), .by_group = TRUE)%>% as.data.frame() ## run on scale.data
      allPosMarkersAdjSig   <- allPosMarkers %>% dplyr::filter(p_val_adj <= pAdjValCutoff) %>% dplyr::mutate(perDiff = pct.1-pct.2) %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)%>% as.data.frame() ## run on data
    }
    ## ---
    write.table(x = allPosMarkers, file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_resFull.txt', resDir, expCondCheckFname, expCondLevel)), quote = F, sep = '\t', row.names = T, col.names = NA)
    if (dim(allPosMarkersAdjSig)[1] > 0) write.table(x = allPosMarkersAdjSig, file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig.txt', resDir, expCondCheckFname, expCondLevel)), quote = F, sep = '\t', row.names = F, col.names = T)
    print(sprintf('A total of %s positively expressed genes identified for experimental condition %s at %s, among them %s are significant up expressed at adjusted p-value significant level of %s', dim(allPosMarkers)[1], expCondCheckFname, expCondLevel, dim(allPosMarkersAdjSig)[1], pAdjValCutoff ))
    ## summarize the no. of significant positively expressed cluster markers
    if (dim(allPosMarkersAdjSig)[1] > 0) {
      allPosMarkersAdjSigNo <- allPosMarkersAdjSig %>% dplyr::group_by(cluster) %>% dplyr::distinct() %>% dplyr::summarise('geneNo' = dplyr::n()) %>% as.data.frame()
      write.table(x = allPosMarkersAdjSigNo, file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig_NoSummary.txt', resDir, expCondCheckFname, expCondLevel)), quote = F, sep = '\t', row.names = F, col.names = T)
    }
    ## -
    expCondPosMarkers[[l]]     <- allPosMarkers
    expCondSigPosMarkers[[l]]  <- allPosMarkersAdjSig
    systime2               <- Sys.time()
    print(sprintf('END %s: finding positive regulated cluster marker genes for experimental condition in %s with computation time: %s %s.', l, expCondLevel, round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
    if (any(table(allPosMarkersAdjSig$cluster) >= topNo)) {
      print('Start: Step 6 making cluster marker genes heatmap plot')
      ## all significant cluster markers heatmap
      if (deMethod == 'DESeq2.bulk') {
        topMarkers               <- allPosMarkersAdjSig %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(FC), .by_group = TRUE) %>% dplyr::top_n(n = topNo, FC) %>% as.data.frame()
        cluterAllsigMarkerheatmap  <- DoHeatmap(seuratObjFinalexpCond, assay = 'RNA', slot = "scale.data", features = topMarkers$gene)
      } else {
        topMarkers               <- allPosMarkersAdjSig %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>% dplyr::top_n(n = topNo, avg_log2FC) %>% as.data.frame()
        cluterAllsigMarkerheatmap  <- DoHeatmap(seuratObjFinalexpCond, assay = 'integrated', slot = "scale.data", features = topMarkers$gene)
      }
      if (topNo<=10) {
        pdf(file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig_heatmap.pdf', resDir, expCondCheckFname, expCondLevel )), width = 25, height = 8)
        print(cluterAllsigMarkerheatmap)
        dev.off()
      } else if (topNo > 10 & topNo <=100) {
        pdf(file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig_heatmap.pdf', resDir, expCondCheckFname, expCondLevel )), width = 25, height = 10)
        print(cluterAllsigMarkerheatmap)
        dev.off()
      } else {
        pdf(file = file.path(sprintf('%s/expCond_%s_%s_allCluster_pos_markers_UpSig_heatmap.pdf', resDir, expCondCheckFname, expCondLevel )), width = 25, height = 25)
        print(cluterAllsigMarkerheatmap)
        dev.off()
      }

    }
    print('END: Step 6 making cluster marker genes heatmap plot')
    print('********************')
    ## ---
  }
  names(expCondPosMarkers)    <- expCondLevels
  names(expCondSigPosMarkers) <- expCondLevels
  ##--------------------------------------------------------------------------------------##
  return(list('fullRes'=expCondPosMarkers, 'sigUp' = expCondSigPosMarkers))
  ##--------------------------------------------------------------------------------------##
}
## -------------------------------------------------------------------------------------- ##
