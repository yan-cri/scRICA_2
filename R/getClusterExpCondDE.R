#' getClusterExpCondDe() Function
#' @details
#' This function is used to identify DEGs in each identified/annotated cell cluster based on the experimental conditions from metadata table column 'sample', 'expCond1', or 'expCond2'
#'
#' @param resDir specify an exiting full path of directory, where results will be saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved with above rdsDir option.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param expCondCheckFname suffix of the directory/folder and file name of the dot plot to be saved, if not defined, the same as the 'expCondCheck' option.
#' @param cellcluster specify the specific cell cluster names for DE analysis, if not specified, all cell clusters will be performed .
#' @param compGroup specify either 1 or 2 group names seperated by '/' for comparisons, these group names is corresponding to experimental condition factor levels. If only 1 group named is specified, it makes comparison between this specified group with all others.
#' @param deMethod DE test method with options: 'wilcox', 't', 'negbinom', 'poisson', 'MAST', 'DESeq2', 'DESeq2.bulk', by default = 'wilcox'.
#' @param pAdjValCutoff adjusted p-value cutoff for significant DEGs detection, by default = 0.05.
#' @param topNo specify the top number of up and down significant DEGs in each identified/annotated cell clusters, by default = 10.
#' @param min.pct only test genes that are detected in this specified minimum fraction of cells in either of these 2 comparison populations, default is 0.1.
#' @param logfc.threshold only test genes that are detected in this specified X-fold difference (log-scale) between these 2 comparison populations cells, default is 0.25 (around 1.19 FC).
#' @param min.cells.group Minimum number of cells in one of the comparison groups.
#' @param deseq2bulk.metaCovariateInput when deMethod = 'DESeq2.bulk', use this option to provide the meta data table for covariate correction. This input should be a data frame object with all corresponding sample items in rows and corresponding group information on columns.
#' @param covariateVarName this option is only used when deMethod = 'MIST' with this specified column for covariate correction.
#' @param norm.method DESeq2 counts across samples normalization method, options are TMM, UQ or DEseq2.
#' @param run.dispersion whether to run DESeq2 counts across genes dispersion normalization, if norm.method='DEseq2', this option is automatically on.
#' @param debug whether to turn on for debug check, by default FALSE (turned off).
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @importFrom Seurat RenameIdents
#' @importFrom Seurat FindMarkers
#' @importFrom utils write.table
#' @importFrom xlsx write.xlsx
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 results
#' @importFrom DESeq2 DESeq
#'
#' @keywords FindMarkers
#' @examples getClusterExpCondDe()
#' @export
#' @return
#' a list item including 3 elements:
#' 1: 'clusterDeMarkers'
#' 2. 'clusterDeResSummary'
#' 3. 'clusterTopDeMarkers'
## ---------------------------------------------------------------------------------------
getClusterExpCondDe <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, expCondCheck='sample',
                                expCondCheckFname = NULL, cellcluster = NULL, compGroup, deMethod = 'wilcox',
                                covariateVarName = NULL, deseq2bulk.metaCovariateInput = NULL, covariateVarLevels = NULL, norm.method = 'UQ', run.dispersion = as.logical(T),
                                min.cells.group = 10, min.pct = 0.1, logfc.threshold = 0.25, pAdjValCutoff = 0.05, topNo = 10, debug = F, outputExcel = as.logical(T)) {
  options(java.parameters = "-Xmx128g")
  ## ---
  if (missing(compGroup)) stop("Please provide option 'compGroup' to specify which 2 groups in your updated experimental condition levels with 'expCondCheck' for comparision.")
  pAdjValCutoff           <- as.numeric(pAdjValCutoff)
  topNo                   <- as.numeric(topNo)
  deMethod                <- as.character(deMethod)
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) print("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  if (!is.null(covariateVarName)) {
    if (!deMethod %in% c('LR', 'negbinom', 'poisson','MAST', 'DESeq2')) stop("'covariateVarName' cannot be applied for the specified 'deMethod'.")
  }
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
    resDir <- paste(resDir, 'results_wNewAnnotation_DEGs', sep = '/')
  } else {
    resDir <- paste(resDir, 'results_wOrgClusterAnnotation_DEGs', sep = '/')
  }
  if (!dir.exists(resDir)) dir.create(resDir)
  print(sprintf('DEG identification on cell clusters will be saved at %s', resDir))
  ## -------------------------------------------------------------------------------------
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
  }
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
  ## update resDir
  resDir                  <- paste(resDir, expCondCheckFname, sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond' and add 'covariateVarName'
  if (expCondCheck == 'sample') {
    seuratObjFinal                     <- seuratObjFinal
    if (covariateVarName == 'sample') stop("DE tests are conducted on 'sample', no need to regress on this variable again with 'covariateVarName'.")
    if (deMethod == 'DESeq2.bulk') stop("DEseq2 bulk RNA tests are combining results on the sample level, no DE analysis can be conducted on the same level.")
  } else {
    if (deMethod == 'MAST') {
      ## covariateVar used for MAST latent.varible option
      if (!is.null(covariateVarName)) {
        if (covariateVarName == 'sample') {
          seuratObjFinal@meta.data$covariateVar = seuratObjFinal@meta.data$orig.ident
        } else {
          seuratObjFinal@meta.data$covariateVar = seuratObjFinal@meta.data[, grep(as.character(covariateVarName), colnames(seuratObjFinal@meta.data))]
        }
      }
      ## ---
    }
    if (deMethod == 'DESeq2.bulk') {
      seuratObjFinal@meta.data$deseq2bulk = seuratObjFinal@meta.data$orig.ident
    }
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  ##--------------------------------------------------------------------------------------##
  Seurat::DefaultAssay(seuratObjFinal)   <- "RNA"
  ##--------------------------------------------------------------------------------------##
  ## if 'cellcluster' provided, only certain cell clusters to conduct DE analysis.
  clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% clusterLevels ) ) stop('Please provide the corresponding cell clusters in identfied idents.')
    print(sprintf("Start step 1: identifing DEGs for selected clusters (%s) based on experimental conditions in rds metadata table column '%s'.", paste(unlist(cellcluster), collapse = ', '), expCondCheck ))
    noClusters                             <- unique(cellcluster)
  } else {
    print(sprintf("Start step 1: identifing DEGs for all identified clusters based on experimental conditions in rds metadata table column '%s'.", expCondCheck))
    noClusters                             <- levels(Seurat::Idents(seuratObjFinal))
  }
  ## -----
  print(sprintf("A total of %s clusters to process", length(noClusters)))
  ## -
  clusterDeMarkers                       <- list()
  normCounts <- list()
  orgCounts  <- list()
  metaTabs   <- list()
  clusterDeResSummary                    <- matrix(data = NA, nrow = length(noClusters), ncol = 4)
  clusterTopDeMarkers                    <- list()
  clusterTopDeMarkersUp                  <- NA
  clusterTopDeMarkersUpCluster           <- NA
  clusterTopDeMarkersDown                <- NA
  clusterTopDeMarkersDownCluster         <- NA
  clusterDeMarkersFname                  <- paste(sprintf('%s/expCondCompDeMarkers_%s', resDir, gsub('/', '-', compGroup) ))
  ## -
  if (is.null(cellcluster)) {
    if (newAnnotation) {
      resFname1                  <- paste(clusterDeMarkersFname, '_full_wNewAnnotation_allClusters', sep = '')
      resFname2                  <- paste(clusterDeMarkersFname, '_adjSig_wNewAnnotation_allClusters', sep = '')
      resFname3                  <- paste(clusterDeMarkersFname, '_adjSig_up_wNewAnnotation_allClusters', sep = '')
      resFname4                  <- paste(clusterDeMarkersFname, '_adjSig_down_wNewAnnotation_allClusters', sep = '')
    } else {
      resFname1                  <- paste(clusterDeMarkersFname, '_full_allClusters', sep = '')
      resFname2                  <- paste(clusterDeMarkersFname, '_adjSig_allClusters', sep = '')
      resFname3                  <- paste(clusterDeMarkersFname, '_adjSig_up_allClusters', sep = '')
      resFname4                  <- paste(clusterDeMarkersFname, '_adjSig_down_allClusters', sep = '')
    }
  } else {
    if (newAnnotation) {
      resFname1                  <- paste(clusterDeMarkersFname, sprintf("_full_wNewAnnotation_%sSelClusters", length(noClusters)), sep = '')
      resFname2                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_wNewAnnotation_%sSelClusters", length(noClusters)), sep = '')
      resFname3                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_up_wNewAnnotation_%sSelClusters", length(noClusters)), sep = '')
      resFname4                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_down_wNewAnnotation_%sSelClusters", length(noClusters)), sep = '')
    } else {
      resFname1                  <- paste(clusterDeMarkersFname, sprintf("_full_%sSelClusters", length(noClusters)), sep = '')
      resFname2                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_%sSelClusters", length(noClusters)), sep = '')
      resFname3                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_up_%sSelClusters", length(noClusters)), sep = '')
      resFname4                  <- paste(clusterDeMarkersFname, sprintf("_adjSig_down_%sSelClusters", length(noClusters)), sep = '')
    }
  }
  ## -
  if (debug) print(sprintf("clusters to be processed are %s", paste(noClusters, collapse = ', ')))
  for (i in 1:length(noClusters)) {
    # for (i in c(1,3,4)) {
    print(sprintf('Start 1.%s processing cluster %s for DE markers identification', i, noClusters[i]))
    seuratObjFinalSubet                  <- subset(seuratObjFinal, idents = as.character(noClusters[i]) )
    print("-=-=-=")
    print("experimental conditions:")
    print(table(seuratObjFinalSubet$expCond))
    print("-=-=-=")
    ## -------------------- ##
    if (deMethod == 'MAST' & !is.null(covariateVarName)) {
      if (!is.null(covariateVarLevels)) {
        if (debug) {
          print("*********")
          print(table(seuratObjFinalSubet@meta.data$covariateVar))
          print(sprintf("input 'covariateVarLevels': %s", paste(covariateVarLevels, collapse = ', ')))
          print("*********")
        }
        if (length(unique(seuratObjFinalSubet@meta.data$covariateVar))<length(covariateVarLevels)) {
          covariateVarLevels <- covariateVarLevels[which(!is.na(match(covariateVarLevels, unique(seuratObjFinalSubet@meta.data$covariateVar))))]
          if (debug) print(sprintf("covariateVarLevels is updated into %s. ", paste(covariateVarLevels, collapse = ', ')))
        }
        seuratObjFinalSubet@meta.data$covariateVar = factor(seuratObjFinalSubet@meta.data$covariateVar, levels = covariateVarLevels)
        # seuratObjFinalSubet@meta.data$covariateVar = relevel(x =  seuratObjFinalSubet@meta.data$covariateVar, ref = covariateVarLevels[length(covariateVarLevels)] )
      } else {
        covariateVarLevels <- names(sort(table(seuratObjFinalSubet@meta.data$covariateVar), decreasing = F))
        if (debug) {
          print("*********")
          print(table(seuratObjFinalSubet@meta.data$covariateVar))
          print(sprintf("sorted levels: %s", paste(covariateVarLevels, collapse = ', ')))
          print("*********")
        }
        seuratObjFinalSubet@meta.data$covariateVar = factor(seuratObjFinalSubet@meta.data$covariateVar, levels = covariateVarLevels)
        # seuratObjFinalSubet@meta.data$covariateVar = relevel(x =  seuratObjFinalSubet@meta.data$covariateVar, ref = covariateVarLevels[length(covariateVarLevels)] )
      }
    }
    ## -------------------- ##
    compGroupName                      <- strsplit(compGroup, split = '/')[[1]]
    if ( all(compGroupName %in% names(table(seuratObjFinalSubet$expCond))) ) {
      ## ---
      if ( length(compGroupName)>2 ) {
        stop("Too many comparision groups provided in 'compGroup', only 2 names shoulb be provided ")
      } else {
        if (length(compGroupName) == 1) {
          if (all(table(seuratObjFinalSubet$expCond) >= min.cells.group) ) {
            if (deMethod=='DESeq2.bulk') {
              deMarkersRes                 <- runBulkDEseq2(object = seuratObjFinalSubet, ident.1 =  compGroupName, min.pct = min.pct, metaCovariateInput = deseq2bulk.metaCovariateInput, debug = debug, min.cells.group = min.cells.group, norm.method = norm.method, run.dispersion = run.dispersion)
              deMarkers                    <- deMarkersRes$deres
              orgCount                     <- deMarkersRes$orgCount
              normCount                    <- deMarkersRes$normCount
              metaTab                      <- deMarkersRes$metaTab
            } else {
              if (is.null(covariateVarName)) {
                deMarkers                  <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName, group.by = 'expCond', test.use = deMethod, min.cells.group = min.cells.group, logfc.threshold = logfc.threshold, min.pct =  min.pct)
              } else {
                if (debug) {
                  print('--------------')
                  print(table(seuratObjFinalSubet@meta.data$covariateVar))
                  print('--------------')
                }
                deMarkers                  <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName, group.by = 'expCond', test.use = deMethod, min.cells.group = min.cells.group, logfc.threshold = logfc.threshold, min.pct =  min.pct, latent.vars = 'covariateVar')
              }
            }
          } else {
            print(sprintf('cluster %s has less than 3 cells in one experimental condiction, no DE markers can be identfied for this cluster', noClusters[i]))
          }
        } else if (length(compGroupName) == 2) {
          if (all(table(seuratObjFinalSubet$expCond)[match(compGroupName, names(table(seuratObjFinalSubet$expCond)))] >= min.cells.group) ) {
            if (deMethod=='DESeq2.bulk') {
              deMarkersRes                 <- runBulkDEseq2(object = seuratObjFinalSubet, ident.1 = compGroupName[1], ident.2 = compGroupName[2], min.pct = min.pct, metaCovariateInput = deseq2bulk.metaCovariateInput, debug = debug, min.cells.group = min.cells.group, norm.method = norm.method, run.dispersion = run.dispersion)
              deMarkers                    <- deMarkersRes$deres
              orgCount                     <- deMarkersRes$orgCount
              normCount                    <- deMarkersRes$normCount
              metaTab                      <- deMarkersRes$metaTab
            } else {
              if (is.null(covariateVarName)) {
                deMarkers                  <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName[1], ident.2 = compGroupName[2], group.by = 'expCond', test.use = deMethod, min.cells.group = min.cells.group, logfc.threshold = logfc.threshold, min.pct =  min.pct)
              } else {
                if (debug) {
                  print('--------------')
                  print(table(seuratObjFinalSubet@meta.data$covariateVar))
                  print('--------------')
                }
                deMarkers                  <- FindMarkers(seuratObjFinalSubet, ident.1 = compGroupName[1], ident.2 = compGroupName[2], group.by = 'expCond', test.use = deMethod, min.cells.group = min.cells.group, logfc.threshold = logfc.threshold, min.pct =  min.pct, latent.vars = 'covariateVar')
              }
            }
          } else {
            print(sprintf('cluster %s has less than 3 cells in one experimental condiction, no DE markers can be identfied for this cluster', noClusters[i]))
          }
        }
        if(debug) print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        if(debug) print('111111111111111')
        clusterDeMarkers[[i]]          <- deMarkers
        if (deMethod!="DESeq2.bulk") {
          # print(sprintf('Maximum p_value is %s, Maximum adjusted p_value is %s', round(max(deMarkers$pvalue), digits = 4), round(max(deMarkers$padj, na.rm = T), digits = 4)))
          deMarkersAdjSig                <- deMarkers %>% dplyr::filter(p_val_adj <= pAdjValCutoff) %>% dplyr::filter(abs(avg_log2FC) >= logfc.threshold) %>% dplyr::mutate(perDiff = pct.1-pct.2) %>% dplyr::mutate(FC = ifelse(avg_log2FC>0, 2^avg_log2FC, -2^(-avg_log2FC)) )
          deMarkersAdjSigUp              <- deMarkersAdjSig %>% dplyr::filter(avg_log2FC > 0) %>% dplyr::arrange(desc(FC))
          deMarkersAdjSigDown            <- deMarkersAdjSig %>% dplyr::filter(avg_log2FC < 0) %>% dplyr::arrange(FC)
        } else {
          if(debug) print('22222222222222')
          normCounts[[i]] <- normCount
          orgCounts[[i]]  <- orgCount
          metaTabs[[i]]   <- metaTab
          # print(sprintf('Maximum p_value is %s, Maximum adjusted p_value is %s', round(max(deMarkers$p_val), digits = 4), round(max(deMarkers$p_val_adj), digits = 4)))
          deMarkersAdjSig                <- deMarkers %>% dplyr::filter(padj <= pAdjValCutoff) %>% dplyr::filter(abs(log2FoldChange) > logfc.threshold) %>% dplyr::mutate(FC = ifelse(log2FoldChange>0, 2^log2FoldChange, -2^(-log2FoldChange)) )
          deMarkersAdjSigUp              <- deMarkersAdjSig %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::arrange(desc(FC))
          deMarkersAdjSigDown            <- deMarkersAdjSig %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::arrange(FC)
        }
        if(debug) print('333333333333')
        if (dim(deMarkersAdjSigUp)[1] > topNo) {
          topNo1 = topNo
        } else {
          if (dim(deMarkersAdjSigUp)[1]!=0) {
            topNo1 = dim(deMarkersAdjSigUp)[1]
          } else {
            topNo1 = 0
          }
        }
        if(debug) print('44444444444444')
        if (dim(deMarkersAdjSigDown)[1] > topNo) {
          topNo2 = topNo
        } else {
          if (dim(deMarkersAdjSigDown)[1]!=0) {
            topNo2 = dim(deMarkersAdjSigDown)[1]
          } else {
            topNo2 = 0
          }
        }
        if(debug) print('55555555555555')
        if (topNo1!=0 & topNo2!=0) {
          clusterTopDeMarkers[[i]]       <- list('up' = rownames(deMarkersAdjSigUp)[1:topNo1],
                                                 'down' = rownames(deMarkersAdjSigDown)[1:topNo2])
        } else if (topNo1!=0 & topNo2==0) {
          clusterTopDeMarkers[[i]]       <- list('up' = rownames(deMarkersAdjSigUp)[1:topNo1])
        } else if (topNo1==0 & topNo2!=0) {
          clusterTopDeMarkers[[i]]       <- list('down' = rownames(deMarkersAdjSigDown)[1:topNo2])
        } else {
          clusterTopDeMarkers[[i]]       <- NA
        }

        if(debug) print('666666666666666')
        if (dim(deMarkersAdjSigUp)[1]!=0) {
          clusterTopDeMarkersUp          <- c(clusterTopDeMarkersUp, rownames(deMarkersAdjSigUp)[1:topNo1])
          clusterTopDeMarkersUpCluster   <- c(clusterTopDeMarkersUpCluster, rep(noClusters[i], length(rownames(deMarkersAdjSigUp)[1:topNo1])) )
        }
        if (dim(deMarkersAdjSigDown)[1]!=0){
          clusterTopDeMarkersDown        <- c(clusterTopDeMarkersDown, rownames(deMarkersAdjSigDown)[1:topNo2])
          clusterTopDeMarkersDownCluster <- c(clusterTopDeMarkersDownCluster, rep(noClusters[i], length(rownames(deMarkersAdjSigDown)[1:topNo2])))
        }
        clusterDeResSummary[i,]        <- c(dim(deMarkers)[1], dim(deMarkersAdjSig)[1], dim(deMarkersAdjSigUp)[1], dim(deMarkersAdjSigDown)[1] )
        print(sprintf('out of %s DE markers, %s are significantly DE at adjusted p-value of %s', dim(deMarkers)[1], dim(deMarkersAdjSig)[1], pAdjValCutoff))
        print(sprintf('out of %s significantly DE markers, %s are positively expressed, %s are negative expressed.', dim(deMarkersAdjSig)[1], dim(deMarkersAdjSigUp)[1], dim(deMarkersAdjSigDown)[1] ))
        ## -
        if(debug) print('777777777777777')
        if (outputExcel) {
          if (i ==1) {
            if (dim(deMarkers)[1] > 0) write.xlsx(x = deMarkers, file = sprintf('%s.xlsx', resFname1), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSig)[1] > 0) write.xlsx(x = deMarkersAdjSig, file = sprintf('%s.xlsx', resFname2), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSigUp)[1] > 0) write.xlsx(x = deMarkersAdjSigUp, file = sprintf('%s.xlsx', resFname3), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = F )
            if (dim(deMarkersAdjSigDown)[1] > 0) write.xlsx(x = deMarkersAdjSigDown, file = sprintf('%s.xlsx', resFname4), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = F )
          } else {
            if (dim(deMarkers)[1] > 0) write.xlsx(x = deMarkers, file = sprintf('%s.xlsx', resFname1), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSig)[1] > 0)write.xlsx(x = deMarkersAdjSig, file = sprintf('%s.xlsx', resFname2), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSigUp)[1] > 0) write.xlsx(x = deMarkersAdjSigUp, file = sprintf('%s.xlsx', resFname3), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = T )
            if (dim(deMarkersAdjSigDown)[1] > 0) write.xlsx(x = deMarkersAdjSigDown, file = sprintf('%s.xlsx', resFname4), sheetName = paste('cluster', gsub('/|:|[ ]|-', '', noClusters[i]), sep = '_'), row.names = T, append = T )
          }
        } else {
          if (dim(deMarkers)[1] > 0) write.table(x = deMarkers, file = sprintf('%s_cluster%s.txt', resFname1, gsub('/|:|[ ]|-', '', noClusters[i]) ), row.names = T, quote = F, col.names = NA, sep = '\t')
          if (dim(deMarkersAdjSig)[1] > 0) write.table(x = deMarkersAdjSig, file = sprintf('%s_cluster%s.txt', resFname2, gsub('/|:|[ ]|-', '', noClusters[i])), row.names = T, quote = F, col.names = NA, sep = '\t')
          if (dim(deMarkersAdjSigUp)[1] > 0) write.table(x = deMarkersAdjSigUp, file = sprintf('%s_cluster%s.txt', resFname3, gsub('/|:|[ ]|-', '', noClusters[i])), row.names = T, quote = F, col.names = NA, sep = '\t')
          if (dim(deMarkersAdjSigDown)[1] > 0) write.table(x = deMarkersAdjSigDown, file = sprintf('%s_cluster%s.txt', resFname4, gsub('/|:|[ ]|-', '', noClusters[i])), row.names = T, quote = F, col.names = NA, sep = '\t')
        }
        if(debug) print('8888888888888888')
      }
      ## ---
    } else {
      stop("Provided comparision group names in 'compGroup' does not match experimental conditions ")
    }
    print(sprintf('End 1.%s processing cluster %s for DE markers identification', i, noClusters[i]))
    print('=========')
  }
  ## -
  clusterTopDeMarkersUpComb          <- data.frame('Gene' = clusterTopDeMarkersUp[-1], 'geneType'= clusterTopDeMarkersUpCluster[-1])
  clusterTopDeMarkersDownComb        <- data.frame('Gene' = clusterTopDeMarkersDown[-1], 'geneType' = clusterTopDeMarkersDownCluster[-1])
  if(debug) print(head(clusterTopDeMarkersUpComb))
  if (dim(clusterTopDeMarkersUpComb[!is.na(clusterTopDeMarkersUpComb$Gene),])[1]!=0) {
    write.xlsx( x = clusterTopDeMarkersUpComb[!is.na(clusterTopDeMarkersUpComb$Gene),], file = sprintf('%s_top%s_upDe.xlsx', clusterDeMarkersFname, topNo1), sheetName = 'up', row.names = F, col.names = T, append = F)
  }
  if (dim(clusterTopDeMarkersDownComb[!is.na(clusterTopDeMarkersDownComb$Gene),])[1] !=0) {
    write.xlsx( x = clusterTopDeMarkersDownComb[!is.na(clusterTopDeMarkersDownComb$Gene),], file = sprintf('%s_top%s_downDe.xlsx', clusterDeMarkersFname, topNo2), sheetName = 'down', row.names = F, col.names = T, append = F)
  }
  ## -
  rownames(clusterDeResSummary)      <- noClusters
  colnames(clusterDeResSummary)      <- c('Total', 'DE', 'Up', 'Down')
  clusterDeResSummaryFname           <- paste(clusterDeMarkersFname, '_NoDEmarkers_summary.xlsx', sep = '')
  write.xlsx( x = clusterDeResSummary, file = clusterDeResSummaryFname, sheetName = 'No. summary', row.names = T, col.names = T, append = F)
  ## ---
  names(clusterDeMarkers)            <- noClusters
  names(clusterTopDeMarkers)         <- noClusters
  if (deMethod=='DESeq2.bulk') {
    names(normCounts)                  <- noClusters
    names(orgCounts)                   <- noClusters
    names(metaTabs)                    <- noClusters
  }
  print(sprintf('End step 1: identifing DEGs for each identified clusters based on sample name experimental conditions'))
  print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
  if (deMethod=='DESeq2.bulk') {
    return(list('orgCounts' = orgCounts, 'normCounts' = normCounts, 'metaTabs' = metaTabs, 'clusterDeMarkers' = clusterDeMarkers, 'clusterDeResSummary' = clusterDeResSummary, 'clusterTopDeMarkers' = clusterTopDeMarkers ))
  }else {
    return(list('clusterDeMarkers' = clusterDeMarkers, 'clusterDeResSummary' = clusterDeResSummary, 'clusterTopDeMarkers' = clusterTopDeMarkers ))
  }
  ## -------------------------------------------------------------------------------------
}

## ----------------- ##
IdentsToCells <- function( object, group.by, ident.1, ident.2 = NULL, cellnames.use) {
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  }
  Seurat::Idents(object = object) <- group.by
  ident.1 <- Seurat::WhichCells(object = object, idents = ident.1)
  if (is.null(ident.2)) {
    ident.2 <- setdiff(x = cellnames.use, y = ident.1)
  } else {
    ident.2 <- Seurat::WhichCells(object = object, idents = ident.2)
  }
  return(list(group1 = ident.1, group2 = ident.2))
}
## ----------------- ##
runBulkDEseq2 <- function(object, ident.1, ident.2=NULL, min.pct, metaCovariateInput = NULL, debug = F, min.cells.group, norm.method = 'TMM', run.dispersion = as.logical(T)) {
  cellGroup = 'expCond' ##updated with 'expCondCheck', compGroupName should be levels shown in 'expCond'
  min.cells.group = as.numeric(min.cells.group)
  if (debug) print(table(object@meta.data[,match(cellGroup, colnames(object@meta.data))]))
  if (is.null(ident.2)) {
    print(sprintf("group1 = %s, group2 = NULL", ident.1))
  } else {
    print(sprintf("group1 = %s, group2 = %s", ident.1, ident.2))
  }
  ## select cells corresponding ident.1/2
  cells      <- IdentsToCells(object = object, group.by = cellGroup, ident.1 = ident.1, ident.2 = ident.2,  cellnames.use = colnames(object) )
  if (debug) print(str(cells))
  # ## calculate normalization factors with respect to all selected cells, either full sets or subsets.
  # cells.meta <- object@meta.data[match(unlist(cells), rownames(object@meta.data)),]
  # cells.no   <- count(cells.meta, deseq2bulk) %>% as.data.frame()
  # cells.no$norm.factor<- cells.no$n/median(cells.no$n)
  ## calculate FC results
  fc.results <- Seurat::FoldChange(object = object, ident.1 = ident.1, ident.2 = ident.2,  group.by = cellGroup )
  if (debug) {
    print(sprintf("Initially, %s genes are used for FC calculation.", dim(fc.results)[1]))
  }
  ## filter based on min.pct option
  fc.results$pmax <- pmax(fc.results$pct.1, fc.results$pct.2)
  fc.results      <- fc.results[fc.results$pmax>min.pct,]
  genes4de        <- rownames(fc.results)
  print(sprintf("After low expression removal with min.pct = %s, %s genes will be used for DE analysis.", min.pct, length(genes4de)))
  ## aggregating counts based on deseq2bulk inherited from 'orig.ident', and meanwhile prepare corresponding metatab.
  countInputs     <- lapply(1:length(cells), function(x) {
    counts             <- Seurat::FetchData(object = object, vars = genes4de, cells = cells[[x]], slot = "count")
    counts$deseq2bulk  <- Seurat::FetchData(object = object, vars = 'deseq2bulk', cells = cells[[x]])$deseq2bulk
    counts.group       <- counts %>% group_by(deseq2bulk) %>% summarise(across(everything(), list(sum))) %>% as.data.frame()
    colnames(counts.group) <- gsub('_1', '', colnames(counts.group))
    ## normalize the sum by median of number of cells weight
    deseq2bulk.no      <- as.data.frame(table(counts$deseq2bulk))
    ## remove samples with low expression cell number at cut-off of 'min.cells.group'
    if (debug) print("before remval")
    if (debug)  print(deseq2bulk.no)
    deseq2bulk.no        <- deseq2bulk.no[deseq2bulk.no$Freq >= min.cells.group,]
    if (debug) print("after remval")
    if (debug) print(deseq2bulk.no)
    if (debug) print('-=-=-=-=-=-=-=-=-')
    deseq2bulk.no$weight <- deseq2bulk.no$Freq/median(deseq2bulk.no$Freq)
    ## match count table with weight table
    counts.group.match   <-counts.group %>% dplyr::filter(deseq2bulk %in% as.character(deseq2bulk.no$Var1) )
    ##match deseq2bulk.no has the same sample order as shown in counts.group.match
    deseq2bulk.no        <- deseq2bulk.no[match(counts.group.match$deseq2bulk, deseq2bulk.no$Var1),]
    if (debug) print(sprintf("min.cells.group = %s", min.cells.group))
    if (debug) print(sprintf("Before: [row=%s %s], After [row=%s %s]", dim(counts.group)[1], dim(counts.group)[2], dim(counts.group.match)[1], dim(counts.group.match)[2] ) )
    counts.group2            <- round(counts.group.match[,-1]/deseq2bulk.no$weight, digits = 0)
    counts.group2$deseq2bulk <- counts.group.match$deseq2bulk
    counts.group             <- counts.group2
    if (debug) print(dim(counts.group))
    if (debug) print(deseq2bulk.no)
    ## -
    counts.group$group <- paste('group', x, sep = '')
    ## match on metaTabPrep
    metaTabPrep               <- data.frame('bulksamp' = deseq2bulk.no %>% dplyr::pull(Var1))
    if (!is.null(metaCovariateInput)) {
      colnames(metaCovariateInput) <- tolower(colnames(metaCovariateInput))
      if (!any(grepl('^bulksamp$', colnames(metaCovariateInput)))) stop("Please provide corresponding 'metaCovariateInput' with column 'bulksamp' to represent bulk pseudo samples. ")
      metaTabPrep <- dplyr::left_join(x = metaTabPrep, y = metaCovariateInput, by = 'bulksamp')
    }
    metaTabPrep$group <- paste('group', x, sep = '')
    ## ----
    return(list('countTab' = counts.group, 'metaTab' = metaTabPrep))
  } )
  ## combining 2 list of cells counts into a combined DF as 'countTab'
  countInputsComb <- rbind(countInputs[[1]]$countTab, countInputs[[2]]$countTab)
  ## ----------------- ##
  if (debug) print("Complete bulk aggregation.")
  if (debug) print(countInputsComb[,1:4])
  countTab <- countInputsComb %>% dplyr::select(-group) %>% dplyr::select(-deseq2bulk) %>% t()
  colnames(countTab) <- countInputsComb %>% dplyr::pull(deseq2bulk)
  if (debug) print("header of countTab")
  if (debug) print(head(countTab))
  if (debug) print(dim(countTab))
  # if (debug) print(sprintf("countTab.txt in '%s'.", paste(getwd(), 'countTab.txt', sep = '/')))
  # if (debug) write.table(x = countTab, file = paste(getwd(), 'countTab.txt', sep = '/'), col.names = NA, row.names = T, quote = F, sep = '\t')
  if (debug) print('-=-=-=-=-=-=-=-=-=-=-')
  ## ----------------- ##
  ## prepare metaTab, if 'metaCovariateInput' != NULL, add corresponding covariates for model DESeq2 model fitting.
  metaTab         <- rbind(countInputs[[1]]$metaTab, countInputs[[2]]$metaTab)
  metaTab2        <- metaTab %>% dplyr::select(-bulksamp) ## remove sample name for model formular establishment
  ## metaTab2 and metaTab are identical, except 'bulksamp' column is removed in metaTab2 for DESeq2 model formular and DESeq2 readin.
  # metaTab2      <- as.data.frame(unclass(metaTab2),stringsAsFactors=TRUE) ##change all character column into factor
  # rownames(metaTab2) <- as.character(metaTab$bulksamp)
  ## ----------------- ##
  if (debug) print("Full metaTab")
  if (debug) print(metaTab)
  print('-=-=-=-=-=-=-=-=-=-=-')
  print(sprintf("A total of %s samples", dim(metaTab)[1]))
  print(table(metaTab$group))
  print('-=-=-=-=-=-=-=-=-=-=-')
  ## ----------------- ##
  design.formula <- formula(paste0(' ~ ', paste(colnames(metaTab2), collapse = '+')))
  print("Design formular is:")
  print(design.formula)
  print('-=-=-=-=-=-=-=-=-=-=-')
  #
  if (norm.method == 'TMM') {
    dge   <- edgeR::DGEList(counts = countTab)
    if (debug) print("using TMM normalization")
    dge           <- edgeR::calcNormFactors(object = dge, method = 'TMM')
    normCount     <- edgeR::cpm(dge)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(edgeR::cpm(dge), digits = 0), colData=S4Vectors::DataFrame(metaTab2), design = design.formula )
  } else if (norm.method == 'UQ'){
    dge   <- edgeR::DGEList(counts = countTab)
    if (debug) print("using upperquartile normalization")
    dge   <- edgeR::calcNormFactors(object = dge, method = 'upperquartile')
    normCount     <- edgeR::cpm(dge)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(edgeR::cpm(dge), digits = 0), colData=S4Vectors::DataFrame(metaTab2), design = design.formula )
  } else {
    normCount     <- edgeR::cpm(countTab)
    dds   <- DESeq2::DESeqDataSetFromMatrix(countData=countTab, colData=S4Vectors::DataFrame(metaTab2), design = design.formula )
  }
  if (norm.method == 'TMM' | norm.method == 'UQ') {
    BiocGenerics::sizeFactors(dds) <- rep(1, length(colnames(edgeR::cpm(dge))))
    if (run.dispersion) {
      dds <- DESeq2::estimateDispersions(dds)
    } else {
      DESeq2::dispersions(dds) <- rep(1, length(rownames(edgeR::cpm(dge))))
    }
    if (debug) print("size factor is")
    if (debug) print(BiocGenerics::sizeFactors(dds))
    if (debug) print("dispersion head is")
    if (debug) print(head(DESeq2::dispersions(dds)))
    dds <- DESeq2::nbinomWaldTest(dds)
  } else {
    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- DESeq2::estimateDispersions(dds)
    dds <- DESeq2::nbinomWaldTest(dds)
    # dds           <- DESeq2::DESeq(dds)
  }
  ## ---
  if (debug) print("resultsNames is")
  if (debug) print(DESeq2::resultsNames(dds))
  ## ---
  res           <- DESeq2::results(dds, contrast=c("group", 'group1', 'group2' ), format="DataFrame")
  res           <- as.data.frame(res)
  tpdeseq2.save <- res[order(res$padj) , ]
  if (debug) print(head(tpdeseq2.save))
  # print('090909090909')
  # print(head(countTab))
  # print('-=-=-=-')
  # print(head(normCount))
  # print('-=-=-=-')
  # print('090909090909')
  return(list(deres = tpdeseq2.save, orgCount = countTab, normCount =  normCount, metaTab = metaTab))
  print("DESeq2 bulk analysis is completed.")
}



