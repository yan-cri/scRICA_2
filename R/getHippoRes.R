#' getHippoRes() Function
#' @details
#' This function is used to run hippo analysis on specified cell clusters and save hippo results into Rdata.
#'
#' @param resDir specify an exiting full path of directory, where hippo clustering analysis results will be saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved with above rdsDir option.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param expCondCheck specify which experimental conditions to be explored, including sample or expCond1/2/....
#' @param cellcluster specify cell clusters (idents) to be conducted with hippo analysis.
#' @param expCond specify the specific experimental conditions to be conducted with hippo analysis.
#' @param hippoResNamePrefix prefix of the hippo analysis results, if not defined, by default = 'hippo_cluster_test'.
#' @param noClusters number of clusters for hippo to identify, default = 3.
#' @param sparseMatrix whether to turn sparse Matrix option on, default off, when turned on, takes longer time to run hippo.
#' @param initial.label.on whether to run hippo starting with initial identified cell clusters information.
#' @param topN = 100 by default, the top number genes plotted in hippo analysis
#'
#' @importFrom Seurat Idents
#' @importFrom SeuratObject DefaultAssay
#' @importFrom utils write.table
#' @importFrom tools file_ext
#' @importFrom utils read.delim
#' @importFrom lightHippo cut_hierarchy
#' @importFrom lightHippo organizing_hippo_features
#' @importFrom lightHippo summarize_current_zero_proportions
#' @importFrom lightHippo cut_hierarchy
#' @importFrom lightHippo summarize_for_feature_dot
#' @importFrom lightHippo makeDotplot
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom ggplot2 ggsave
#' @keywords getHippoRes
#' @examples getHippoRes(rds,  expCondCheck='sample/expCond*', cellcluster)
#' @export
#'
#' @return
#' a results directory/folder with saved Rdata hippo analysis results inside
#'
#' ## ------------------------------------------------------------------------------------ ##

getHippoRes <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL, expCondCheck='all',
                        cellcluster = NULL , expCond = NULL, noClusters = 3, sparseMatrix = F, initial.label.on = F,
                        hippoResNamePrefix = 'hippo_cluster_test', topN = 100, debug = as.logical(F)) {
  ##--------------------------------------------------------------------------------------##
  if (is.null(cellcluster)) stop("Please provide 'cellcluster' for hippo analysis")
  sel.expConds <- expCond ## to avoid 'expCond' otpion with metadata 'expCond' rename this paratmer into sel.expConds
  ## ----
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) stop("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
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
  resDir                <- paste(resDir, 'hippo_results', sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  resDir                <- paste(resDir, hippoResNamePrefix, sep = '/')
  if (!dir.exists(resDir)) dir.create(resDir)
  ##--------------------------------------------------------------------------------------##
  # if (expCondCheck == 'sample') {
  #   if (is.null(expCondCheckFname)) {
  #     expCondCheckFname        <- 'expCond_sample'
  #   } else {
  #     expCondCheckFname        <- expCondCheckFname
  #   }
  # } else {
  #   if (is.null(expCondCheckFname)) {
  #     expCondCheckFname        <- as.character(expCondCheck)
  #   } else {
  #     expCondCheckFname        <- expCondCheckFname
  #   }
  # }
  ##--------------------------------------------------------------------------------------##
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    source(as.character(newAnnotationRscriptName))
    print('-=-=-=-')
    print('updated Idents are as below:')
    print(table(Idents(seuratObjFinal)))
    print('-=-=-=-')
  }
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond'
  if (expCondCheck == 'all') {
    seuratObjFinal                     <- seuratObjFinal
  } else {
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  ##--------------------------------------------------------------------------------------##
  ## if provided, subset on 'cellcluster'
  orgClusterLevels <- levels(Seurat::Idents(seuratObjFinal))
  if (debug) print(sprintf("%s idents levels are %s", length(orgClusterLevels), paste(orgClusterLevels, collapse = ', ')))
  if (!is.null(cellcluster)) {
    if (any(!cellcluster %in% orgClusterLevels ) ) stop('Please provide the corresponding cell clusters in identfied idents for hippo analysis.')
    print(sprintf('Subsetting %s specific cell clusters: %s', length(cellcluster), paste(cellcluster, collapse = ',')))
    seuratObjFinal        <- subset(seuratObjFinal, idents = as.character(cellcluster) )
  }
  ## in provided, subset on 'expCond'
  expCondLevels           <- levels(factor(seuratObjFinal@meta.data$expCond))
  if (debug) print(sprintf("%s expCond levels are %s", length(expCondLevels), paste(expCondLevels, collapse = ', ')))
  if (!is.null(sel.expConds)) {
    if (any(!sel.expConds %in% expCondLevels ) ) stop("Please provide the corresponding experimental condition levels specified in 'expCondCheck' option.")
    print(sprintf('Subsetting %s specific expCond: %s', length(sel.expConds), paste(sel.expConds, collapse = ',')))

    if (length(sel.expConds)==1) {
      seuratObjFinal          <- subset(seuratObjFinal, expCond == sel.expConds)
    } else {
      for (i in 1:length(sel.expConds)) {
        seuratObjFinalPrep    <- subset(seuratObjFinal, subset = expCond == as.character(sel.expConds[i]))
        if (i ==1) {
          seuratObjFinalPrep2 = seuratObjFinalPrep
        } else {
          seuratObjFinalPrep2 <- merge(seuratObjFinalPrep2, seuratObjFinalPrep)
        }
      }
      seuratObjFinal          <- seuratObjFinalPrep2
    }
  }
  # print(sprintf("hippo input is as below"))
  # seuratObjFinal
  ##--------------------------------------------------------------------------------------##
  inputDataPrep <- seuratObjFinal
  if (sparseMatrix) {
    print("lighthippo input is a sparse matrix")
    inputData               <- inputDataPrep@assays$RNA@counts ## sparse matrix
  } else {
    print("lighthippo input is a dense matrix")
    inputData               <- as.matrix(inputDataPrep@assays$RNA@counts) ##dense matrix
  }
  # print(head(inputData[,1:4]))
  # print(dim(inputData))
  # print('7574655647578282-=-=-=-=-=-=-=-=-=-')
  print(sprintf('lightHippo: %s cells in combined cell cluster (%s) with a total of %s genes expressed', dim(inputData)[2], paste(cellcluster, collapse = '; '), dim(inputData)[1] ))
  print('-=-=-=-=-=')
  ##--------------------------------------------------------------------------------------##
  ## 1. run light hippo, returned
  systime1         <- Sys.time()
  if (initial.label.on) {
    inputData.cluster  <- as.numeric(Idents(inputDataPrep))
    print(sprintf("Processing analysis with 'initial.labels' on for 'K.round' = %s", noClusters))
    set.seed(1234)
    K.round            <- noClusters - length(unique(inputData.cluster)) + 1
    lightHippoRes      <- lightHippo::lightHIPPO(dat = inputData, K.round = K.round, initial.labels = inputData.cluster)
    final_feature_list <- lightHippo::organizing_hippo_features(lightHippoRes)
    save(inputData, lightHippoRes, inputData.cluster, noClusters, cellcluster, hippoResNamePrefix, file = file.path(resDir, sprintf('lightHIPPOres_%s_k%s.Rdata', hippoResNamePrefix, noClusters) ) )
  } else {
    print(sprintf("Processing analysis with 'initial.labels' off for 'K.round' = %s", noClusters))
    set.seed(1234)
    lightHippoRes      <- lightHippo::lightHIPPO(dat = inputData, K.round = noClusters, initial.round = 0)
    final_feature_list <- lightHippo::organizing_hippo_features(lightHippoRes)
    save(inputData, lightHippoRes, noClusters, cellcluster, hippoResNamePrefix, file = file.path(resDir, sprintf('lightHIPPOres_%s_k%s.Rdata', hippoResNamePrefix, noClusters) ) )
  }
  systime2         <- Sys.time()
  print(sprintf("lighthippo complete used %s %s.", round(difftime(systime2, systime1), digits = 2), attr(difftime(systime2, systime1), "units") ))
  print("lighthippo next round ID is shown as below:")
  print(table(lightHippoRes$next_round_IDs))
  print('-=-=-=-')
  ##--------------------------------------------------------------------------------------##
  ## 2. make diagnostic plot
  total.num.gene  <- nrow(inputData)
  set.seed(20200610)
  randomIDs       <- sample(1:total.num.gene, 5000)
  summarizing_dat <- lightHippo::summarize_current_zero_proportions(inputData[randomIDs, ], lightHippoRes$next_round_IDs)
  plot_dat_per_cluster_inflation <- lightHippo::visualize_current_zero_proportions(summarizing_dat)
  if (noClusters <= 3) {
    pdf(file = file.path(resDir, sprintf('diagnosticPlot_cellCluster_%s_k%s.pdf', hippoResNamePrefix, noClusters) ), width = 2*(noClusters+1), height = ceiling((noClusters+1)/4)*3)
  } else {
    pdf(file = file.path(resDir, sprintf('diagnosticPlot_cellCluster_%s_k%s.pdf', hippoResNamePrefix, noClusters) ), width = 8, height = ceiling(noClusters/3.9)*3)
  }
  # pdf("lightHIPPO_counts_inflation_check.png")
  print(plot_dat_per_cluster_inflation)
  dev.off()
  ##--------------------------------------------------------------------------------------##
  ttID       <- lightHippo::cut_hierarchy(lightHippoRes, K = noClusters)
  print(sprintf("%s clusters table is:", noClusters))
  print(table(ttID))
  print('-=-=-=-')
  check_these <- final_feature_list$ID[[noClusters-1]]
  tt.selected <- lightHippo::summarize_for_feature_dot(inputData[check_these , ], ttID)
  p <-  lightHippo::makeDotplot(input = tt.selected, topN = topN)
  if (noClusters <= 4) {
    ggsave(filename = file.path(resDir, sprintf('topFeatures_cellCluster_%s_k%s_dotplot.pdf', hippoResNamePrefix, noClusters) ), plot = p, width = ceiling(topN/4), height = 5)
  } else if (noClusters>4 & noClusters <=6) {
    ggsave(filename = file.path(resDir, sprintf('topFeatures_cellCluster_%s_k%s_dotplot.pdf', hippoResNamePrefix, noClusters) ), plot = p, width = ceiling(topN/4), height = noClusters)
  } else {
    ggsave(filename = file.path(resDir, sprintf('topFeatures_cellCluster_%s_k%s_dotplot.pdf', hippoResNamePrefix, noClusters) ), plot = p, width = ceiling(topN/4), height = ceiling(noClusters/1.5))
  }
  ##--------------------------------------------------------------------------------------##
  viewRes  <- lightHippo::cut_hierarchy(lightHippoRes, K = noClusters, cut_sequence = T)
  pdf(file = file.path(resDir, sprintf('%s_k%s_clustersTree.pdf', hippoResNamePrefix, noClusters) ), width = 5, height = 4 )
  lightHippo::visualize_hippo_hierarchy(viewRes)
  dev.off()
  ##--------------------------------------------------------------------------------------##
  if (initial.label.on) {
    noStart = length(unique(inputData.cluster))
  } else {
    noStart = 2
  }
  for (j in noStart:noClusters) {
    res          <- lightHippo::cut_hierarchy(lightHippoRes, K = j)
    cellHippoPer <- as.data.frame( table(res))
    colnames(cellHippoPer) <- c('hippo', 'No cells')
    fname <- file.path(resDir, sprintf('Per_cellCluster_%s_k%s.xlsx', hippoResNamePrefix, noClusters) )
    # print(cellHippoPer)
    if (j ==1) {
      xlsx::write.xlsx(x = cellHippoPer, file = fname, sheetName = paste('k', j, sep = '_'), row.names = F, append = F)
    } else {
      xlsx::write.xlsx(x = cellHippoPer, file = fname, sheetName = paste('k', j, sep = '_'), row.names = F, append = T)
    }
  }
  print(sprintf("light Hippo analysis complete for k=%s", noClusters))
  print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
  ##--------------------------------------------------------------------------------------##
}
