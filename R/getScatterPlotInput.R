# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(gridExtra)
## -------------------------------------------------------------------------------------- ##
#' getExpCondClusterMarkers() Function
#' @details
#' This function is used to identify positively expressed cluster markers for all originally identified/annotated cell clusters from the experimental conditions specified in the metadata table.
#'
#' @param resDir specify a exiting full path of integration results analysis are saved.
#' @param rds provide integrated RDS object, user can also provide the full path of the RDS where integrated RDS object is saved.
#' @param newAnnotation logical option, whether to add the new cell types annotation for identified cell clusters from provided integrated RDS file.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify the full path of the R script where new cell annotations are defined.
#' @param subsetOn whether to subset certain cell clusters when it is 'idents' or conditions when it is 'expCond1/2...'.
#' @param subset define the specific cell cluster names when it is 'subsetOn=idents' or expConds inheritated from metadata table 'subsetOn=expCond1/2...' to extract.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param selectedGroups specify 2 groups in a list to be plotted with scatter plot, it can be either cell ident names presented in 'Idents' of RDS object or experimental condition levels shown in the 'expCondCheck' condition.
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
#' @importFrom easylabel easylabel
#'
#' @keywords getScatterPlot
#' @examples getScatterPlot(rds, expCondCheck='Idents', selectedGroups = list('group1'=c('Ident1', 'Ident2'), 'group2' = c('Ident4', 'Ident5')))
#' @export
#' @return an interactive scatter plot window
## -------------------------------------------------------------------------------------- ##
getScatterPlot <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL,
                           subsetOn = NULL, subset = NULL, subsetGenes = NULL, expCondCheck='idents', selectedGroups, interactive = T) {
  if(length(selectedGroups)!=2) stop("Please provide corresponding x/y groups used for scatter plot")
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
    resDir <- paste(resDir, 'scatterPlot_wNewAnnotation', sep = '/')
  } else {
    resDir <- paste(resDir, 'scatterPlot_wOrgClusterAnnotation', sep = '/')
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
  if (!is.null(subsetOn)) {
    if (subsetOn=='idents') {
      ## ------
      ## subset on 'cellcluster'
      cellcluster = subset
      clusterLevels <- levels(Seurat::Idents(seuratObjFinal))
      if (!is.null(cellcluster)) {
        if (any(!cellcluster %in% clusterLevels ) ) stop('Please provide the corresponding cell clusters in identfied idents.')
        seuratObjFinal <- subset(seuratObjFinal, idents = as.character(cellcluster) )
        print(sprintf("Note: Only these cell clusters (%s) can be explored", paste(cellcluster, collapse = ', ')))
      } else {
        print("Note: All cell clusters definided in idents can be explored.")
      }
      ## ------
    } else {
      ## ------
      sel.expConds = subset
      if (!subsetOn%in%colnames(seuratObjFinal@meta.data)) {
        stop("ERROR: 'subsetOn' does not exist in your 'rds' metadata.")
      } else {
        seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(as.character(subsetOn), colnames(seuratObjFinal@meta.data))]
      }
      ## ------
      ## 'expCond' if provided, subset on 'expCond'
      expCondLevels           <- levels(factor(seuratObjFinal@meta.data$expCond))
      if (!is.null(sel.expConds)) {
        if (any(!sel.expConds %in% expCondLevels ) ) stop("Please provide the corresponding experimental condition levesl specified in 'subset' option.")
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
      ## ------
    }
  }
  ##--------------------------------------------------------------------------------------##
  ## update 'seuratObjFinal@meta.data$expCond'
  if (expCondCheck == 'idents') {
    seuratObjFinal@meta.data$expCond   <- Seurat::Idents(seuratObjFinal)
  } else {
    if (!expCondCheck%in%colnames(seuratObjFinal@meta.data)) {
      stop("ERROR: 'expCondCheck' does not exist in your 'rds' metadata.")
    } else {
      seuratObjFinal@meta.data$expCond <- seuratObjFinal@meta.data[, grep(sprintf('^%s$', as.character(expCondCheck)), colnames(seuratObjFinal@meta.data))]
    }
  }
  DefaultAssay(seuratObjFinal) <- "RNA"
  ##--------------------------------------------------------------------------------------##
  selected.groups <- as.character(unlist(selectedGroups))
  print("-=-=-=-=-=-=-=-=-")
  print("Start to extract expressions for scatter plot")
  res4plotPrep <- sapply(selected.groups, FUN = function(x) {
    subset.res <- subset(seuratObjFinal, expCond ==  x)
    norm.res <- subset.res@assays$RNA@data
    print(sprintf("%s cells for %s mean calculations", dim(norm.res)[2], x))
    norm.res.mean <- Matrix::rowMeans(x = norm.res)
    return(norm.res.mean)
  })
  ##--------------------------------------------------------------------------------------##
  ## subset genes if not null 'subsetGenes'
  if (!is.null(subsetGenes)) {
    res4plotPrep2 <- res4plotPrep[match(subsetGenes, rownames(res4plotPrep)),]
    res4plotPrep  <- res4plotPrep2[!is.na(rownames(res4plotPrep2)),]
    print(sprintf("subsetGenes is on, %s gene expression are extracted for scatter plot.", dim(res4plotPrep)[1]))
  } else {
    print(sprintf("subsetGenes is off, %s gene expression are extracted for scatter plot.", dim(res4plotPrep)[1]))
  }
  print("END to extract expressions for scatter plot")
  print("-=-=-=-=-=-=-=-=-")
  ## ---------
  res4plot <- sapply(selectedGroups, function(x) {
    if (length(x)==1) {
      res <- res4plotPrep[,match(x, colnames(res4plotPrep))]
    } else {
      res <- rowMeans(res4plotPrep[,match(x, colnames(res4plotPrep))])
    }
    return(res)
  })
  colnames(res4plot) <- gsub(pattern = '-|[ ]|/', replacement = '_', colnames(res4plot))
  if (is.null(colnames(res4plot))) {
    colnames(res4plot) <- c('Group1', 'Group2')
  }
  ##--------------------------------------------------------------------------------------##
  if (interactive) {
    easylabel::easylabel(res4plot, x = colnames(res4plot)[1] , y = colnames(res4plot)[2], output_shiny = T)
  } else {

  }
  ##--------------------------------------------------------------------------------------##
  return(res4plot)
}
