#' getCellChatRes() Function
#' @details
#' This function is used to conduct cellChat analysis on the selected expCond levels.
#'
#' @param resDir full path of integration results analysis are saved, where RDS file is saved inside the 'RDS_Dir'. This path is also returned by getClusterMarkers() execution.
#' @param rds User also can provide the full path of RDS file instead of 'resDir' where RDS file is saved in. If this option is used, please also provide 'resDir' to specify where the analysis results will be saved.
#' @param newAnnotation logical value to indicate whether to add the annotation for identified cell clusters from getClusterMarkers() integration analysis.
#' @param newAnnotationRscriptName if 'newAnnotation = T', please specify here for the full path of the R script where cell clusters are defined.
#' @param genomeSpecies specify sample's genome species, by default 'human', currently supporting human, mouse, and rate.
#' @param selected.DB whether to conduct cellChat analysis on the selected pathway, a total of 3 options to choose, they are 'Cell-Cell Contact', 'ECM-Receptor' and 'Secreted Signaling'.
#' @param expCondCheck specify which experimental conditions to be explored, including sample, idents, or expCond1/2/....
#' @param expCondCheckFname prefix of saved rds file name.
#' @param expCond specify the specific experimental condition for cellchat analysis, if not specified, all experimental conditions will be performed .
#' @param min.cells.filtering the minimum number of cells in the group used for analysis.
#'
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 element_line
#' @importFrom Seurat Idents
#' @importFrom Seurat NoLegend
#' @importFrom Seurat RotatedAxis
#' @importFrom SeuratObject DefaultAssay
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom utils write.table
#' @importFrom tools file_ext
#' @importFrom utils read.delim
#' @importFrom xlsx read.xlsx
#' @importFrom cowplot plot_grid
#' @importFrom cowplot get_legend
#' @importFrom grid unit
#'
#' @keywords getCellChatRes
#' @examples getCellChatRes(rds, expCondCheck='sample/idents/expCond*')
#' @export
#' @return cellChat analysis result
## ---------------------------------------------------------------------------------------
# library(ggplot2)
# library(Seurat)
# library(grDevices)
# library(tools)
# library(xlsx)
getCellChatRes <- function(resDir=NULL, rds=NULL, newAnnotation=F, newAnnotationRscriptName=NULL,
                           genomeSpecies = 'human', selected.DB = NULL,
                           expCondCheck='sample', expCondCheckFname = NULL, expCond = NULL,
                           min.cells.filtering = 10, debug = F ){
  ## ---
  newAnnotation           <- as.logical(newAnnotation)
  if (newAnnotation & is.null(newAnnotationRscriptName)) print("Option 'newAnnotation' is on, please provide corresponding option 'newAnnotationRscriptName'.")
  sel.expConds <- expCond ## to avoid 'expCond' otpion with metadata 'expCond' rename this paratmer into sel.expConds
  ## ---
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
  ## ------
  ## update results directory if new annotation is used
  if (newAnnotation) {
    cellchatResDir        <- paste(resDir, 'results_wNewAnnotation_cellChat', sep = '/')
  } else {
    cellchatResDir        <- paste(resDir, 'results_wOrgClusterAnnotation_cellChat', sep = '/')
  }
  if (!dir.exists(cellchatResDir)) dir.create(cellchatResDir)
  ##--------------------------------------------------------------------------------------##
  if (newAnnotation) {
    ## Assign cell type identity to clusters
    ## redefine the level of Idents on the y-axis can be adjusted here by inputting order for cell annotation
    print("Before adding new annotation")
    print(table(Idents(seuratObjFinal)))
    source(newAnnotationRscriptName)
    print("After adding new annotation")
    print(table(Idents(seuratObjFinal)))
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
  if (!is.null(sel.expConds)) {
    expCondLevels.org = levels(factor(seuratObjFinal@meta.data$expCond))
    if (any(!sel.expConds %in% expCondLevels.org ) ) stop("Please provide the corresponding experimental condition levesl specified in 'expCondCheck' option.")
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
  ##--------------------------------------------------------------------------------------##
  if (debug) print(table(seuratObjFinal@meta.data$expCond))
  expCondLevels = levels(factor(seuratObjFinal@meta.data$expCond))
  ##--------------------------------------------------------------------------------------##
  cellChatResList <- list()
  for (i in 1:length(expCondLevels)) {
    print(sprintf("%s. Start to process cellChat analysis on '%s'.", i, expCondLevels[i]))
    seuratObjFinal.subset  <- subset(x = seuratObjFinal, expCond == expCondLevels[i])
    data.cells             <- seuratObjFinal.subset@assays$RNA@data
    print(sprintf("A total %s cells with %s gene expression used for cellChat analysis", dim(data.cells)[2], dim(data.cells)[1]))
    if (debug) print(sprintf("[gene number, cell number] = [%s, %s]", dim(data.cells)[1], dim(data.cells)[2]))
    meta.input            <- data.frame('expCondCheck' = seuratObjFinal.subset@meta.data[, grep(as.character(expCondCheck), colnames(seuratObjFinal.subset@meta.data))])
    meta.input$labels     <- as.character(Seurat::Idents(seuratObjFinal.subset))
    if (debug) print(head(meta.input))
    if (debug) print(table(meta.input$labels))
    if (debug) unique(meta.input$labels)
    ## --------------------------------
    cellchat <- CellChat::createCellChat(object = data.cells, meta = meta.input, group.by = "labels")
    cellchat <- CellChat::addMeta(cellchat, meta = meta.input)
    cellchat <- CellChat::setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
    if (debug) print(levels(cellchat@idents)) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    ## --------------------------------
    print(sprintf("genome species is %s.", genomeSpecies))
    if (genomeSpecies == 'human') {
      CellChatDB <- CellChat::CellChatDB.human
    } else if (genomeSpecies == 'mouse') {
      CellChatDB <- CellChat::CellChatDB.mouse
    }
    # showDatabaseCategory(CellChatDB)
    if (debug) table(CellChatDB$interaction$annotation)
    ## --------------------------------
    if (!is.null(selected.DB)) {
      CellChatDB.use <- CellChat::subsetDB(CellChatDB, search = as.character(selected.DB))
    } else {
      CellChatDB.use <- CellChatDB
    }
    cellchat@DB <- CellChatDB.use
    if (debug) table(CellChatDB.use$interaction$pathway_name)
    pathway.names <- names(table(CellChatDB.use$interaction$pathway_name))
    ## --------------------------------
    cellchat <- CellChat::subsetData(cellchat)
    # future::plan("multiprocess", workers = 4)
    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
    cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
    ## --------------------------------
    cellchat <- CellChat::computeCommunProb(cellchat)
    cellchat <- CellChat::filterCommunication(cellchat, min.cells = min.cells.filtering)
    cellchat <- CellChat::computeCommunProbPathway(cellchat)
    cellchat <- CellChat::aggregateNet(cellchat)
    ## --------------------------------
    ## Compute the network centrality scores
    cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    print(sprintf("%s. Complete cellChat analysis on '%s'.", i, expCondLevels[i]))
    saveRDS(cellchat, file = file.path(cellchatResDir, sprintf("%s_level_%s_cellchatRes.rds", expCondCheckFname, expCondLevels[i])))
    cellChatResList[[i]] = cellchat
    print('-=-=-=--=-=-=-==-=-=-=-=-=-')
  }
  names(cellChatResList) <- expCondLevels
  return(cellChatResList)
  ## --------------------------------
  ##--------------------------------------------------------------------------------------##
}
## ---------------------------------------------------------------------------------------
