## R functions used to identify doublets/multiplets with DoubletDecon                   ##
## Developed by Yan Li, Jan, 2021                                                       ##
##--------------------------------------------------------------------------------------##
#' meata2list() Function
#'
#' This function allows you to find and estimate doublets with DoubletDecon for the provide iput 'cellrangerResList'
#' @param metadata list including full path to cellranger analysis results for different samples.
#' @keywords meata2list
#' @export
#' @examples meata2list()
#' @return
#'
##----------------------------------------------------------------------------------------
meata2list <- function(metadata) {
  cellrangerResList    <- list()
  for (i in 1:length(metadata$sample)) {
    cellrangerResList[[i]] <- as.character(metadata$path[i])
  }
  names(cellrangerResList) <- metadata$sample
  return(cellrangerResList)
}
##----------------------------------------------------------------------------------------
