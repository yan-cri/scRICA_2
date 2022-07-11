##--------------------------------------------------------------------------------------##
#' Improved_Seurat_Pre_Process() Function
#'
#' This function allows you to find and estimate doublets with DoubletDecon for the provide input 'cellrangerResList', taken from DoubletDecon directly due to function calling problems.
#' @param seuratObject list including full path to cellranger analysis results for different samples.
#' @param num_genes num_genes
#' @param write_files to be added
#' @param data_type to be added
#'
#' @importFrom SeuratObject UpdateSeuratObject
#' @importFrom Seurat FindAllMarkers
#' @importFrom Seurat Idents
#' @importFrom utils packageVersion
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr top_n
#'
#' @keywords DoubletDecon
#' @export
#' @examples Improved_Seurat_Pre_Process()
#'
##----------------------------------------------------------------------------------------
Improved_Seurat_Pre_Process <- function (seuratObject, num_genes = 50, write_files = FALSE, data_type = "counts") {

  version = packageVersion("Seurat")
  # seuratObject = UpdateSeuratObject(object = seuratObject)
  if (data_type == "counts") {
    expression = as.data.frame(seuratObject@assays[["RNA"]]@counts)
  }
  else if (data_type == "data") {
    expression = as.data.frame(seuratObject@assays[["RNA"]]@data)
  }
  else if (data_type == "scaled.data") {
    expression = as.data.frame(seuratObject@assays[["RNA"]]@scale.data)
  }
  seuratObject.markers = FindAllMarkers(object = seuratObject,
                                        only.pos = TRUE, min.pct = 0.25)
  if (version >= package_version(x = "3.9.9")) {
    genes = seuratObject.markers %>% group_by(cluster) %>%
      top_n(n = num_genes, wt = avg_log2FC)
  }
  else if (version >= package_version(x = "3.0.0") && version <
           package_version(x = "3.9.9")) {
    genes = seuratObject.markers %>% group_by(cluster) %>%
      top_n(n = num_genes, wt = avg_logFC)
  }
  else {
    print("This function only works with Seurat 3 or 4. Please update Seurat.")
  }
  clusters = as.data.frame(Idents(object = seuratObject))
  colnames(expression) = gsub("-", ".", colnames(expression))
  if (class(clusters[, 1]) == "character") {
    clusters[, 1] = gsub("-", ".", clusters[, 1])
  }
  else {
    row.names(clusters) = gsub("-", ".", row.names(clusters))
  }
  if (class(genes$cluster) == "factor") {
    if (min(as.numeric(as.character(genes$cluster))) == 0) {
      genes$cluster = as.numeric(as.character(genes$cluster)) +
        1
      clusters[, 1] = as.numeric(as.character(clusters[,
                                                       1])) + 1
    }
    else if (min(as.numeric(as.character(genes$cluster))) ==
             1) {
      genes$cluster = as.numeric(as.character(genes$cluster))
      clusters[, 1] = as.numeric(as.character(clusters[,
                                                       1]))
    }
    else {
      print("Unexpected cluster numbering scheme. Cluster numbers are expected to be continuous numbers starting from either 0 or 1. Please check conversion for correctness following this function.")
    }
  }
  clusters2 = clusters[order(clusters[, 1]), , drop = FALSE]
  expression = expression[row.names(clusters2)]
  genes2 = genes[order(genes$cluster), ]
  allgenes = expression
  expression = expression[row.names(expression) %in% as.character(genes$gene),
                          ]
  geneOrder = intersect(genes2$gene, as.character(row.names(expression)))
  expression = expression[match(geneOrder, row.names(expression)),
                          ]
  allgenes = rbind(clusters2[, 1], allgenes)
  expression = rbind(clusters2[, 1], expression)
  row.names(allgenes)[1] = "column_clusters-flat"
  row.names(expression)[1] = "column_clusters-flat"
  genes3 = genes2[match(geneOrder, genes2$gene), ]
  rowToAdd = c(NA, genes3$cluster)
  rowToAdd2 = rep(NA, nrow(allgenes))
  expression = cbind(rowToAdd, expression)
  allgenes = cbind(rowToAdd2, allgenes)
  colnames(expression)[1] = "row_clusters-flat"
  colnames(allgenes)[1] = "row_clusters-flat"
  groups = cbind(as.numeric(expression[1, 2:ncol(expression)]),
                 as.numeric(expression[1, 2:ncol(expression)]))
  row.names(groups) = as.character(colnames(expression)[2:ncol(expression)])
  if (write_files == TRUE) {
    write.table(expression, "ICGS_expression.txt", sep = "\t")
    write.table(allgenes, "ICGS_fullExpression.txt", sep = "\t")
    write.table(groups, "ICGS_groups.txt", sep = "\t", col.names = F)
  }
  return(list(newExpressionFile = expression, newFullExpressionFile = allgenes,
              newGroupsFile = groups))

}
