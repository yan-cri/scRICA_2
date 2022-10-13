rm(list = ls())
## -------------------------------------------------------------------------------------##
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scRICA))
parser              <- ArgumentParser()
## -------------------------------------------------------------------------------------##
## options
parser$add_argument("--outputDir", type="character",
                    default = paste(sprintf('%s/lightHippo_results', getwd())),
                    help="Full path to the lighthippo analysis results directory",
                    metavar='')

parser$add_argument("--rdsFname", type="character",
                    help="rds file name",
                    metavar='')

parser$add_argument("--cellAnnotation", type="character",
                    help="cluster annotation Rscript file name",
                    metavar='')

parser$add_argument("--cellclusters", type="character",
                    default = 'all',
                    help="Selected cluster levels based on the idents of input rds file",
                    metavar='')

parser$add_argument("--expCondCheck", type="character",
                    help = "options shown in RDS column names: e.g. 'sample', 'expCond1', or 'expCond2'...",
                    metavar='')

parser$add_argument("--expCondCheckFname", type="character",
                    help="The DE analysis results file name prefix",
                    metavar='')

parser$add_argument("--topNo", type="integer",
                    default = '10',
                    help="Output thses top number of differentiall expressed genes from each cell cluster.",
                    metavar='')

parser$add_argument("--deMethod", type="character",
                    default = 'wilcox',
                    help="The DE analysis method, options are 'wilcox', 't', 'negbinom', 'poisson', 'MAST', 'DESeq2', 'DESeq2.bulk'.",
                    metavar='')

parser$add_argument("--deseq2bulkSampleMeta", type="character",
                    default = NULL,
                    help="if deMethod is 'DESeq2.bulk', if specified DE analysis will be regressed on the provided variables in the meta table.",
                    metavar='')

parser$add_argument("--pAdjValCutoff", type="double",
                    default = '0.05',
                    help="The threshold of adjust p-values for DEGs identification.",
                    metavar='')

parser$add_argument("--min.cells.group", type="integer",
                    default = '10',
                    help="The minimum number of cells presenting in each comparison groups.",
                    metavar='')

parser$add_argument("--min.pct", type="double",
                    default = '0.1',
                    help="The minimum percentage of cells presenting in each comparison groups.",
                    metavar='')

parser$add_argument("--logfc.threshold", type="double",
                    default = '0.25',
                    help="The FC (log scale) threshold between defined comparison groups.",
                    metavar='')

parser$add_argument("--debug", type="character",
                    default = 'F',
                    help="Whether to print extra analysis processing information for debug, by default False",
                    metavar='')
## ---
args                <- parser$parse_args()
print('-=-=-=-=-=-=-=-=-=-')
print("All input options")
print(args)
print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
## -------------------------------------------------------------------------------------##
debug   = as.logical(args$debug)
print(sprintf("START expCond cluster markers identification anlaysis "))
if (!is.null(args$deseq2bulkSampleMeta)){
  sampMeta <- read.delim(file = as.character(args$deseq2bulkSampleMeta), header = T, sep = '\t')
} else {
  sampMeta <- NULL
}
## ------
print("START read in RDS file")
rds1             <- readRDS(as.character(args$rdsFname))
setwd(as.character(args$outputDir))
print("Complete read in RDS file")
print('-=-=-=-=-=-=-=-=-=-')
## ---
if (is.null(args$cellAnnotation)) {
  if (args$cellclusters == 'all') {
    if (args$deMethod=='DESeq2.bulk') {
      deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = F,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        deMethod = as.character(args$deMethod),
                                        deseq2bulk.metaCovariateInput = sampMeta,
                                        topNo = as.numeric(args$topNo),
                                        min.cells.group = as.numeric(args$min.cells.group),
                                        min.pct =  as.numeric(args$min.pct),
                                        logfc.threshold = as.numeric(args$logfc.threshold),
                                        pAdjValCutoff = as.numeric(args$pAdjValCutoff), debug = debug)
    } else {
      deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = F,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        deMethod = as.character(args$deMethod),
                                        topNo = as.numeric(args$topNo),
                                        min.cells.group = as.numeric(args$min.cells.group),
                                        min.pct =  as.numeric(args$min.pct),
                                        logfc.threshold = as.numeric(args$logfc.threshold),
                                        pAdjValCutoff = as.numeric(args$pAdjValCutoff), debug = debug)
    }

  } else {
    sel.cellClusters <- paste(unlist(strsplit(args$cellclusters, split = ', ')), sep = ', ' )
    if (args$deMethod=='DESeq2.bulk') {
      deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = F,
                                        cellcluster = sel.cellClusters,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        deMethod = as.character(args$deMethod),
                                        deseq2bulk.metaCovariateInput = sampMeta,
                                        topNo = as.numeric(args$topNo),
                                        min.cells.group = as.numeric(args$min.cells.group),
                                        min.pct =  as.numeric(args$min.pct),
                                        logfc.threshold = as.numeric(args$logfc.threshold),
                                        pAdjValCutoff = as.numeric(args$pAdjValCutoff), debug = debug)
    } else {
      deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = F,
                                        cellcluster = sel.cellClusters,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        deMethod = as.character(args$deMethod),
                                        topNo = as.numeric(args$topNo),
                                        min.cells.group = as.numeric(args$min.cells.group),
                                        min.pct =  as.numeric(args$min.pct),
                                        logfc.threshold = as.numeric(args$logfc.threshold),
                                        pAdjValCutoff = as.numeric(args$pAdjValCutoff), debug = debug)
    }

  }
  print(sprintf("END expCond cluster markers identification anlaysis "))
  save(deRes, file = file.path(args$outputDir, sprintf("clusterMarkerGenes_results_wOrgClusterAnnotation/%s.Rdata", args$expCondCheckFname)))
} else {
  if (args$cellclusters == 'all') {
    if (args$deMethod=='DESeq2.bulk') {
      deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        deMethod = as.character(args$deMethod),
                                        deseq2bulk.metaCovariateInput = sampMeta,
                                        topNo = as.numeric(args$topNo),
                                        min.cells.group = as.numeric(args$min.cells.group),
                                        min.pct =  as.numeric(args$min.pct),
                                        logfc.threshold = as.numeric(args$logfc.threshold),
                                        pAdjValCutoff = as.numeric(args$pAdjValCutoff), debug = debug)
    } else {
      deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        deMethod = as.character(args$deMethod),
                                        topNo = as.numeric(args$topNo),
                                        min.cells.group = as.numeric(args$min.cells.group),
                                        min.pct =  as.numeric(args$min.pct),
                                        logfc.threshold = as.numeric(args$logfc.threshold),
                                        pAdjValCutoff = as.numeric(args$pAdjValCutoff), debug = debug)
    }

  } else {
    sel.cellClusters <- paste(unlist(strsplit(args$cellclusters, split = ', ')), sep = ', ' )
    if (args$deMethod=='DESeq2.bulk') {
      deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                        cellcluster = sel.cellClusters,
                                        newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        deMethod = as.character(args$deMethod),
                                        deseq2bulk.metaCovariateInput = sampMeta,
                                        topNo = as.numeric(args$topNo),
                                        min.cells.group = as.numeric(args$min.cells.group),
                                        min.pct =  as.numeric(args$min.pct),
                                        logfc.threshold = as.numeric(args$logfc.threshold),
                                        pAdjValCutoff = as.numeric(args$pAdjValCutoff), debug = debug)
    } else {
      deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                        cellcluster = sel.cellClusters,
                                        newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        deMethod = as.character(args$deMethod),
                                        topNo = as.numeric(args$topNo),
                                        min.cells.group = as.numeric(args$min.cells.group),
                                        min.pct =  as.numeric(args$min.pct),
                                        logfc.threshold = as.numeric(args$logfc.threshold),
                                        pAdjValCutoff = as.numeric(args$pAdjValCutoff), debug = debug)
    }

  }
  print(sprintf("END expCond cluster markers identification anlaysis "))
  save(deRes, file = file.path(args$outputDir, sprintf("clusterMarkerGenes_results_wNewAnnotation/%s.Rdata", args$expCondCheckFname)))

}
## -------------------------------------------------------------------------------------##
print('END==========END==========END==========END==========END')
## -------------------------------------------------------------------------------------##
