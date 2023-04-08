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

parser$add_argument("--hippoResNamePrefix", type="character",
                    help="name prefix of analysis output folder",
                    metavar='')

parser$add_argument("--cellAnnotation", type="character",
                    help="cluster annotation Rscript file name",
                    metavar='')

parser$add_argument("--expCondCheck", type="character",
                    default = 'all',
                    help = "all or 3 options: 'sample', 'expCond1', or 'expCond2'",
                    metavar='')

parser$add_argument("--expCond", type="character",
                    default = 'all',
                    help = "when expCondCheck is not all, define selected experimental conditions to implment hippo analysis when expCondCheck is no 'all'.",
                    metavar='')

parser$add_argument("--rdsFname", type="character",
                    help="rds file name",
                    metavar='')

parser$add_argument("--cellclusters", type="character",
                    help="Selected cluster levels based on the idents of input rds file",
                    metavar='')


parser$add_argument("--initialLable", type="character",
                    default = 'False',
                    help="Logical, where to run lighthippo starting with inital cluster",
                    metavar='')

parser$add_argument("--sparseMatrix", type="character",
                    default = 'False',
                    help="Logical, where to run lighthippo with sparse matrix, if False, run with dense matrix",
                    metavar='')

parser$add_argument("--noCluster", type="character",
                    default = '3',
                    help="List of numeric values for the number of cluster lists, default=3",
                    metavar='')

## ---
args                <- parser$parse_args()
print(args)
## -------------------------------------------------------------------------------------##
rds1             <- readRDS(as.character(args$rdsFname))
setwd(as.character(args$outputDir))
sel.cellClusters <- paste(unlist(strsplit(args$cellclusters, split = ', ')), sep = ', ' )
sel.expConds     <- paste(unlist(strsplit(args$expCond, split = ', ')), sep = ', ' )

## ---
if (is.null(args$cellAnnotation)) {
  if (args$expCondCheck == 'all') {
    getHippoRes(rds = rds1, newAnnotation = F,
                cellcluster  = sel.cellClusters, noClusters = as.numeric(args$noCluster),
                sparseMatrix = as.logical(args$sparseMatrix), initial.label.on = as.logical(args$initialLable),
                hippoResNamePrefix = args$hippoResNamePrefix)
  } else {
    getHippoRes(rds = rds1, newAnnotation = F,
                expCondCheck = as.character(args$expCondCheck),
                cellcluster  = sel.cellClusters, expCond = sel.expConds, noClusters = as.numeric(args$noCluster),
                sparseMatrix = as.logical(args$sparseMatrix), initial.label.on = as.logical(args$initialLable),
                hippoResNamePrefix = args$hippoResNamePrefix)
  }

} else {
  if (args$expCondCheck == 'all') {
    getHippoRes(rds = rds1, newAnnotation=T, newAnnotationRscriptName=args$cellAnnotation,
                cellcluster  = sel.cellClusters, noClusters = as.numeric(args$noCluster),
                sparseMatrix = as.logical(args$sparseMatrix), initial.label.on = as.logical(args$initialLable),
                hippoResNamePrefix = args$hippoResNamePrefix)
  } else {
    getHippoRes(rds = rds1, newAnnotation=T, newAnnotationRscriptName=args$cellAnnotation,
                expCondCheck = as.character(args$expCondCheck),
                cellcluster  = sel.cellClusters, expCond = sel.expConds, noClusters = as.numeric(args$noCluster),
                sparseMatrix = as.logical(args$sparseMatrix), initial.label.on = as.logical(args$initialLable),
                hippoResNamePrefix = args$hippoResNamePrefix)
  }
}
## -------------------------------------------------------------------------------------##
print('END==========END==========END==========END==========END')
## -------------------------------------------------------------------------------------##



