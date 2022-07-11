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

## ---
args                <- parser$parse_args()
print('-=-=-=-=-=-=-=-=-=-')
print("All input options")
print(args)
print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
## -------------------------------------------------------------------------------------##
print("START read in RDS file")
rds1             <- readRDS(as.character(args$rdsFname))
setwd(as.character(args$outputDir))
print("Complete read in RDS file")
print('-=-=-=-=-=-=-=-=-=-')
## ---
print(sprintf("START expCond cluster markers identification anlaysis "))
if (is.null(args$cellAnnotation)) {
  deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                    newAnnotation = F,
                                    expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                    topNo = as.numeric(args$topNo))
  print(sprintf("END expCond cluster markers identification anlaysis "))
  save(deRes, file = file.path(args$outputDir, sprintf("clusterMarkerGenes_results_wOrgClusterAnnotation/%s.Rdata", args$expCondCheckFname)))
} else {
  deRes <- getExpCondClusterMarkers(rds = rds1, resDir = as.character(args$outputDir),
                                    newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                    expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                    topNo = as.numeric(args$topNo))
  print(sprintf("END expCond cluster markers identification anlaysis "))
  save(deRes, file = file.path(args$outputDir, sprintf("clusterMarkerGenes_results_wNewAnnotation/%s.Rdata", args$expCondCheckFname)))

}
## -------------------------------------------------------------------------------------##
print('END==========END==========END==========END==========END')
## -------------------------------------------------------------------------------------##
