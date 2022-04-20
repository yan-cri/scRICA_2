## R script used to conduct scRNA-Seq comparative integration analysis based scRICA fns ##
##--------------------------------------------------------------------------------------##
# remove.packages('scRICA')
# .rs.restartR()
rm(list = ls())
## ---
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scRICA))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(dplyr))
currentDir          <- getwd()
## ---------------------------------------------------------------------------------------
parser              <- ArgumentParser()
## add parser argument options
parser$add_argument("-m", "--metadata", nargs = 1,
                    type="character", help="Required: Full path to metadata table",
                    metavar = '')
parser$add_argument("-o", "--outputDirName", type="character",
                    default = 'output_results',
                    help="Full path to save analysis results directory name, default 'output_results' ",
                    metavar = '')
## ------
## options used for countReadin()
## common option used for different functions
parser$add_argument("-g", "--genomeSpecies", type="character",
                    default = 'human',
                    help="genome species, options are human, mice, rat",
                    metavar = '')
parser$add_argument("-mc", "--minCells", type="integer",
                    default = '3',
                    help="customized option used to include features detected in at least this many cells, by default 3",
                    metavar = '')
parser$add_argument("-mf", "--minFeatures", type="integer",
                    default = '200',
                    help="customized option used to include cells where at least this many features are detected, by default 200",
                    metavar = '')
parser$add_argument("-mtF", "--mtFiltering", type="character",
                    default = 'T',
                    help="logical, whether to filter out mitochondrial content for downstream analysis",
                    metavar = '')
parser$add_argument("-mtP", "--mtPerCutoff", type="integer",
                    default = '20',
                    help="if mtFilter option is on, set customized mitochondrial content filtering percentage cutoff level, by default 20",
                    metavar = '')
parser$add_argument("--multiomics", type="character",
                    default = 'F',
                    help="logical, whether is this multiomics data",
                    metavar = '')
parser$add_argument("--extraFilter", type="character",
                    default = 'T',
                    help="logical, whether extra cells filtering implemented with information in the metadata table column 'filterFname'.",
                    metavar = '')
parser$add_argument("--integration.method", type="character",
                    default = 'CCA',
                    help="2 options, CCA or RPCA as integration method.",
                    metavar = '')
parser$add_argument("--integration.nfeatures", type="integer",
                    default = '2000',
                    help="Number of features used for integration, by default 2000",
                    metavar = '')
## ------
args <- parser$parse_args()
print(args)
print(sprintf(" mtFiltering = %s ", as.logical(args$mtFiltering)))
print('---===---')
## ---------------------------------------------------------------------------------------
## 0. prepare and output dir to save analysis results
outputDir               <- paste(getwd(), as.character(args$outputDirName), sep = '/')
print(sprintf('Analysis results are saved in %s', outputDir))
## ---
metadata                <- read.delim2(file = as.character(args$metadata), header = T)
print(head(metadata))
print(sprintf("%s samples in provided metadata table to be processed", dim(metadata)[1]))
## ---------------------------------------------------------------------------------------
## 1. parameters/options used for step1: doublets identification analysis with findDoublets()
seuratProcssedObjList   <- processQC(metadata = metadata,
                                     multiomics = as.logical(args$multiomics),
                                     extraFilter = as.logical(args$extraFilter),
                                     resDirName = args$outputDirName,
                                     genomeSpecies = args$genomeSpecies,
                                     minCells = args$minCells, minFeatures = args$minFeatures,
                                     mtFiltering = as.logical(args$mtFiltering), mtPerCutoff = args$mtPerCutoff)
save(seuratProcssedObjList, file = file.path(seuratProcssedObjList$resDir, "seuratReadin_qcProcessed.Rdata"))
print("END-=-=-=-=-=-END-=-=-=-=-=-END-=-=-=-=-=-END")
## ---------
## 2. integration
seuratIntegratedRes     <- getClusterMarkers(qcProcessedResults = seuratProcssedObjList,
                                             integrationMethod = as.character(args$integration.method),
                                             nfeatures = as.numeric(args$integration.nfeatures),
                                             topN = 50)
str(seuratIntegratedRes)
## ---------
setwd(currentDir)
## END---------------------------------------END---------------------------------------END
