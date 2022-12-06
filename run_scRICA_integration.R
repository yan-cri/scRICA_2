## R script used to conduct scRNA-Seq comparative integration analysis based on scRICA fn getClusterMarkers()  ##
##--------------------------------------------------------------------------------------##
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
parser$add_argument("--seuratObjListInput", nargs = 1,
                    type="character", help="Required: seurat object list, which contains a list of seurat object to be integrated",
                    metavar = '')
parser$add_argument("--resDirName", type="character",
                    default = 'output_results',
                    help="Analysis results directory name, default 'output_results' ",
                    metavar = '')
parser$add_argument("--integration.nfeatures", type="integer",
                    default = '2000',
                    help="number of features used in the integration step.",
                    metavar = '')
parser$add_argument("--integration.method", type="character",
                    default = 'CCA',
                    help="2 options, CCA or RPCA as integration method.",
                    metavar = '')
parser$add_argument("--integration.k.weight", type="integer",
                    default = '100',
                    help="Number of neighbors to consider when weighting anchors during the integration step.",
                    metavar = '')
## ------
args <- parser$parse_args()
print(args)
print('---===---')
## ---------------------------------------------------------------------------------------
## 0. prepare and output dir to save analysis results
outputDir               <- paste(getwd(), as.character(args$resDirName), sep = '/')
print(sprintf('Analysis results are saved in %s', outputDir))
## ---------------------------------------------------------------------------------------
print("load input data")
load(as.character(args$seuratObjListInput))
print(seuratObjList)
print("STEP1-=-=-=-=-=-STEP1-=-=-=-=-=-STEP1-=-=-=-=-=-STEP1")
## ---------
# ?getClusterMarkers() changed into getIntegrationClusterMarkers()
seuratIntegratedRes     <- getClusterMarkers(qcProcessedResults = seuratObjList, int.k.weight = as.numeric(args$integration.k.weight),
                                             resDirName = as.character(args$resDirName),
                                             integrationMethod = as.character(args$integration.method),
                                             nfeatures = as.numeric(args$integration.nFeatures), topN = 50)
str(seuratIntegratedRes)
## ---------
setwd(currentDir)
print("END-=-=-=-=-=-END-=-=-=-=-=-END-=-=-=-=-=-END")
## END---------------------------------------END---------------------------------------END
