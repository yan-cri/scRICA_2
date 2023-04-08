## scRICA findDoublets(): identify/estimate doublets                                    ##
##--------------------------------------------------------------------------------------##
rm(list = ls())
## ---
suppressPackageStartupMessages(library(argparse))
# source('/Users/yanli/Desktop/842_multiomics_t2d/R/findDoublets2.R') ## this is used for batch1 (5+19samples), filter files colnames = 'ATAC_filter'
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scRICA))
currentDir          <- getwd()
## ---------------------------------------------------------------------------------------
parser              <- ArgumentParser()
## add parseroptions
parser$add_argument("-m", "--metadata", type="character", help="Required: Full path to metadata table",
                    metavar = '')
parser$add_argument("-o", "--resDirName", type="character",
                    default = 'output_results',
                    help="Full path to save analysis results directory name, default 'output_results' ",
                    metavar = '')
parser$add_argument("-g", "--genomeSpecies", type="character",
                    default = 'human',
                    help="genome species, options are human, mouse, human_mouse",
                    metavar = '')
parser$add_argument("--multiomics", type="character",
                    default = 'F',
                    help="logical, whether is this multiomics data",
                    metavar = '')
parser$add_argument("-rh", "--doubletDeconRhop", type="double",
                    default = '0.5',
                    help="doubletDecon Rhop, default 0.5",
                    metavar = '')
parser$add_argument("-p", "--doubletDeconPMF", type="logical",
                    default = 'F',
                    help="doubletDecon PMF, default is F, do not use step2 unique gene expression in doublet determination",
                    metavar = '')
parser$add_argument("-f", "--extraFilter", type="character",
                    default = 'F',
                    help="extraFilter, default is F, do not use extraFilter",
                    metavar = '')
parser$add_argument("-c", "--doubletDeconNoCore", type="double",
                    default = '-1',
                    help="customized doubletDeconNoCore, by default automatic detect",
                    metavar = '')
                    # , metavar='customized doubletDeconNoCore')
args <- parser$parse_args()
print(args)
print('---===---')
## ---------------------------------------------------------------------------------------
print(sprintf('metadataFname is %s', args$metadata))
print(sprintf('resDirName is %s', args$resDirName))
# ?findDoublets
## 1. parameters/options used for step1: doublets identification analysis with findDoublets()
metadata                <- read.delim2(file = as.character(args$metadata), header = T)
genomeSpecies           <- args$genomeSpecies
doubletDeconRhop        <- args$doubletDeconRhop
doubletDeconPMF         <- args$doubletDeconPMF
doubletDeconNoCore      <- args$doubletDeconNoCore
resFilename             <- args$resDirName
extraFilter             <- as.logical(args$extraFilter)
# print(str(extraFilter))
# print(sprintf("'extraFilter' is %s", extraFilter ))
# ## ---
# metadata                <- read.delim2(file = '/Users/yanli/Desktop/842_multiomics_t2d/metadata_batch1_5sampl_findDoubletsOnly.txt', header = T)
# metadata                <- read.delim2(file = '/Users/yanli/Desktop/842_multiomics_t2d/metadata_batch1_5sampl_findDoubletsOnly_wExtraFilter.txt', header = T)
# genomeSpecies           <- 'human'
# doubletDeconRhop        <- 0.5
# doubletDeconPMF         <- as.logical(T)
# doubletDeconNoCore      <- -1
# resFilename             <- 'findDoublets_results_scRICA_batch1_5samp'
## ---
doubletDeconResDir      <- findDoublets(metadata = metadata,
                                        multiomics = as.logical(args$multiomics),
                                        extraFilter = extraFilter,
                                        genomeSpecies = genomeSpecies,
                                        doubletDeconRhop =  doubletDeconRhop,
                                        doubletDeconPMF = doubletDeconPMF,
                                        doubletDeconNoCore = doubletDeconNoCore,
                                        resFilename = resFilename)
print(doubletDeconResDir)
setwd(currentDir)
## ---------------------------------------------------------------------------------------
