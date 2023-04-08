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

parser$add_argument("--cellclusters", type="character",
                    default = 'all',
                    help="Selected cluster levels based on the idents of input rds file",
                    metavar='')

parser$add_argument("--expCond", type="character",
                    default = 'all',
                    help = "Define selected experimental conditions for trajectory analysis, the conditions should be consistent with the levels in defined in the 'expCondCheck'.",
                    metavar='')

parser$add_argument("--slingshotclusterLabels", type="character",
                    default = 'GMM',
                    help = "Define slingshotclusterLabels option, by default 'GMM'.",
                    metavar='')

parser$add_argument("--topFeatureNo", type="integer",
                    default = '2000',
                    help="Number of top feature genes used in normalization and scaling precedure during pseudo-time trajectory analysis.",
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
print(sprintf("START pseudo-time trajectory analysis."))
if (is.null(args$cellAnnotation)) {
  if (args$expCond == 'all') {
    if (args$cellclusters == 'all') {
      print("Analysis is conducted on all experimental conditions specified in 'expCondCheck' for all cell types. ")
      pseudoRes <- getExpCondClusterPseudotime(rds = rds1, resDir = as.character(args$outputDir),
                                               newAnnotation = F,
                                               expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                               slingshotclusterLabels = as.character(args$slingshotclusterLabels),
                                               topFeatureNo = as.numeric(args$topFeatureNo) )
    } else {
      print("Analysis is conducted on all experimental conditions specified in 'expCondCheck' the selected cell types in 'cellclusters'. ")
      sel.cellClusters <- paste(unlist(strsplit(args$cellclusters, split = ', ')), sep = ', ' )
      pseudoRes <- getExpCondClusterPseudotime(rds = rds1, resDir = as.character(args$outputDir),
                                               newAnnotation = F,
                                               expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                               cellcluster = sel.cellClusters,
                                               slingshotclusterLabels = as.character(args$slingshotclusterLabels),
                                               topFeatureNo = as.numeric(args$topFeatureNo) )
    }
  } else {
    sel.expConds     <- paste(unlist(strsplit(args$expCond, split = ', ')), sep = ', ' )
    if (args$cellclusters == 'all') {
      print("Analysis is conducted on the selected experimental conditions specified in 'expCond' for all cell types. ")
      pseudoRes <- getExpCondClusterPseudotime(rds = rds1, resDir = as.character(args$outputDir),
                                               newAnnotation = F,
                                               expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                               expCond = sel.expConds,
                                               slingshotclusterLabels = as.character(args$slingshotclusterLabels),
                                               topFeatureNo = as.numeric(args$topFeatureNo) )
    } else {
      print("Analysis is conducted on the selected experimental conditions specified in 'expCond' for the selected cell types in 'cellclusters'. ")
      pseudoRes <- getExpCondClusterPseudotime(rds = rds1, resDir = as.character(args$outputDir),
                                               newAnnotation = F,
                                               expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                               expCond = sel.expConds,
                                               slingshotclusterLabels = as.character(args$slingshotclusterLabels),
                                               topFeatureNo = as.numeric(args$topFeatureNo),
                                               cellcluster = sel.cellClusters )
    }
  }
} else {
  if (args$expCond == 'all') {
    if (args$cellclusters == 'all') {
      print("Analysis is conducted on all experimental conditions specified in 'expCondCheck' for all cell types. ")
      pseudoRes <- getExpCondClusterPseudotime(rds = rds1, resDir = as.character(args$outputDir),
                                               newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                               expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                               slingshotclusterLabels = as.character(args$slingshotclusterLabels),
                                               topFeatureNo = as.numeric(args$topFeatureNo) )
    } else {
      print("Analysis is conducted on all experimental conditions specified in 'expCondCheck' the selected cell types in 'cellclusters'. ")
      sel.cellClusters <- paste(unlist(strsplit(args$cellclusters, split = ', ')), sep = ', ' )
      pseudoRes <- getExpCondClusterPseudotime(rds = rds1, resDir = as.character(args$outputDir),
                                               newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                               expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                               cellcluster = sel.cellClusters,
                                               slingshotclusterLabels = as.character(args$slingshotclusterLabels),
                                               topFeatureNo = as.numeric(args$topFeatureNo) )
    }
  } else {
    sel.expConds     <- paste(unlist(strsplit(args$expCond, split = ', ')), sep = ', ' )
    if (args$cellclusters == 'all') {
      print("Analysis is conducted on the selected experimental conditions specified in 'expCond' for all cell types. ")
      pseudoRes <- getExpCondClusterPseudotime(rds = rds1, resDir = as.character(args$outputDir),
                                               newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                               expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                               expCond = sel.expConds,
                                               slingshotclusterLabels = as.character(args$slingshotclusterLabels),
                                               topFeatureNo = as.numeric(args$topFeatureNo) )
    } else {
      print("Analysis is conducted on the selected experimental conditions specified in 'expCond' for the selected cell types in 'cellclusters'. ")
      sel.cellClusters <- paste(unlist(strsplit(args$cellclusters, split = ', ')), sep = ', ' )
      pseudoRes <- getExpCondClusterPseudotime(rds = rds1, resDir = as.character(args$outputDir),
                                               newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                               expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                               expCond = sel.expConds,
                                               slingshotclusterLabels = as.character(args$slingshotclusterLabels),
                                               topFeatureNo = as.numeric(args$topFeatureNo),
                                               cellcluster = sel.cellClusters )
    }
  }
}
print(sprintf("END pseudo-time trajectory analysis."))
save(pseudoRes, file = file.path(args$outputDir, sprintf("results_wNewAnnotation_pseudoTime/%s.Rdata", args$expCondCheckFname)))
## -------------------------------------------------------------------------------------##
print('END==========END==========END==========END==========END')
## -------------------------------------------------------------------------------------##
