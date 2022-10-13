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

parser$add_argument("--selected.DB", type="character",
                    default = NULL,
                    help="Conduct cellChat analysis on the selected pathways with 3 options, they are 'Cell-Cell Contact', 'ECM-Receptor' and 'Secreted Signaling'.",
                    metavar='')

parser$add_argument("--expCond", type="character",
                    default = 'all',
                    help = "Selecte levels of samples attributes to conducte cellchat analysis, should be consistent with the levels in defined in the 'expCondCheck'.",
                    metavar='')

parser$add_argument("--min.cells.group", type="integer",
                    default = '10',
                    help="The minimum number of cells presenting in each comparison groups.",
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
print(sprintf("START cellChat anlaysis "))
## ------
print("START read in RDS file")
rds1             <- readRDS(as.character(args$rdsFname))
setwd(as.character(args$outputDir))
print("Complete read in RDS file")
print('-=-=-=-=-=-=-=-=-=-')
## ---
if (is.null(args$cellAnnotation)) {
  if (is.null(args$selected.DB)) {
    if (args$expCond == 'all') {
      cellChatResList <- getCellChatRes(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = F,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        min.cells.filtering = as.numeric(args$min.cells.group),
                                        debug = debug)
    } else {
      sel.expConds     <- paste(unlist(strsplit(args$expCond, split = ', ')), sep = ', ' )
      cellChatResList <- getCellChatRes(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = F,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        expCond = sel.expConds,
                                        min.cells.filtering = as.numeric(args$min.cells.group),
                                        debug = debug)
    }

  } else {
    if (args$expCond == 'all') {
      cellChatResList <- getCellChatRes(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = F,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        selected.DB = as.character(args$selected.DB),
                                        min.cells.filtering = as.numeric(args$min.cells.group),
                                        debug = debug)
    }else {
      sel.expConds     <- paste(unlist(strsplit(args$expCond, split = ', ')), sep = ', ' )
      cellChatResList <- getCellChatRes(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = F,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        expCond = sel.expConds,
                                        selected.DB = as.character(args$selected.DB),
                                        min.cells.filtering = as.numeric(args$min.cells.group),
                                        debug = debug)
    }

  }
  print(sprintf("END cellChat anlaysis "))
  save(cellChatResList, file = file.path(args$outputDir, sprintf("results_wOrgClusterAnnotation_cellChat/%s_full.Rdata", args$expCondCheckFname)))
} else {
  if (is.null(args$selected.DB)) {
    if (args$expCond == 'all') {
      cellChatResList <- getCellChatRes(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        min.cells.filtering = as.numeric(args$min.cells.group),
                                        debug = debug)
    }else {
      sel.expConds     <- paste(unlist(strsplit(args$expCond, split = ', ')), sep = ', ' )
      cellChatResList <- getCellChatRes(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        expCond = sel.expConds,
                                        min.cells.filtering = as.numeric(args$min.cells.group),
                                        debug = debug)
    }

  } else {
    if (args$expCond == 'all') {
      cellChatResList <- getCellChatRes(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        selected.DB = as.character(args$selected.DB),
                                        min.cells.filtering = as.numeric(args$min.cells.group),
                                        debug = debug)
    }else {
      sel.expConds     <- paste(unlist(strsplit(args$expCond, split = ', ')), sep = ', ' )
      cellChatResList <- getCellChatRes(rds = rds1, resDir = as.character(args$outputDir),
                                        newAnnotation = T, newAnnotationRscriptName = args$cellAnnotation,
                                        expCondCheck = as.character(args$expCondCheck), expCondCheckFname = as.character(args$expCondCheckFname),
                                        expCond = sel.expConds,
                                        selected.DB = as.character(args$selected.DB),
                                        min.cells.filtering = as.numeric(args$min.cells.group),
                                        debug = debug)
    }

  }
  print(sprintf("END cellChat anlaysis "))
  save(cellChatResList, file = file.path(args$outputDir, sprintf("results_wNewAnnotation_cellChat/%s_full.Rdata", args$expCondCheckFname)))
}
## -------------------------------------------------------------------------------------##
print('END==========END==========END==========END==========END')
## -------------------------------------------------------------------------------------##
