packages <- c("easylabel")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

BCpackages <- c("MAST", "scater", "destiny", "slingshot", "SingleCellExperiment")
if (length(setdiff(BCpackages, rownames(installed.packages()))) > 0) {
  # source("http://bioconductor.org/biocLite.R")
  # biocLite(setdiff(BCpackages, rownames(installed.packages())))
  BiocManager::install(setdiff(BCpackages, rownames(installed.packages())))
}

## install below github packages manually
githubPacakges <- c("lightHippo", "DoubletDecon")
## devtools::install_github("ChenMengjie/lightHippo")
## devtools::install_github('EDePasquale/DoubletDecon')

sapply(c(packages, BCpackages, githubPacakges), require, character.only=T)

print(sapply(c(packages, BCpackages, githubPacakges), require, character.only=T))
