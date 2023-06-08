# scRICA: **s**ingle-**c**ell **R**NA-Seq **I**ntegrative **C**omparative **A**nalysis 

### 1. What is scRICA
It is a systematic workflow R package for integrative and comparative scRNA-seq analysis. Various parameter options are inputted from a metadata table and then inherited through the entire analysis workflow. scRICA can significantly improves efficiency of programming for comparative analyses with multiple sample attributes. It is distributed as an R package and also offers a command line execution option.

scRICA categorizes the workflow of integrative and comparative analysis using multiple-sample scRNA-seq data into four steps as below:
+ Step 1: Pre-processing and quality control;
+ Step 2: Multi-sample integration;
+ Step 3: Visualization by attribute groups;
+ Step 4: Downstream analysis including differential expression (DE) analysis, pseudo-time trajectory analysis, and cell clusters identification.

### 2. pacakge installation and vignettes
  * Prerequisite packages
```r
BiocManager::install("MAST")
BiocManager::install("scater")
BiocManager::install("destiny")
install.packages("mclust")
install.packages("easylabel")
install.packages("corrplot")
devtools::install_github('EDePasquale/DoubletDecon')
devtools::install_github("ChenMengjie/lightHippo")
BiocManager::install("slingshot")
BiocManager::install("DESeq2")
BiocManager::install("SingleCellExperiment")
BiocManager::install("tradeSeq")
```

  * github installation

```r
devtools::install_github(repo = 'yan-cri/scRICA', build_vignettes = F, force = T)
library(scRICA)
```
  
  * local download installation
     + Download package to your local computer to the local computer via `git clone https://github.com/yan-cri/scRICA.git`
     + Install downloaded scRICA package via:
```r
devtools::install('Path/to/Downloaded/scRICA', build_vignettes = F)
library(scRICA)
```

### 3. Demonstration data
#### 3.1 Demonstration data for multi-sample pre-processing and integration
As shown in below, a total of 6 samples single cell gene expression count matrix results from cell ranger analysis can be downloaded from https://drive.google.com/drive/folders/1wwxmi4cvASRDFtiFKs6MbflFDOOxw3AY?usp=sharing using as the demonstration data for this package analysis implementations. 

```
-- 3041A
   |__barcodes.tsv
   |__genes.tsv
   |__matrix.mtx
-- 3041F
   |__barcodes.tsv
   |__genes.tsv
   |__matrix.mtx
-- 3041I
   |__barcodes.tsv
   |__genes.tsv
   |__matrix.mtx
-- 3396A
   |__barcodes.tsv.gz
   |__features.tsv.gz
   |__matrix.mtx.gz
-- 3396F
   |__barcodes.tsv.gz
   |__features.tsv.gz
   |__matrix.mtx.gz
-- 3396I
   |__barcodes.tsv.gz
   |__features.tsv.gz
   |__matrix.mtx.gz
```

#### 3.2 Demonstration data example of an integrated RDS

This RDS file can be downloaded from https://drive.google.com/drive/folders/1wwxmi4cvASRDFtiFKs6MbflFDOOxw3AY?usp=sharing, which can be used as demonstration data for scRICA step 3 visualization and step 4 downstream analysis implementations.

This RDS includes an integration study from our human cell atlas study of Fallopian Tubes, where a total of 18 scRNA-seq sequencing samples from 8 donors fallopian anatomical sites, including isthmus (I), ampule (A) and fambiriae (F). As shown in the table head, these experimental attributes, including donor and anatomical sites, were defined in the original metadata table column ’’ and ’’ which has been inherited into the integrated RDS object column ’’ and ’’ respectively.


### 4. Usage implementation
The implementation vignette can be seen at https://rpubs.com/yli_cri/1025790.

## Feedback
If you have further questions or suggestions regarding this package, please contact Yan Li at yli22@bsd.uchicago.edu from the bioinformatics core at the Center for Research Bioinformatics (CRI), biological science division (BSD), University of Chicago.



