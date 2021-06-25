# scRICA: **s**ingle-**c**ell **R**NA-Seq **I**ntegrative **C**omparative **A**nalysis 

## what is scRICA
It is a R workflow package which can be used to perform downstream integrative, comparative analyses and visualization of a batch of scRNA-Seq count matrix from different experimental conditions efficiently and reproducibly. This package includes 2 types of fucntions: 1). workflow analysis functions and 2). visualization functions:

#### 1. Workflow analysis functions:
  
  * **`findDoublets()`**, this function is based on [DoubletDecon](https://github.com/EDePasquale/DoubletDecon) to perform either or both medoids and centroids deconvolution algorithms (https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31286-0) for doublet cells detection.
  * **`countReadin()`**, read the count matrix from different experiments defined in a metadata table into R for upstream quality control evaluation and mitochondrial content filtering process.
  * **`getClusterMarkers()`**, perform integrative analysis on all filtered cells from different experiments to identify different clustering cell types, and further identify conserved gene markers from each identified clustering cell types.
  * **`getExpCondClusterMarkers()`**, identify top expressed gene markers for samples from different experimental conditions with respect to each identified or annotated clustering cell types.
  * **`getClusterExpCondDe()`**, identify differential expressed gene markers for samples from 2 different experimental conditions with respect to each identified or annotated clustering cell types.
  * **`getClusterPseudo()`**, perform pseudotime trajectory functional analysis on the specified clustering cells with 3 different methods, including  principal components analysis (PCA) based on [Scater](https://bioconductor.org/packages/release/bioc/html/scater.html), Diffusion Maps based on [density](https://bioconductor.org/packages/release/bioc/html/destiny.html), and [slingshot](https://www.bioconductor.org/packages/release/bioc/html/slingshot.html).

#### 2. Visualization functions:

  * **`getClusterSummaryReplot()`**: summarize the number of cell in each identified or annotated clustering cell types, and generate corresponding tSNE and UMAP plots from different experimental conditions.
  * **`getGoiDotplot()`**: generate dot-plots of provided gene markers on the specified experimental conditions.
  * **`getGoiFeatureplot()`**: generate feature-plots of provided gene markers on the specified experimental conditions.
  * **`plotObjPseudotime()`**: generate pseudotime functional trajectory plots from 3 different analysis methods.

### Input metadata table file

As a workflow package, this package has an inherited structure, this package only replies on an intitial metadata table, where defines all samples names and corresponding count matrix results access locations, such as:

sample  | path  | doubletsRmMethod | doubletsResDir
------------- | -------------  | -------------  | ------------- 
sample1_condA_cond1  | /FullPath/to/CountMatrix/ | OL/centroids/medoids/none | /FullPath/to/result/dir/findDoublets()/
sample2_condA_cond2  | /FullPath/to/CountMatrix/ | OL/centroids/medoids/none | /FullPath/to/result/dir/findDoublets()/
sample3_condA_cond3  | /FullPath/to/CountMatrix/ | OL/centroids/medoids/none | /FullPath/to/result/dir/findDoublets()/
sample4_condB_cond1  | /FullPath/to/CountMatrix/ | OL/centroids/medoids/none | /FullPath/to/result/dir/findDoublets()/
sample5_condB_cond2  | /FullPath/to/CountMatrix/ | OL/centroids/medoids/none | /FullPath/to/result/dir/findDoublets()/

where if the column `'doubletsRmMethod' = 'none'`, there is no need to provide information for column `'doubletsResDir'`.

### Analysis workflow

#### 1. Step1: doublet cells detection
  
#### 2. Step2: integrative analysis
  
#### 3. Step3: gene marker identfication
  
#### 4. Step4: visualization

## Feedback
If you have further questions or suggestions regarding this pacakge, please contact Yan Li at yli22@bsd.uchicago.edu from the bioinformatics core at the Center for Research Informatics (CRI), biological science division (BSD), University of Chicago.
