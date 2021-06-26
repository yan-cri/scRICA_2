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

### Analysis workflow implementation example

The entire scRNA-Seq integrative and comparative analysis can be implemented with below workflow steps based on the input metadata table and other associated options for each step's analysis. User can skip doublet detection step starting from step2 to perform integration analysis directly by setting up metadata table column `'doubletsRmMethod'` as `none`, otherwise, start the analysis from step1 `findDoublets()` with first 2 columns of above metadata table information. 

We provided an ovrian caner study of 1 patient samples from 3 different fallopian tissue types as this package demonstration dataset, which can be accessible at package 'inst/extdata' folder with command `list.files(system.file('extdata', package = 'scRICA', mustWork = T))`. For the demonstration, we also include the corresponding metadata table for this experiment accessible at `system.file('extdata', 'metadata.txt', package = 'scRICA', mustWork = T)`. Please be aware column `path` information in this demonstration metadata table is based on the Mac system, if you are using other operating systerm, pelase make your own metadata table.

#### 1. Step1: doublet cells detection: 
```{r}
metadata                <- read.delim2(file = system.file('extdata', 'metadata.txt', package = 'scRICA', mustWork = T), header = T) 
doubletDeconResDir      <- findDoublets(metadata = metadata, 
                                        genomeSpecies = 'human', 
                                        doubletDeconRhop =  0.5, 
                                        doubletDeconPMF = 'F', 
                                        resFilename = 'scRICA_test_doublets_checking')
```
This step will conduct cells doublets detection and save the entire analysis results into the defined `resFilename` under current operating directory `getwd()`, it will return a full path of directory to indict where the results were saved at. 

#### 2. Step2: quality control and integrative analysis

This steps analysis will use 2 functions: `countReadin()` and `getClusterMarkers()`.

Based on the previous step's analysis, we need to update our metadata data by including column `doubletsRmMethod` and `doubletsResDir` as shown in above metadata table input file section, where `doubletsRmMethod` allow users to specify which doublet cells detection analysis method to be used for doublet cells filtering, 3 extra options are provided here, 'medoids', 'centroids', or 'OL' (doublet cells detected by both medoids and centroids deconvolution algorithms); `doubletsResDir` specify where the previous analysis results are located. The metadata table can be updated as below and used for this step's integrative analysis via function `countReadin()`.

```{r}
metadataUpdate                  <- metadata
metadataUpdate$doubletsRmMethod <- rep('OL', length(metadataUpdate$sample))
metadataUpdate$doubletsResDir   <- rep(as.character(doubletDeconResDir), length(metadataUpdate$sample))
## ---
seuratProcssedObjList           <- countReadin(metadata = metadataUpdate, 
                                               resDirName = 'scRICA_test', 
                                               genomeSpecies = 'human', 
                                               minCells = 3, minFeatures = 200, 
                                               mtFiltering = 'T', mtPerCutoff = 20)
```

This function will conduct integrative analysis for all samples provided in the metadata table. If user would like to filter out the mitochondrial content, as shown here, please specify `mtFiltering = 'T'` together with the mitochondrial content filtering percentage with option `mtPerCutoff`. This function will 1) prepare a Seurat object by reading the count matrix from different samples specified in the metadata table into R, and 2) perform the quality control evaluations. 

```{r}
seuratIntegratedRes             <- getClusterMarkers(qcProcessedSeuratObjList = seuratProcssedObjList$qcProcessObj, 
                                                     anchorIntegrate = T, 
                                                     resDirName = 'scRICA_test' )
                                             
```
In addition to integrative analysis for all samples, this function will also conduct the top expressed markers detection with respect to each cell clusters. The entire analysis results will be saved in the specified `resDirName`, the final integrated RDS file will be saved in `'RDS_Dir'` folder inside `resDirName` directory. If no `resDirName` is specified, a folder named as `integration_analysis_results` will be created to save this function's analysis results.

#### 3. Step3: gene marker identfication
2 types of gene markers identification with respect to experimental conditions can be performed with corresponding functions `getClusterExpCondDe()` and `getClusterExpCondDe()`, where as specified above, `getExpCondClusterMarkers()` can identify top expressed gene markers for samples from different experimental conditions with respect to each identified or annotated clustering cell types; and `getClusterExpCondDe()` can identify differential expressed gene markers for samples from 2 different experimental conditions with respect to each identified or annotated clustering cell types.

#### 4. Step4: clustering and gene marker visualizations
The integrative cell clustering analysis results can be visualized via function `getClusterSummaryReplot()`, and additional gene markers exploration can be done via function `getGoiDotplot()` and `getGoiFeatureplot()`.

#### 5. Step5: pseudotime functional trajectory analysis and visualization

Pseudotime functional trajectory analysis can be performed at different specfied experimental cell clusers via function `getClusterPseudo()`, and the corresponding visualization can be performed via function `plotObjPseudotime()`.

## Feedback
If you have further questions or suggestions regarding this pacakge, please contact Yan Li at yli22@bsd.uchicago.edu from the bioinformatics core at the Center for Research Informatics (CRI), biological science division (BSD), University of Chicago.
