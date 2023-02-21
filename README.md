# scRICA: **s**ingle-**c**ell **R**NA-Seq **I**ntegrative **C**omparative **A**nalysis 

### 1. What is scRICA
It is a systematic workflow R pacakge for integrative and comparative scRNA-seq analysis. Various parameter options are inputted from a metadata table and then inherited through the entire analysis workflow. scRICA can significantly improves efficiency of programming for comparative analyses with multiple sample attributes. It is distributed as an R package and also offers a command line execution option.

scRICA categorizes the workflow of integrative and comparative analysis using multiple-sample scRNA-seq data into four steps as shown below picture: ‘Step 1’ for pre-processing and quality control; ‘Step 2’ for multi-sample integration; ‘Step 3’ for visualization by attribute groups; and ‘Step 4’ for downstream analysis including differential expression (DE) analysis, pseudo-time trajectory analysis, and cell clusters identification.


### 2. pacakge installation and vignettes
  * Prerequisite packages
```r
devtools::install_github('EDePasquale/DoubletDecon')
devtools::install_github("ChenMengjie/lightHippo")
BiocManager::install("MAST")
BiocManager::install("scater")
BiocManager::install("destiny")
install.packages("mclust")
install.packages("easylabel")
```

  * github installation

```r
devtools::install_github(repo = 'yan-cri/scRICA', build_vignettes = T, force = T)
library(scRICA)
```
  
  * local download installation
     + Download package to your local computer to the local computer via `git clone https://github.com/yan-cri/scRICA.git`
     + Install downloaded scRICA package via:
```r
devtools::install('Path/to/Downloaded/scRICA', build_vignettes = T)
library(scRICA)
```
  
  * The package usage information can be seen from package vignettes via:
```r
browseVignettes(package = 'scRICA')
```

### 3. Demonstration data

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

### 4. Input metadata table file

This workflow package has its own inherited structure, and requires an initial metadata table to initiate the entire scRNA-Seq workflow analyses. 4 columns are required in the metadata table, they are `sample`, `path`, `doubletsRmMethod` and `expCond1` for samples names, full path to sample's count matrix tables directory (cellranger analysis results), doublets detective methods with four options, and experimental condition levels respectively; up to 2 experimental conditions specified in column `expCond1` and `expCond2` can be explored with this package. If no doublets removal is needed for samples, please specify 'none' for that sample in the column `doubletsRmMethod`.

sample  | path  | expCond1 | expCond2 | doubletsRmMethod
------------- | -------------  | -------------  | ------------- | ------------- 
sample1_condA_cond1  | /FullPath/to/CountMatrix/ | condA | cond1 | OL/centroids/medoids/none 
sample2_condA_cond2  | /FullPath/to/CountMatrix/ | condA | cond2 | OL/centroids/medoids/none 
sample3_condA_cond3  | /FullPath/to/CountMatrix/ | condA | cond3 | OL/centroids/medoids/none
sample4_condB_cond1  | /FullPath/to/CountMatrix/ | condB | cond1 | OL/centroids/medoids/none 
sample5_condB_cond2  | /FullPath/to/CountMatrix/ | condB | cond2 | OL/centroids/medoids/none 

To figure out the path where the package demonstration data are located, this can be find out with the command `system.file('extdata', package = 'scRICA', mustWork = T)` shown as below. An example of input metadata table can be find in the `'misc'` folder, where we include 4 package demonstration samples, 3396A, 3396F, 3041A, and 3041I. Please update the corresponding path in the metadata table column 'path'. 

```
print(system.file('extdata', package = 'scRICA', mustWork = T))
[1] "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/scRICA/extdata"
```

```
metadata <- read.delim2(file = 'path/to/metadata.txt', header = T) 
```

```
  sample expCond1 expCond2                                                                                 path doubletsRmMethod
1  A3396     3396        A /Library/Frameworks/R.framework/Versions/4.1/Resources/library/scRICA/extdata/3396A/               OL
2  F3396     3396        F /Library/Frameworks/R.framework/Versions/4.1/Resources/library/scRICA/extdata/3396F/               OL
3  A3041     3041        A  /Library/Frameworks/R.framework/Versions/4.1/Resources/library/scRICA/extdata/3041A             none
4  I3041     3041        I  /Library/Frameworks/R.framework/Versions/4.1/Resources/library/scRICA/extdata/3041I             None
```

### 5. Analysis workflow implementations

This package separate the scRNA-Seq integrative and comparative analysis into 3 categories, including integrative analysis with QC assessment for doublets and mitochondrial content, comparative analysis from different experimental condition groups, and visualizations.

#### 5.1 integrative analysis with QC assessment
 * QC assessment without mitochondrial content filtering
 * QC assessment without mitochondrial content filtering: turn on `mtFiltering` option and setup the mitochondrial content percentage cut-off values with option `mtPerCutoff`

```r
library(scRICA)
qcResult <- processQC(metadata = metadata, resDirName = 'scRICA_test_result', genomeSpecies = 'human')
qcResult <- processQC(metadata = metadata, resDirName = 'scRICA_test_result', genomeSpecies = 'human', mtFiltering = T, mtPerCutoff = 20)
```

```
## -- scRICA_test_result
##    |__doublet_results
##       |__ ...
##       |__doubletDecon_preProcessed_results
##          |__ ...
##    |__MT_percentage_summary.txt
##    |__No_filtered_cells_summary.txt
##    |__org_doubletsRemoval_cellNoSummary.txt
##    |__QC_plots
##       |__featureScatter_A3041.pdf
##       |__ ...
##       |__featureViolin_A3041.pdf
##       |__ ...
##       |__topVariableFeature_A3041.pdf
##       |__ ...
```

   The analyses results will be saved in the defined 'resDirName' option under the current working directory shown as below data structure, where includes 2 directories and 3 text files. 
   + Directory 'doublet_results' include the doublets detection results;
   + Directory 'QC_plots' contains all quality assessement results with respect to each sample in the metadata table;
   + Text file 'MT_percentage_summary.txt' presents a summary of mitochondrial content persentage for all samples;
   + Text file 'org_doubletsRemoval_cellNoSummary.txt' shows a summary of number of cells before and after doublets removal;
   + Text file 'No_filtered_cells_summary.txt' presents number of cells before and after mitochondrial content filtering, if otpion ` mtFiltering = T`.

#### 5.2 Integration analyses
```r
results <- getClusterMarkers(qcProcessedResults = qcResult)
```

```
## -- scRICA_test_result
##    |__allCluster_pos_markers_no.txt
##    |__allCluster_pos_markers_top10.txt
##    |__allCluster_pos_markers.txt
##    |__cluster_heatmap_top10PosMarkers.pdf
##    |__pca_plot.pdf
##    |__RDS_Dir
##       |__scRICA_test_result.rds
##    |__tsne_plot_samplSep.pdf
##    |__tsne_plot.pdf
##    |__umap_plot_samplSep.pdf
##    |__umap_plot.pdf
```

The integration results will be saved in the same dirctory defined in the `processQC()` execution. In addition to the previous `processQC()`, we can see that the integration analysis results are shown as below data structure.
 + Integration analysis results are saved in RDS file inside directory 'RDS_Dir';
 + Two types of clustering results are saved in file name starting with 'tsne*' and 'umap*' respectively;
 + 'allCluster_pos_markers.txt' includes the identified gene markers with resepct to each cell cluster;
 + 'allCluster_pos_markers_no.txt' summarize the number of identified gene markers in each cell cluster;
 + 'allCluster_pos_markers_top10.txt' presents the top N (defined in the option `topN` of `getClusterMarkers()`) identified gene markers with resepct to each cell cluster;
 + 'cluster_heatmap_top10PosMarkers.pdf' shows a heatmap of the top N (defined in the option `topN` of `getClusterMarkers()`) identified gene markers from each cell cluster.

#### 5.3 Integrated results summarization and visualization

The integrated analyses result can be summarized on the different experimental conditions with function `getClusterSummaryReplot()` shown as below, where user can define which experimental conditions to be summarized. If user would like to add cluster cells annotation, this can be doen by turn on the option `newAnnotation = T` and provide correpsonding Rscript name with option `newAnnotationRscriptName` wher the corresponding cluster cells annotations are defined.

```r
getClusterSummaryReplot(resDir = results$resDir, newAnnotation = F, expCondCheck = 'sample', expCondSepName = 'sample_org')
```

```
## -- scRICA_test_result
##    |__results_wOrgClusterAnnotation
##       |__expCond_sample_org
##          |__cellNo_summary_orgClusterAnnotation_expCond_sample_org.pdf
##          |__cellNo_summary_orgClusterAnnotation_expCond_sample_org.txt
##          |__new_tSNE_plot_expCond_sample_org
##             |__tsne_plot_noLabel_integrate_orgAnnotation.pdf
##             |__tsne_plot_wLabel_integrate_orgAnnotation.pdf
##             |__tsne_plot_wLabel_orgAnnotation_expCond_sample_org.pdf
##          |__new_UMAP_plot_expCond_sample_org
##             |__UMAP_plot_noLabel_integrate_orgAnnotation.pdf
##             |__UMAP_plot_wLabel_integrate_orgAnnotation.pdf
##             |__UMAP_plot_wLabel_orgAnnotation_expCond_sample_org.pdf
```
    
The option `expCondCheck` of `sample`, `expCond1`, or `expCond2` summarized the identified clusters cell numbers with respect to each sample defined in the metadata table column `sample`, experimental condition levels defined in the metadata table column `expCond1` or `expCond2` respectively. The summarized results are saved inside the directory 'results_wOrgClusterAnnotation' or 'results_wOrgClusterAnnotation' under the directory name of option `resDir`.

```r
getClusterSummaryReplot(resDir = results$resDir, newAnnotation = F, expCondCheck = 'expCond1', expCondSepName = 'expCond1_org' )
getClusterSummaryReplot(resDir = results$resDir, newAnnotation = F, expCondCheck = 'expCond2', expCondSepName = 'expCond2_org' )
```

```
## -- scRICA_test_result
##    |__results_wOrgClusterAnnotation
##       |__expCond_expCond1_org
##          |__ ...
##       |__expCond_expCond2_org
##          |__ ...
##       |__expCond_sample_org
##          |__ ...
```

#### 5.4 markers genes explorations

  * dot plot exploration
  
```r
getGoiDotplot(resDir = results$resDir, goiFname = 'scRICA/misc/marker_genes.xlsx', expCondCheck = 'expCond1', expCondSepName = 'expCond1_org')
```

```
## -- scRICA_test_result
##    |__results_wOrgClusterAnnotation
##       |__dotplot_selected_markers_expCond_expCond1_org
##          |__goiDotplots_markerGenes_dotplot_expCond_expCond1_org.pdf
```

  * feature plot exploration
```r
getGoiFeatureplot(resDir = results$resDir, goiFname = 'scRICA/misc/marker_genes.xlsx', expCondCheck  = 'expCond1')
```

```
## -- scRICA_test_result
##    |__results_wOrgClusterAnnotation
##       |__featurePlot_selected_markers_expCond1
##          |__goiFeaturePlot
##             |__umap_expMaxValSummary_expCondexpCond1.txt
##             |__umap_featurePlot_CAPS_expCondexpCond1.pdf
##             |__umap_featurePlot_CD79A_expCondexpCond1.pdf
##             |__ ...
```

```r   
getGoiFeatureplot(resDir = results$resDir, goiFname = 'scRICA/misc/marker_genes_set2.xlsx', expCondCheck  = 'expCond1', featurePlotFnamePrefix='goiFeaturePlot_set2')
```

```
## -- scRICA_test_result
##    |__results_wOrgClusterAnnotation
##       |__featurePlot_selected_markers_expCond1
##          |__goiFeaturePlot
##             |__umap_featurePlot_COL1A1_expCondexpCond1.pdf
##             |__umap_featurePlot_CSPG4_expCondexpCond1.pdf
##             |__ ...
##          |__goiFeaturePlot_set2
##             |__umap_expMaxValSummary_expCondexpCond1.txt
##             |__umap_featurePlot_CAPS_expCondexpCond1.pdf
##             |__ ...
```

## Feedback
If you have further questions or suggestions regarding this package, please contact Yan Li at yli22@bsd.uchicago.edu from the bioinformatics core at the Center for Research Bioinformatics (CRI), biological science division (BSD), University of Chicago.



