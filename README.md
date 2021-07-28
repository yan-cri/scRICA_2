# scRICA: **s**ingle-**c**ell **R**NA-Seq **I**ntegrative **C**omparative **A**nalysis 

### 1. What is scRICA
It is a R workflow package which can be used to perform scRNA-Seq downstream integrative, comparative analyses and visualization. This package can process a batch of scRNA-Seq count matrix from different experimental conditions fro integration and comparative analysis efficiently and reproducible. The package functions can be categorized into: 1). analysis workflow functions and 2). visualization functions:

### 2. pacakge installation and vignettes

  * github installation
```
devtools::install_github(repo = 'yan-cri/scRICA', build_vignettes = T, force = T)
library(scRICA)
```
  
  * local download installation
     + Download package to your local computer to the local computer via `git clone https://github.com/yan-cri/scRICA.git`
     + Install downloaded scRICA package via:
```
devtools::install('Path/to/Downloaded/scRICA', build_vignettes = T)
library(scRICA)
```
  
  * The package usage information can be seen from package vignettes via:
```
browseVignettes(package = 'scRICA')
```

## 3. Feedback
If you have further questions or suggestions regarding this package, please contact Yan Li at yli22@bsd.uchicago.edu from the bioinformatics core at the Center for Research Bioinformatics (CRI), biological science division (BSD), University of Chicago.



