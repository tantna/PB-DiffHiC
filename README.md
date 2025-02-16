# Introduction
`PB-DiffHiC` is a new optimized parametric statistical framework that directly analyzes the non-imputed pseudo-bulk Hi-C data at 10 Kb resolution. PB-DiffHiC incorporates Gaussian convolution, stability of short-range interactions, and Poisson distribution to enable joint normalization and detection of significant differential chromatin interactions between conditions.

`PB-DiffHiC` provides a unified framework for two primary setups and consists of two key steps:
- Applying Gaussian convolution to enhance interaction signals in pseudo-bulk Hi-C data.
- Optimizing parametric hypothesis testing by combining P-value calculation with the estimation of scaling factors for normalizing contact matrices across conditions.

This project includes three main R scripts:
- `PBdata_preprocess.R`: Applies Gaussian convolution to the raw pseudo-bulk Hi-C data and flattens interactions within the selected range.
- `PB-DiffHiC_merged.R`: Detects Differential Chromatin Interactions under the merged-replicate setup.
- `PB-DiffHiC_merged.R`: Detects Differential Chromatin Interactions under the two-replicate setup.
  
# Installation
Running `PB-DiffHiC` requires the following dependency packages：
```r
library(SCBN)
library(Matrix)
library(mvtnorm)
library(edgeR)
library(statmod)
```
[SCBN](https://bioconductor.org/packages/release/bioc/html/SCBN.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) click the link to see how to download. Other dependent packages can be downloaded directly using `install.packages()`.

# Usage
## Input
