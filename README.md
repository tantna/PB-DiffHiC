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
In the analysis of differential chromatin interactions, we only focused on intra-chromosomal interactions, excluding inter-chromosomal interactions. To speed up file reading and modeling, we chose to create separate files for each chromosome. Since the Hi-C contact matrix is symmetric, we only store the interaction information in the upper triangular part of the matrix, meaning for bin pairs (i, j), we store interactions where i ≤ j.

`PB-DiffHiC` has no specific requirements for input files, as long as R can read them. The input files must store information about the positions of interactions and the contact counts. The input file format can be either a dataframe or a matrix. If the file is a dataframe (.csv, .txt, etc.), each dataframe for a chromosome should contain four columns: (chr, x1, y1, ifs):
- `chr` - Chromosome of the interaction.
- `x1` - a start coordinate (in 10kb) of the first bin.
- `y1` - a start coordinate (in 10kb) of the second bin.
- `ifs` - the interaction frequency between bin pair (IFs).

```r
chr	 x1	 y1	ifs
chr19	1000	1002	7
chr19	1000	1003	6
chr19	1000	1004	2
chr19	1000	1005	2
chr19	1000	1006	1
chr19	1000	1007	2
chr19	1000	1008	2
```
You can also convert the Hi-C contact matrix data into a matrix and directly input it into `PB-DiffHiC`.

## Gaussian Convolution
The purpose of this step is to apply Gaussian convolution to enhance the interaction signals between interactions. After Gaussian convolution, the matrix is flattened into a vector, and interactions with a contact count of zero are also retained. The function to achieve this is implemented in `PBdata_preprocess.R`.
```r
source('PBdata_preprocess.R',encoding = 'UTF-8')
gauss2vec(data,test_dis,msize,keepdis,ksize,mu=0,sigma=1,data2mat=TRUE,gauss=TRUE)
```
- `data` - The input data containing the position and contact count information of interactions.
- `test_dis` - The gene distance range for interactions to be tested. For example, when the resolution is 10kb, if only interactions with a gene distance ≤ 1MB (1,000,000) are to be tested, then `test_dis` = 1000000 / 10000 + 1 = 101.
- `msize` - The size of the contact matrix. The chromosome length file is located in the `ext` folder. The size of the matrix is calculated by dividing the chromosome length by the resolution, rounding down, and then adding 1 to the result.
- `keepdis` - The size of the housekeeping gene distance. The housekeeping genome will be used for subsequent scaling factors calculations. Here, keepdis refers to how many diagonals to select. For example, `keepdis` = 1 refers to the main diagonal.
- `ksize` - The size of the Gaussian convolution kernel. In the PB-DiffHiC paper, `ksize` = 3 was used.
- `mu`, `sigma` - The mean and variance of the Gaussian function, with default values of 0 and 1.
- `data2mat` - Whether to convert the data into a matrix. If the input data is already a matrix, set `data2mat` = FALSE.
- `gauss` - Whether to perform Gaussian convolution. If Gaussian convolution is not needed, set `gauss` = FALSE.

The output obtained from the above function:

- `vec` - The contact count vector obtained after flattening the interactions within the testdis range, which will be used as input for the subsequent hypothesis testing.
- `h_keep` - The number of short-range interactions. In the vector vec, the range from 1 to h_keep is considered short-range interactions and is used to calculate the scaling factors.
- `from`, `to` - The positions of the interactions, equivalent to x1 and y1 in the input data.

## Hypothesis testing
A unified testing framework is provided under both the merged-replicate setup and the two-replicate setup.

After data processing, the hypothesis testing framework is used to detect differential chromatin interactions. Since there are differences in contact counts between conditions, scaling factors are first calculated to reduce the variation across conditions. Based on the computed scaling factors, hypothesis testing is then performed to obtain the P-value for each interaction.
```r
#merged-replicate setup
source('PB-DiffHiC_merged.R',encoding = 'UTF-8')
result=PB_merged(Datalist,Hkind,scale_factor,Scale=TRUE)
```
```r
#two-replicate setup
source('PB-DiffHiC_two.R',encoding = 'UTF-8')
result=PB_two(Datalist,Hkind,scale_factor,Scale=TRUE)
```
- `Datalist` - The data list obtained after Gaussian convolution, where the data from different conditions are stored in the list. The two-replicate setup is a general term, and the number of samples per condition can be ≥2.
- `Hkind` - The number of short-range interactions.
- `scale_factor` - The normalization factor. This parameter can be used to specify the value of the scale factor. In the merged-replicate setup, scale_factor is a single number (e.g., scale_factor = 1); in the two-replicate setup, it is a vector with the length equal to the total number of samples (e.g., in `PB-DiffHiC`'s two-replicate setup, you can set `scale_factor` = rep(1,4)).
- `Scale` - Whether to compute the scaling factors. If `Scale` = TRUE, the scale factor will be computed based on short-range interactions, and the specified `scale_factor` value will be ignored.

This function will perform hypothesis testing on all interactions in the input data list (datalist), and the results (such as P-values) will be output in the order of datalist. The final outputs are:
- `pv` - The P-value for each interaction.
- `qv` - The P-value for each interaction adjusted using the Benjamini-Hochberg method (BH).
- `scale` - The P-values computed using short-range interactions.

If you only need the scaling factors, you can directly call the function to compute the scaling factors.
```r
#merged-replicate setup
source('PB-DiffHiC_merged.R',encoding = 'UTF-8')
scale_result=cal_scale_merged(Datalist,Hkind)
```
```r
#two-replicate setup
source('PB-DiffHiC_two.R',encoding = 'UTF-8')
scale_result=cal_scale_merged(Datalist,Hkind)
```
The output will be the scaling factors, which can then be used in hypothesis testing by specifying the `scale_factor` parameter and setting the `Scale` parameter to FALSE.

# Example
In this example, we will use chromosome 19 from the mouse ESC and NPC dataset [(GSE210585)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210585) published by [Lee et al.](https://pubmed.ncbi.nlm.nih.gov/37649383/), with a resolution of 10kb. To generate pseudo-bulk Hi-C data, the `bin-step` processed single-cell Hi-C data from [SnapHiC-D](https://pubmed.ncbi.nlm.nih.gov/37649383/) will be combined to create pseudo-bulk Hi-C data. We will focus only on interactions with a gene distance within 1MB(`test_dis`=101), using the first five diagonals of the Hi-C contact matrix as short-range interactions(`keepdis`=5). The size of the Gaussian convolution kernel is set to 3(`ksize`=3). The example pseudo-bulk Hi-C data is located in the `example` folder.

First, we read the chromosome length file and prepare the data:
```r
binsize=10000
chrom_size=read.csv('ext/mouse_chromsize.txt',sep='\t',header=FALSE)
msize=chrom_size[19,2]%/%binsize+1
```
```r
#merged-replicate setup
filelist=c('example/two_samples/ESC_chr19_pseudo.csv','example/two_samples/ESC_chr19_pseudo.csv')
```
```r
#two-replicate setup
filelist=c('example/two_samples/ESC_chr19_s1.csv','example/two_samples/ESC_chr19_s2.csv',
        'example/two_samples/NPC_chr19_s1.csv','example/two_samples/NPC_chr19_s2.csv')
```
Then apply Gaussian convolution to the data and perform straightening and filling:
```r
datalist=list()
for (j in seq_along(filelist)){
  print(sprintf('file %s is doing!!!',j))
  df=read.csv(filelist[j])
  Hicvec=gauss2vec(df,101,msize,5,3,mu=0,sigma=1,sp2mat=TRUE,gauss=TRUE)
  datalist[[j]]=Hicvec$vec
}
HkeepCount=Hicvec$h_keep #the count of short-range interaction
```
Subsequently, calculate the scaling factors and perform hypothesis testing to obtain the p-value:
```r
#merged-replicate setup
source('PB-DiffHiC_merged.R',encoding = 'UTF-8')
result=PB_merged(datalist,HkeepCount,1,Scale=TRUE)
```
```r
#two-replicate setup
source('PB-DiffHiC_two.R',encoding = 'UTF-8')
result=PB_two(datalist,HkeepCount,rep(1,4),Scale=TRUE)
```
Finally, integrate the results:
```r
result_frame=data.frame(chr='chr19',x1=hicvec$from,y1=hicvec$to,pv=result$pv)
```
