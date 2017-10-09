# iCNV
Integrated copy number variation detection toolset

## Author
Zilu Zhou, Nancy R. Zhang

## Maintainer
Zilu Zhou <zhouzilu@upenn.edu>
Please comment on the *Issues* section for addtional questions.

## Description
**iCNV** is a normalization and copy number variation detection procedure for mutiple study designs: WES only, WGS only, SNP array only, or any combination of SNP and sequencing data. **iCNV** applies platform specific normalization, utilizes allele specific reads from sequencing and integrates matched NGS and SNP-array data by a Hidden Markov Model (HMM).

## Installation
* Install from GitHub
```r
# Install dependent packages first
install.packages("fields")
install.packages("truncnorm")
install.packages("ggplot2")
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("CODEX")

# Install iCNV
install.packages("devtools")
library(devtools)
install_github("zhouzilu/iCNV")
```

## Demo code & Vignettes
* [Vignettes](https://github.com/zhouzilu/iCNV/blob/master/vignettes/iCNV-vignette.Rmd)
* [Demo](https://github.com/zhouzilu/iCNV/tree/master/demo)

## Utils
There are a couple of useful functions in [Utils](https://github.com/zhouzilu/iCNV/tree/master/utils). See each function for more information.
