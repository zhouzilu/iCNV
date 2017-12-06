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

## Update
**iCNV** has made a lot of changes on 10/31/2017 for stability, bug fixing and computation power. We strongly recommend you update iCNV to the newest version using the following command.
* Update instruction
```r
# Remove iCNV
remove.packages('iCNV')

# reinstall iCNV
install.packages("devtools")
library(devtools)
install_github("zhouzilu/iCNV")
```

## Workflow overview
Number in the parentheses referring to different section in [Vignettes](https://github.com/zhouzilu/iCNV/blob/master/vignettes/iCNV-vignette.Rmd) and function in the parenthesese can be found in [Utils](https://github.com/zhouzilu/iCNV/tree/master/utils).
```
                    NGS                                                 |                   Array
       BAM              BED (UCSC for WES or bed_generator.R for WGS)   |               SNP Intensity(in standard format)
        |----------------|                                              |                    |
        |                |                                              |                    |icnv_array_input.R (2.4)
        |SAMTools(2.3)   |CODEX(2.2)                                    |                    |
        |                |                                              |              |-----------|
Variants BAF(vcf)       PLR                                             |         Array LRR   Array BAF
        |                |                                              |              |           |
        |                |                                              |              |PCA(2.4)   |
        |                |                                              |              |           |
        |                |                                              |      Normalized LRR      |
        |                |                                              |              |           |
        --------------------------------------------------------------------------------------------
                                                        |
                                                        |iCNV_detection
                                                        |
                                                   CNV calling
```
## Demo code & Vignettes
* [Vignettes](https://github.com/zhouzilu/iCNV/blob/master/vignettes/iCNV-vignette.Rmd)
* [Demo](https://github.com/zhouzilu/iCNV/tree/master/demo)

## Utils
There are a couple of useful functions in [Utils](https://github.com/zhouzilu/iCNV/tree/master/utils). See each function for more information.
