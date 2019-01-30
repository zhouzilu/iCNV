# iCNV
Integrated copy number variation detection toolset

## Author
Zilu Zhou, Nancy R. Zhang

## Maintainer
Zilu Zhou <zhouzilu@pennmedicine.upenn.edu>
Please comment on the *Issues* section for addtional questions.

## Description
**iCNV** is a normalization and copy number variation detection procedure for mutiple study designs: WES only, WGS only, SNP array only, or any combination of SNP and sequencing data. **iCNV** applies platform specific normalization, utilizes allele specific reads from sequencing and integrates matched NGS and SNP-array data by a Hidden Markov Model (HMM).

## Installation
* Install from GitHub (Recommended)
```r
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("CODEX")

# Install iCNV
install.packages("devtools")
library(devtools)
install_github("zhouzilu/iCNV")
```

## Update
**iCNV** has made a lot of changes on 1/29/2019 for stability, bug fixing, computation power and functionality. We strongly recommend you update iCNV to the newest version using the following command.
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
Number in the parentheses referring to different section in [Vignettes](http://rawgit.com/zhouzilu/iCNV/master/vignettes/iCNV-vignette.html) and function details can be found [here](https://bioconductor.org/packages/release/bioc/manuals/iCNV/man/iCNV.pdf)
```
        NGS                                             |           Array
BAM    BED(UCSC for WES or bed_generator.R for WGS 2.2) |    SNP Intensity(in standard format)
 |                |                                     |             |
 |----------------|                                     |             |
 |                |                                     |             |icnv_array_input (2.4)
 |SAMTools(2.3)   |CODEX(2.2)                           |             |
 |                |                                     |             |-----------|
Variants BAF(vcf) PLR                                   |        Array LRR   Array BAF
 |                |                                     |             |           |
 |                |                                     |             |SVD(2.4)   |
 |                |                                     |             |           |
 |                |                                     |     Normalized LRR      |
 |                |                                     |             |           |
 -----------------------------------------------------------------------------------
                                          |
                                          |iCNV_detection(2.5-2.6)
                                          |
                                     CNV calling
                                          |
                                          |icnv_output_to_gb()
                                          |
                              Genome Browser input
```
## Demo code & Vignettes with rich displays
* [Vignettes](http://rawgit.com/zhouzilu/iCNV/master/vignettes/iCNV-vignette.html)
* [Vignettes.RMD](https://github.com/zhouzilu/iCNV/blob/master/vignettes/iCNV-vignette.Rmd)
* [Function manual](https://bioconductor.org/packages/release/bioc/manuals/iCNV/man/iCNV.pdf)

## Citation
Zilu Zhou, Weixin Wang, Li-San Wang, Nancy Ruonan Zhang; Integrative DNA copy number detection and genotyping from sequencing and array-based platforms, Bioinformatics, Volume 34, Issue 14, 15 July 2018, Pages 2349–2355, https://doi.org/10.1093/bioinformatics/bty104
