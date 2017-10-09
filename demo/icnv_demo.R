library(iCNV)
library(truncnorm)
library(fields)
library(ggplot2)

##### IMPORTANT #####
# For detailed notation and pipeline, please visit https://github.com/zhouzilu/iCNV/blob/master/vignettes/iCNV-vignette.Rmd

####################################################
####################################################
##########                                ##########
##########           data input           ##########
##########                                ##########
####################################################
####################################################

# We load our input data. There are 9 list, each list have 10 samples
# sequencing plr and position: ngs_plr, ngs_plr.pos
# sequencing baf and position: ngs_baf, ngs_baf.pos
# snp array lrr and position: snp_lrr, snp_lrr.pos
# snp array snp and position: snp_baf, snp_baf.pos
# sample id: sample_id

load('demo_data.rda')
ls()
str(ngs_plr) # List of n vector, each one is the PLR for an exon
str(ngs_plr.pos) # List of n matrix (p x 2), each one is the start and end location for an exon
str(ngs_baf) # List of n vector, each one is the variants BAF from .bam
str(ngs_baf.pos) # List of n vector, each one is the variants BAF position
str(snp_lrr) # List of n vector, each one is the normalized LRR for a SNP
str(snp_lrr.pos) # List of n vector, each one is a SNP position
str(snp_baf) # List of n vector, each one is the BAF for a SNP
str(snp_baf.pos) # List of n vector, each one is the SNP BAF position

####################################################
####################################################
##########                                ##########
##########    iCNV without genotype       ##########
##########                                ##########
####################################################
####################################################

projname='icnv.demo'
icnv_res0=iCNV_detection(ngs_plr,snp_lrr,
                         ngs_baf,snp_baf,
                         ngs_plr.pos,snp_lrr.pos,
                         ngs_baf.pos,snp_baf.pos,
                         projname=projname,CN=0,mu=c(-3,0,2),cap=T,visual = 1)

####################################################
####################################################
##########                                ##########
##########       single sample plot       ##########
##########                                ##########
####################################################
####################################################

plotindi(ngs_plr,snp_lrr,
         ngs_baf,snp_baf,
         ngs_plr.pos,snp_lrr.pos,
         ngs_baf.pos,snp_baf.pos,
         icnv_res0,1)

####################################################
####################################################
##########                                ##########
##########        iCNV with genotype      ##########
##########                                ##########
####################################################
####################################################

projname='icnv.demo.geno'
icnv_res1=iCNV_detection(ngs_plr,snp_lrr,
                         ngs_baf,snp_baf,
                         ngs_plr.pos,snp_lrr.pos,
                         ngs_baf.pos,snp_baf.pos,
                         projname=projname,CN=1,mu=c(-3,0,2),cap=T,visual = 1)

####################################################
####################################################
##########                                ##########
########## iCNV with genotype (all plot)  ##########
##########                                ##########
####################################################
####################################################

projname='icnv.demo.geno.allplot'
icnv_res2=iCNV_detection(ngs_plr,snp_lrr,
                         ngs_baf,snp_baf,
                         ngs_plr.pos,snp_lrr.pos,
                         ngs_baf.pos,snp_baf.pos,
                         projname=projname,CN=1,mu=c(-3,0,2),cap=T,visual = 2)


####################################################
####################################################
##########                                ##########
##########     generate output tables     ##########
##########                                ##########
####################################################
####################################################

icnv.output = output_list(icnv_res2,sample_id,CN=1)
head(icnv.output)

####################################################
####################################################
##########                                ##########
##########       iCNV with NGS only       ##########
##########                                ##########
####################################################
####################################################

projname='icnv.demo.ngs'
icnv_res0=iCNV_detection(ngs_plr=ngs_plr,
                         ngs_baf=ngs_baf,
                         ngs_plr.pos=ngs_plr.pos,
                         ngs_baf.pos=ngs_baf.pos,
                         projname=projname,CN=0,mu=c(-3,0,2),cap=T,visual = 0)

projname='icnv.demo.geno.ngs'
icnv_res1=iCNV_detection(ngs_plr=ngs_plr,
                         ngs_baf=ngs_baf,
                         ngs_plr.pos=ngs_plr.pos,
                         ngs_baf.pos=ngs_baf.pos,
                         projname=projname,CN=1,mu=c(-3,0,2),cap=T,visual = 1)

projname='icnv.demo.geno.allplot.ngs'
icnv_res2=iCNV_detection(ngs_plr=ngs_plr,
                         ngs_baf=ngs_baf,
                         ngs_plr.pos=ngs_plr.pos,
                         ngs_baf.pos=ngs_baf.pos,
                         projname=projname,CN=1,mu=c(-3,0,2),cap=T,visual = 2)

icnv.output = output_list(icnv_res2,sample_id,CN=1)
head(icnv.output)

# Note for WGS
# You could easily generate BED file of all location using bed_generator function in iCNV/utils/bed_generator.R
bed_generator(chr=22,hg=19)
bed_generator(chr=22,hg=19,start=5001,end=10000,by=500)

####################################################
####################################################
##########                                ##########
##########       iCNV with SNP only       ##########
##########                                ##########
####################################################
####################################################

projname='icnv.demo.snp'
icnv_res0=iCNV_detection(snp_lrr=snp_lrr,
                         snp_baf=snp_baf,
                         snp_lrr.pos=snp_lrr.pos,
                         snp_baf.pos=snp_baf.pos,
                         projname=projname,CN=0,mu=c(-3,0,2),cap=T,visual = 0)


projname='icnv.demo.geno.snp'
icnv_res1=iCNV_detection(snp_lrr=snp_lrr,
                         snp_baf=snp_baf,
                         snp_lrr.pos=snp_lrr.pos,
                         snp_baf.pos=snp_baf.pos,
                         projname=projname,CN=1,mu=c(-3,0,2),cap=T,visual = 1)

projname='icnv.demo.geno.allplot.snp'
icnv_res2=iCNV_detection(snp_lrr=snp_lrr,
                         snp_baf=snp_baf,
                         snp_lrr.pos=snp_lrr.pos,
                         snp_baf.pos=snp_baf.pos,
                         projname=projname,CN=1,mu=c(-3,0,2),cap=T,visual = 2)

icnv.output = output_list(icnv_res2,sample_id,CN=1)
head(icnv.output)

