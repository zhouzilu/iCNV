library(iCNV)
library(truncnorm)
library(fields)
library(ggplot2)

##### IMPORTANT #####
# For detailed notation and pipeline, please visit https://github.com/zhouzilu/iCNV


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
str(ngs_plr)
str(ngs_plr.pos)
str(ngs_baf)
str(ngs_baf.pos)
str(snp_lrr)
str(snp_lrr.pos)
str(snp_baf)
str(snp_baf.pos)

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
                         projname=projname,CN=0,mu=c(-3,0,2),cap=T,visual = 2)

####################################################
####################################################
##########                                ##########
##########       single sample plot       ##########
##########                                ##########
####################################################
####################################################

icnv_call = icnv_res0[[1]]
plotindi(ngs_plr,snp_lrr,
         ngs_baf,snp_baf,
         ngs_plr.pos,snp_lrr.pos,
         ngs_baf.pos,snp_baf.pos,
         icnv_call,1)

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
                         projname=projname,CN=1,mu=c(-3,0,2),cap=T,visual = 2)

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
                         projname=projname,CN=1,mu=c(-3,0,2),cap=T,visual = 1)


####################################################
####################################################
##########                                ##########
##########     generate output tables     ##########
##########                                ##########
####################################################
####################################################

icnv.output = output_list(icnv_res2,sample_id,CN=1)
head(icnv.output)
