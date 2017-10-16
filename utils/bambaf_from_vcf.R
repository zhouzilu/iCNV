#' If your vcf follow the format in the example, you could use this function to extract NGS baf from vcf files. Remember to load library before hands.
#' Save 6 lists, each list has N entry. N = # of individuals (or vcf file)
#' bambafnm: name of the bamfiles; bambafchrls: the chromosome; bambafposls: the position of the variants; 
#' bambafls: the BAF of the variants; bambafidls: the ID of the variants; filenm:the file name
#' @param vcf_list All the vcf names stored in vcf.list; could use command:"ls *.vcf > vcf.list" to generate.
#' @param chr Specify the chromosome you want to generate. Must be of int from 1-22. If not specify, this function will generate all chromosomes.
#' @return void
#' @examples
#' setwd('PATH/TO/FOLDER')
#' library(tidyr);library(data.table);library(dplyr)
#' bambaf_from_vcf('example_vcf.list')
#' bambaf_from_vcf('example_vcf.list',chr=22)
#' load('bambaf_22.rda')
#' ngs_baf = bambafls
#' ngs_baf.pos = bambafposls
#' @export
bambaf_from_vcf = function(vcf_list,chr=NULL){

	filenm=read.table(vcf_list,header=F,as.is=T)[[1]]

	# Read all the VCF data
	baf.all=sapply(filenm,bam.baf,simplify=F)

	# Seperate the data into different chromosome
	if(is.null(chr)){
		for (chr in 1:22){
			baf.all.chr=sapply(baf.all,function(x){return(x %>% filter(`#CHROM`==chr))},simplify=F)
			bambafnm=list()
			bambafchrls=list()
			bambafposls=list()
			bambafidls=list()
			bambafls=list()
			for (i in seq(1:length(baf.all.chr))){
				x=baf.all.chr[[i]]
				bambafnm[[i]]=names(x)[4]
				bambafchrls[[i]]=x$`#CHROM`
				bambafposls[[i]]=x$POS
				bambafidls[[i]]=x$ID
				bambafls[[i]]=x[[4]]
			}
			save(bambafnm,bambafchrls,bambafposls,bambafls,bambafidls,filenm,file=paste0('bambaf_',chr,'.rda'))
		}
	}else{
		baf.all.chr=sapply(baf.all,function(x){return(x %>% filter(`#CHROM`==chr))},simplify=F)
		bambafnm=list()
		bambafchrls=list()
		bambafposls=list()
		bambafidls=list()
		bambafls=list()
		for (i in seq(1:length(baf.all.chr))){
			x=baf.all.chr[[i]]
			bambafnm[[i]]=names(x)[4]
			bambafchrls[[i]]=x$`#CHROM`
			bambafposls[[i]]=x$POS
			bambafidls[[i]]=x$ID
			bambafls[[i]]=x[[4]]
		}
		save(bambafnm,bambafchrls,bambafposls,bambafls,bambafidls,filenm,file=paste0('bambaf_',chr,'.rda'))
	}
}

bam.baf=function(filenm){
	cat('load ',filenm,'...\n')
	dt = fread(filenm,skip='#CHROM')
	nm=names(dt)[length(dt)]
	cat('individual ',nm,'\n')
	dt <- dt %>% filter(QUAL!='.') %>% select(-REF,-ALT,-QUAL,-FILTER,-FORMAT,-INFO)
	bafi <- sapply(strsplit(dt[[4]],':'),function(x){a=as.numeric(strsplit(x[2],',')[[1]]);return(a[1]/(a[1]+a[2]))})
	dt[nm] <- bafi
	return(dt)
}
