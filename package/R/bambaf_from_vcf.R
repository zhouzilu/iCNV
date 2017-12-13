#' If your vcf follow the format in the example, you could use this function to extract NGS baf from vcf files. Remember to load library before hands.
#' Save 6 lists, each list has N entry. N = # of individuals (or vcf file)
#' ngs_baf.nm: name of the bamfiles; ngs_baf.chr: the chromosome; ngs_baf.pos: the position of the variants; 
#' ngs_baf: the BAF of the variants; ngs_baf.id: the ID of the variants; filenm:the file name
#' @param dir The directory to all the vcf stored; default is right in this folder.
#' @param vcf_list All the vcf names stored in vcf.list; could use command:"ls *.vcf > vcf.list" to generate.
#' @param chr Specify the chromosome you want to generate. Must be of int from 1-22. If not specify, this function will generate all chromosomes.
#' @return void 
#' @examples
#' dir='PATH/TO/FOLDER'
#' bambaf_from_vcf(dir,'example_vcf.list')
#' bambaf_from_vcf(dir,'example_vcf.list',chr=22)
#' load('bambaf_22.rda')
#' str(ngs_baf)
#' str(ngs_baf.pos)
#' @export
bambaf_from_vcf = function(dir='.',vcf_list,chr=NULL){
	`%>%`=tidyr::`%>%`
	filenm=file.path(dir,read.table(file.path(dir,vcf_list),header=F,as.is=T)[[1]])

	# Read all the VCF data
	baf.all=sapply(filenm,bam.baf,simplify=F)

	# Seperate the data into different chromosome
	if(is.null(chr)){
		for (chr in 1:22){
			baf.all.chr=sapply(baf.all,function(x){x1=x %>% dplyr::filter(`#CHROM`==chr); return(x1)},simplify=F)
			ngs_baf.nm=list()
			ngs_baf.chr=list()
			ngs_baf.pos=list()
			ngs_baf.id=list()
			ngs_baf=list()
			for (i in seq(1:length(baf.all.chr))){
				x=baf.all.chr[[i]]
				ngs_baf.nm[[i]]=names(x)[4]
				ngs_baf.chr[[i]]=x$`#CHROM`
				ngs_baf.pos[[i]]=x$POS
				ngs_baf.id[[i]]=x$ID
				ngs_baf[[i]]=x[[4]]
			}
			save(ngs_baf.nm,ngs_baf.chr,ngs_baf.pos,ngs_baf,ngs_baf.id,filenm,file=file.path(dir,paste0('bambaf_',chr,'.rda')))
		}
	}else{
		baf.all.chr=sapply(baf.all,function(x){return(x %>% dplyr::filter(`#CHROM`==chr))},simplify=F)
		ngs_baf.nm=list()
		ngs_baf.chr=list()
		ngs_baf.pos=list()
		ngs_baf.id=list()
		ngs_baf=list()
		for (i in seq(1:length(baf.all.chr))){
			x=baf.all.chr[[i]]
			ngs_baf.nm[[i]]=names(x)[4]
			ngs_baf.chr[[i]]=x$`#CHROM`
			ngs_baf.pos[[i]]=x$POS
			ngs_baf.id[[i]]=x$ID
			ngs_baf[[i]]=x[[4]]
		}
		save(ngs_baf.nm,ngs_baf.chr,ngs_baf.pos,ngs_baf,ngs_baf.id,filenm,file=file.path(dir,paste0('bambaf_',chr,'.rda')))
	}
}

bam.baf=function(filenm){
	`%>%`=tidyr::`%>%`
	cat('load ',filenm,'...\n')
	dt = data.table::fread(filenm,skip='#CHROM')
	nm=names(dt)[length(dt)]
	cat('individual ',nm,'\n')
	dt <- dt %>% dplyr::filter(QUAL!='.') %>% dplyr::select(-REF,-ALT,-QUAL,-FILTER,-FORMAT,-INFO)
	bafi <- sapply(strsplit(dt[[4]],':'),function(x){a=as.numeric(strsplit(x[2],',')[[1]]);return(a[1]/(a[1]+a[2]))})
	dt[nm] <- bafi
	return(dt)
}
