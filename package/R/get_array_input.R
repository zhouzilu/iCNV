#' Get array information from given format
#' 
#' If your array input file follow the format in the example, you could use this function to extract array LRR and baf. Remember to load library before hands.
#' Save 4*[# of chr] lists, each list has N entry. N = # of individuals
#' snp_lrr: SNP LRR intensity; snp_lrr.pos: the position of the SNPs
#' snp_baf: the BAF of the SNPs; snp_baf.pos: the position of the SNPs
#' 
#' @param dir A string. The directory path to the folder where store signal intensity file according to chr
#' @param pattern A string. The pattern of all the intensity file
#' @param chr Specify the chromosome you want to generate. Must be of int from 1-22. If not specify, this function will generate files for all chromosomes.
#' @param projectname Name of the project
#' @return void 
#' @examples
#' load('icnv.demoarray_lrrbaf_22.rda')
#' str(snp_lrr)
#' str(snp_lrr.pos)
#' str(snp_baf)
#' str(snp_baf.pos)
#' \dontrun{
#' dir='PATH/TO/FOLDER'
#' pattern=paste0('*.csv.arrayicnv$')
#' icnv_array_intput(dir,pattern,chr=22,projectname='icnv.demo')
#' }
#' @export

get_array_intput=function(dir,pattern,chr=NULL,projectname=''){
	if(is.null(chr)){
		for (chr in 1:22){
			cat(chr,'\n')
			dirPath=paste0(dir,"/chr",chr)
			signalFile= list.files(dirPath, pattern = pattern)
			signaldir <- file.path(dirPath, signalFile)
			signal.DT <- lapply(signaldir, data.table::fread, sep=",")
			snp_lrr = sapply(signal.DT,function(t){t[[4]]},simplify = FALSE)
			names(snp_lrr)=signalFile
			snp_baf = sapply(signal.DT,function(t){t[[5]]},simplify = FALSE)
			names(snp_baf)=signalFile
			snp_lrr.pos = lapply(seq_len(length(snp_lrr)), function(i) signal.DT[[1]]$POS)
			snp_baf.pos = lapply(seq_len(length(snp_baf)), function(i) signal.DT[[1]]$POS)
			save(snp_lrr,snp_baf,snp_lrr.pos,snp_baf.pos, file=file.path(dirPath,paste0(projectname,'array_lrrbaf_',chr,'.rda')))
		}
	}else{
		cat(chr,'\n')
		dirPath=paste0(dir,"/chr",chr)
		signalFile= list.files(dirPath, pattern = pattern)
		signaldir <- file.path(dirPath, signalFile)
		signal.DT <- lapply(signaldir, data.table::fread, sep=",")
		snp_lrr = sapply(signal.DT,function(t){t[[4]]},simplify = FALSE)
		names(snp_lrr)=signalFile
		snp_baf = sapply(signal.DT,function(t){t[[5]]},simplify = FALSE)
		names(snp_baf)=signalFile
		snp_lrr.pos = lapply(seq_len(length(snp_lrr)), function(i) signal.DT[[1]]$POS)
		snp_baf.pos = lapply(seq_len(length(snp_baf)), function(i) signal.DT[[1]]$POS)
		save(snp_lrr,snp_baf,snp_lrr.pos,snp_baf.pos, file=file.path(dir,paste0(projectname,'array_lrrbaf_',chr,'.rda')))
	}

}


