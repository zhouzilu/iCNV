#' Generate BED file for WGS dataset. 
#' 
#' Default position generated from USCS genome browser
#' 
#' @param chr Specify the chromosome you want to generate. Must be of int from 1-22
#' @param hg Specify the coordinate you want to generate from. Start and end position of hg19 and hg38 have been pre-implemented.
#' @param start The start position of your BED file.
#' @param end The end position of your BED file.
#' @param by The chunk of your DNA for each bin. Default 1kb.
#' @return void
#' @examples
#' bed_generator(chr=22,hg=38)
#' bed_generator(22,38,5001,10000,by=500)
#' @export
bed_generator <- function(chr, hg, start=NULL, end=NULL, by=1000){
	hg19_loc = matrix(c(0,249251000,1,243200000,1,198023000,1,191155000,1,180916000,1,171116000,1,159139000,1,146364000,1,141214000,1,135535000,
		0,135007000,1,133852000,1,115170000,1,107350000,1,102532000,1,90355000,1,81196000,1,78078000,1,59129000,1,63026000,
		0,48130000,1,51305000),ncol=2,byrow = TRUE)
	hg38_loc = matrix(c(0,248957000,1,242194000,1,198296000,1,190215000,1,181539000,1,170806000,1,159346000,1,145139000,1,138395000,1,133798000,
		0,135087000,1,133276000,1,114365000,1,107044000,1,101992000,1,90339000,1,83258000,1,80374000,1,58617000,1,64445000,
		0,46710000,1,50818468),ncol=2,byrow = TRUE)
	print(paste0('Utilize hg',hg))
	if(hg==19){
		if(is.null(start)){
			start=hg19_loc[chr,1]
		}
		if(is.null(end)){
			end=hg19_loc[chr,2]
		}
	}else if(hg==38){
		if(is.null(start)){
			start=hg38_loc[chr,1]
		}
		if(is.null(end)){
			end=hg38_loc[chr,2]
		}
	}else{
		if(is.null(start) & is.null(end)){
			print('Need to specify start and end position')
			return(0)
		}
	}
	print(paste0('chr',chr,'; start ',start,';end ',end,';by ',by))
	s=seq(start-1,end,by)
	bed=matrix(NA,ncol=3,nrow=length(s)-1)
	bed[,1]=chr
	bed[,3]=s[-1]
	bed[,2]=s[-length(s)]+1
	options(scipen=10)
	write.table(bed,file=paste0('chr',chr,'.hg',hg,'.',start,'.',end,'.wgs.bed'),col.names = FALSE,row.names = FALSE,sep='\t',quote=FALSE)
	options(scipen=0)
}
