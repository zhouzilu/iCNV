#' Get BAM baf information from vcf
#' 
#' If your vcf follow the format in the example, you could use this function to extract NGS baf from vcf files. Remember to load library before hands.
#' Save 6 lists, each list has N entry. N = # of individuals (or vcf file)
#' ngs_baf.nm: name of the bamfiles; ngs_baf.chr: the chromosome; ngs_baf.pos: the position of the variants; 
#' ngs_baf: the BAF of the variants; ngs_baf.id: the ID of the variants; filenm:the file name
#' 
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics axis grid legend par plot points
#' @importFrom stats aggregate dnorm dunif kmeans sd
#' @importFrom utils read.table write.table
#' @param dir The directory to all the vcf stored; default is right in this folder. Defualt '.'
#' @param vcf_list All the vcf names stored in vcf.list; could use command:"ls *.vcf > vcf.list" to generate.
#' @param chr Specify the chromosome you want to generate. Must be of int from 1-22. If not specify, this function will generate all chromosomes. Defualt NULL
#' @param projectname Name of the project. Default ''
#' @return void 
#' @examples
#' dir <- system.file("extdata", package="iCNV")
#' bambaf_from_vcf(dir,'bam_vcf.list',projectname='icnv.demo.')
#' bambaf_from_vcf(dir,'bam_vcf.list',chr=22,projectname='icnv.demo.')
#' @export
bambaf_from_vcf = function(dir='.',vcf_list,chr=NULL,projectname=''){
    `%>%`=tidyr::`%>%`
    filenm=file.path(dir,read.table(file.path(dir,vcf_list),header=FALSE,as.is=TRUE)[[1]])

    # Read all the VCF data
    baf.all=lapply(filenm,bam.baf)

    # Seperate the data into different chromosome
    if(is.null(chr)){
        for (chr in seq_len(22)){
            baf.all.chr=lapply(baf.all,function(x){x1=x %>% dplyr::filter(`#CHROM`==chr); return(x1)})
            ngs_baf.nm=list()
            ngs_baf.chr=list()
            ngs_baf.pos=list()
            ngs_baf.id=list()
            ngs_baf=list()
            for (i in seq_along(baf.all.chr)){
                x=baf.all.chr[[i]]
                ngs_baf.nm[[i]]=names(x)[4]
                ngs_baf.chr[[i]]=x$`#CHROM`
                ngs_baf.pos[[i]]=x$POS
                ngs_baf.id[[i]]=x$ID
                ngs_baf[[i]]=x[[4]]
            }
            save(chr,ngs_baf.nm,ngs_baf.chr,ngs_baf.pos,ngs_baf,ngs_baf.id,filenm,file=paste0(projectname,'bambaf_',chr,'.rda'))
        }
    }else{
        baf.all.chr=lapply(baf.all,function(x){return(x %>% dplyr::filter(`#CHROM`==chr))})
        ngs_baf.nm=list()
        ngs_baf.chr=list()
        ngs_baf.pos=list()
        ngs_baf.id=list()
        ngs_baf=list()
        for (i in seq_along(baf.all.chr)){
            x=baf.all.chr[[i]]
            ngs_baf.nm[[i]]=names(x)[4]
            ngs_baf.chr[[i]]=x$`#CHROM`
            ngs_baf.pos[[i]]=x$POS
            ngs_baf.id[[i]]=x$ID
            ngs_baf[[i]]=x[[4]]
        }
        save(ngs_baf.nm,ngs_baf.chr,ngs_baf.pos,ngs_baf,ngs_baf.id,filenm,file=paste0(projectname,'bambaf_',chr,'.rda'))
    }
}

bam.baf=function(filenm){
    `%>%`=tidyr::`%>%`
    cat('load ',filenm,'...\n')
    dt = data.table::fread(filenm,skip='#CHROM')
    nm=names(dt)[length(dt)]
    cat('individual ',nm,'\n')
    dt <- dt %>% dplyr::filter(QUAL!='.') %>% dplyr::select(-REF,-ALT,-QUAL,-FILTER,-INFO)
    ind_AD=match('AD',strsplit(dt[[4]],':')[[1]])
    ind_DP=match('DP',strsplit(dt[[4]],':')[[1]])
    bafi <- vapply(strsplit(dt[[5]],':'),function(x){a=as.numeric(strsplit(x[ind_AD],',')[[1]]);return(a[1]/(a[1]+a[2]))},0.0)
    dt[nm] <- bafi
    dt <- dt %>% dplyr::select(-FORMAT)
    return(dt)
}
