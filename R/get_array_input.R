#' Get array information from given format
#' 
#' If your array input file follow the format in the example, you could use this
#'  function to extract array LRR and baf. Remember to load library before hands.
#' Recent modification also will automatically sort the location of the SNPs.
#' Save 4*[# of chr] lists, each list has N entry. N = # of individuals
#' snp_lrr: SNP LRR intensity; snp_lrr.pos: the position of the SNPs
#' snp_baf: the BAF of the SNPs; snp_baf.pos: the position of the SNPs
#' 
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics axis grid legend par plot points
#' @importFrom stats aggregate dnorm dunif kmeans sd
#' @importFrom utils read.table write.table
#' @param dir A string. The directory path to the folder where store signal intensity file according to chr. Type character 
#' @param pattern A string. The pattern of all the intensity file. Type character 
#' @param chr Specify the chromosome you want to generate. Must be of int from 1-22. If not specify, this function will generate files for all chromosomes. Default NULL
#' @param projname Name of the project. Type character 
#' @return void 
#' @examples
#' dir <- system.file("extdata", package="iCNV")
#' pattern <- paste0('*.csv.arrayicnv$')
#' get_array_input(dir,pattern,chr=22,projname='icnv.demo.')
#' @export

get_array_input <- function(dir=character(),pattern=character(),chr=NULL,projname=''){
    stopifnot(is.character(dir))
    stopifnot(is.character(pattern))
    stopifnot(is.character(projname))
    if(is.null(chr)){
        for (chr in seq_len(22)){
            cat(chr,'\n')
            dirPath <- dir
            signalFile <-  list.files(dirPath, pattern = pattern)
            signaldir <- file.path(dirPath, signalFile)
            signal.DT <- lapply(signaldir, data.table::fread, sep=",")
            snp_lrr <- lapply(signal.DT,function(t){t[[4]][t$Chr==chr]})
            names(snp_lrr) <- signalFile
            snp_baf <- lapply(signal.DT,function(t){t[[5]][t$Chr==chr]})
            names(snp_baf) <- signalFile
            snp_lrr.pos <- lapply(seq_along(snp_lrr), 
                function(i) signal.DT[[1]]$POS[signal.DT[[1]]$Chr==chr])
            snp_baf.pos <- lapply(seq_along(snp_baf), 
                function(i) signal.DT[[1]]$POS[signal.DT[[1]]$Chr==chr])

            # Sort snp locations
            snp_lrr=mapply(function(snp_lrri,snp_lrr.posi){
                ord=order(snp_lrr.posi)
                return(snp_lrri[ord])
            },snp_lrr,snp_lrr.pos,SIMPLIFY=FALSE)
            names(snp_lrr) <- signalFile

            snp_baf=mapply(function(snp_bafi,snp_baf.posi){
                ord=order(snp_baf.posi)
                return(snp_bafi[ord])
            },snp_baf,snp_baf.pos,SIMPLIFY=FALSE)
            names(snp_baf) <- signalFile
            
            snp_lrr.pos=mapply(function(snp_lrr.posi){
                ord=order(snp_lrr.posi)
                return(snp_lrr.posi[ord])
            },snp_lrr.pos,SIMPLIFY=FALSE)
            names(snp_lrr.pos) <- signalFile
            
            snp_baf.pos=mapply(function(snp_baf.posi){
                ord=order(snp_baf.posi)
                return(snp_baf.posi[ord])
            },snp_baf.pos,SIMPLIFY=FALSE)
            names(snp_baf.pos) <- signalFile
            
            save(snp_lrr,snp_baf,snp_lrr.pos,snp_baf.pos, 
                file=paste0(projname,'array_lrrbaf_',chr,'.rda'))
        }
    }else{
        cat(chr,'\n')
        dirPath <- dir
        signalFile <- list.files(dirPath, pattern = pattern)
        signaldir <- file.path(dirPath, signalFile)
        signal.DT <- lapply(signaldir, data.table::fread, sep=",")
        snp_lrr <- lapply(signal.DT,function(t){t[[4]]})
        names(snp_lrr) <- signalFile
        snp_baf <- lapply(signal.DT,function(t){t[[5]]})
        names(snp_baf) <- signalFile
        snp_lrr.pos <- lapply(seq_along(snp_lrr), 
            function(i) signal.DT[[1]]$POS)
        snp_baf.pos <- lapply(seq_along(snp_baf), 
            function(i) signal.DT[[1]]$POS)

        # Sort snp locations
        snp_lrr=mapply(function(snp_lrri,snp_lrr.posi){
            ord=order(snp_lrr.posi)
            return(snp_lrri[ord])
        },snp_lrr,snp_lrr.pos,SIMPLIFY=FALSE)
        names(snp_lrr) <- signalFile

        snp_baf=mapply(function(snp_bafi,snp_baf.posi){
            ord=order(snp_baf.posi)
            return(snp_bafi[ord])
        },snp_baf,snp_baf.pos,SIMPLIFY=FALSE)
        names(snp_baf) <- signalFile
        
        snp_lrr.pos=mapply(function(snp_lrr.posi){
            ord=order(snp_lrr.posi)
            return(snp_lrr.posi[ord])
        },snp_lrr.pos,SIMPLIFY=FALSE)
        names(snp_lrr.pos) <- signalFile
        
        snp_baf.pos=mapply(function(snp_baf.posi){
            ord=order(snp_baf.posi)
            return(snp_baf.posi[ord])
        },snp_baf.pos,SIMPLIFY=FALSE)
        names(snp_baf.pos) <- signalFile
        
        save(snp_lrr,snp_baf,snp_lrr.pos,snp_baf.pos, 
            file=paste0(projname,'array_lrrbaf_',chr,'.rda'))
    }
}


