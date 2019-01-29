#' Generate ouput list.
#' 
#' Generate human readable output from result calculated by iCNV_detection function
#' 
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics axis grid legend par plot points
#' @importFrom stats aggregate dnorm dunif kmeans sd
#' @importFrom utils read.table write.table
#' @param icnv_res CNV inference result. Output from iCNV_detection()
#' @param sampleid the name of the sample, same order as the input
#' @param CN An indicator variable with value {0,1} for whether exact copy number inferred in iCNV_detection. 0 no exact CN, 1 exact CN. Type integer. Default 0.
#' @param min_size A integer which indicate the minimum length of the CNV you are interested in. This could remove super short CNVs due to noise. Type integer. Default 0. Recommend 1000.
#' @param min_gap A integer which indicate the minimum gap of CNV you are likely to merge. The can remove gaps due to unstable estimation. However, the gap size have to be smaller than min_size to prevent error. Type integer. Default 0.
#' @return output CNV list of each individual
#' @examples
#' icnv.output <- output_list(icnv_res=icnv_res0,sampleid=sampname_qc, CN=0)
#' @export
output_list <- function(icnv_res,sampleid=NULL,CN=0,min_size=0,min_gap=0){
  stopifnot(is.numeric(CN))
  stopifnot(is.numeric(min_size))
  if (CN!=0){
    testres <- icnv_res[[1]]
    result <- lapply(icnv_res[[2]],function(x){x})
    Lpos <- lapply(testres,function(x){x[[2]]})
  }
  else{
    testres <- icnv_res[[1]]
    result <- lapply(testres,function(x){x[[1]]})
    Lpos <- lapply(testres,function(x){x[[2]]})
  }
  if(is.null(sampleid)){
      sampleid <- seq(1,length(testres))
  }
  # extract CNV region
  res <- mapply(function(I,pos){
    It <- I
    post <- pos
    i <- 2
    while(i<=length(It)){
      if (It[i]!=2 & It[i]==It[i-1]){
        It <- It[-i]
        temp <- post[i,2]
        post[i-1,2] <- temp
        post <- post[-i,]
      }
      else{
        i <- i+1
      }
    }
    return(list(It,post))
  },result,Lpos,SIMPLIFY = FALSE)
  
  # Generate calls in easy format
  icnv_res <- mapply(function(res){
    cnv <- res[[1]]
    pos <- res[[2]]
    ind <-  cnv!=2
    cnv <- cnv[ind]
    pos <- matrix(pos[ind,],ncol=2)
    return(cbind(cnv,pos))
  },res,SIMPLIFY = TRUE)
  names(icnv_res) <- sampleid
  
  # Merge min_gap
  icnv_res <- mapply(function(res){
    if(nrow(res)>1){
      not.merged=TRUE
      while(not.merged){
        not.merged=FALSE
        sel=which((res[-1,2]-res[-nrow(res),3])<min_gap)
        if (length(sel)>0){
          for(i in sel){
            if(res[i,1]==res[i+1,1]){
              not.merged=TRUE
              res[i,3]=res[i+1,3]
              res=res[-i,]
              break
            }
          }
        }
      }
      return(res)
    }else{
      return(res)
    }
  },icnv_res)
  # filt min_size
  icnv_res <- mapply(function(res){
    if(nrow(res)>0){
      filt <- ((res[,3]-res[,2])>=min_size)
      res=res[filt,,drop=FALSE]
      return(res)
    }else{
      return(res)
    }
  },icnv_res)
  
  return(icnv_res)
}