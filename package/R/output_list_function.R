#' Generate ouput list.
#' 
#' Generate human readable output from result calculated by iCNV_detection function
#' 
#' @param icnv_res CNV inference result. Output from iCNV_detection()
#' @param sampleid the name of the sample, same order as the input
#' @param CN An indicator variable with value {0,1} for whether exact copy number inferred in iCNV_detection. 0 no exact CN, 1 exact CN. Default 0.
#' @param min_size A integer which indicate the minimum length of the CNV you are interested in. This could remove super short CNVs due to noise. Default 0. Recommend 1000.
#' @return output CNV list of each individual
#' @examples
#' icnv.output = output_list(icnv_res=icnv_res,sampleid=sampname_qc, CN=0)
#' @export
output_list=function(icnv_res,sampleid=NULL,CN=0,min_size=0){
  if (CN!=0){
    testres=icnv_res[[1]]
    result=lapply(icnv_res[[2]],function(x){x})
    Lpos=lapply(testres,function(x){x[[2]]})
  }
  else{
    testres=icnv_res[[1]]
    result=lapply(testres,function(x){x[[1]]})
    Lpos=lapply(testres,function(x){x[[2]]})
  }
  if(is.null(sampleid)){
      sampleid = seq(1,length(testres))
  }
  # extract CNV region
  res=mapply(function(I,pos){
    It=I
    post=pos
    i=2
    while(i<=length(It)){
      if (It[i]!=2 & It[i]==It[i-1]){
        It=It[-i]
        temp=post[i,2]
        post[i-1,2]=temp
        post=post[-i,]
      }
      else{
        i=i+1
      }
    }
    return(list(It,post))
  },result,Lpos,SIMPLIFY = FALSE)

  icnv_res = mapply(function(res){
    cnv=res[[1]]
    pos=res[[2]]
    ind= cnv!=2
    cnv=cnv[ind]
    pos=matrix(pos[ind,],ncol=2)
    filt=(pos[,2]-pos[,1]>=min_size)
    cnv=cnv[filt]
    pos=matrix(pos[filt,],ncol=2)
    return(cbind(cnv,pos))
  },res,SIMPLIFY = TRUE)

  names(icnv_res)=sampleid
  return(icnv_res)
}