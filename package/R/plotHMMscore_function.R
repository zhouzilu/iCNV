#' Plot CNV inference score.
#' 
#' Plot out CNV inference score. Each row is a sample, each column is a SNP or, exon (WES) or bin (WGS). Red color indicate score
#' favor duplication whereas blue favor deletion.
#' 
#' @param icnv_res CNV inference result. Result from iCNV_detection() (i.e. iCNV_detection(...))
#' @param h start position of this plot. Default Start of the whole chromosome
#' @param t end position of this plot. Default End of the whole chromosome
#' @param title of this plot. Default "score plot"
#' @param output generated from output_list_function. If it isn't null, only CNVs in output file will be highlighted. Default NULL
#' @return void
#' @examples
#' pdf(file=paste0(projname,'.pdf'),width=13,height = 10)
#' plotHMMscore(icnv_res,h=20000000, t=30000000, title='my favorite subject')
#' dev.off()
#' @export
plotHMMscore=function(icnv_res,h=NULL,t=NULL,title="score plot",output=NULL){
  if(is.null(h)){
    h=min(icnv_res[[1]][[1]][[2]])
  }
  if(is.null(t)){
    t=max(icnv_res[[1]][[1]][[2]])
  }
  HMMcalls=icnv_res[[1]]
  CNV=icnv_res[[2]]
  sel=(h<=HMMcalls[[1]][[2]][,1] & HMMcalls[[1]][[2]][,2]<=t)
  Lpos=HMMcalls[[1]][[2]][sel,]
  scores=lapply(HMMcalls,function(x){x[[7]]})
  scores=t(sapply(scores,function(x)x, simplify = TRUE))[,sel]
  result=lapply(HMMcalls,function(x){x[[1]]})
  result=t(sapply(result,function(x)x, simplify = TRUE))[,sel]
  toplot=scores
  l=1
  fields::image.plot(x=seq(1,ncol(toplot)),y=seq(1,nrow(toplot)),z=t(pmin(pmax(toplot,-l),l)),zlim=c(-l,l),main = title,ylab='sample',xlab='')
  if (!is.null(output)){
    addCNVtoplot2(output,Lpos)
  }
  else if (is.null(CNV)){
    addCNVtoplot1(result)
    legend("topright",c("del", "dup"),
    col = c('white','black'),text.col = "green4", pch = c(1,20),cex = 0.75)
  }
  else{
    CNV=t(sapply(CNV,function(x)x, simplify = TRUE))[,sel]
    addCNVtoplot(CNV)
    legend("topright",c("0", "1", "3", "4"),
    col = c('white','grey','magenta','black'),text.col = "green4", pch = c(1,20,20,20),cex = 0.75)
  }
}

addCNVtoplot2=function(output,Lpos){
  for (i in 1:length(output)){
    outputi=output[[i]]
    if(length(outputi)>0){
      apply(outputi,1,function(x){
        if(x[1]==1){
          del=which(Lpos[,1]>=x[2]&Lpos[,2]<=x[3])
          points(x=del,y=rep(i,length(del)),col='white',pch=16,cex=0.5)
        }
        if(x[1]==3){
          dup=which(Lpos[,1]>=x[2]&Lpos[,2]<=x[3])
          points(x=dup,y=rep(i,length(dup)),col='black',pch=16,cex=0.5)
        }
      })
      # cat(i,'del:',length(del),' dup:',length(dup),'\n')
    }
  }
}

addCNVtoplot1=function(result){
  for (i in 1:nrow(result)){
    del = which(result[i,]==1)
    dup = which(result[i,]==3)
    # cat(i,'del:',length(del),' dup:',length(dup),'\n')
    sel=del
    points(x=sel,y=rep(i,length(sel)),col='white',pch=16,cex=0.5)
    sel=dup
    points(x=sel,y=rep(i,length(sel)),col='black',pch=16,cex=0.5)
  }
}

addCNVtoplot=function(result){
  for (i in 1:nrow(result)){
    hemidel = which(result[i,]==1)
    homodel = which(result[i,]==0)
    cp1dup = which(result[i,]==3)
    cp2dup = which(result[i,]==4)
    # cat(i,'del:',length(del),' dup:',length(dup),'\n')
    sel=hemidel
    points(x=sel,y=rep(i,length(sel)),col='grey',pch=16,cex=0.5)
    sel=homodel
    points(x=sel,y=rep(i,length(sel)),col='white',pch=16,cex=0.5)
    sel=cp1dup
    points(x=sel,y=rep(i,length(sel)),col='magenta',pch=16,cex=0.5)
    sel=cp2dup
    points(x=sel,y=rep(i,length(sel)),col='black',pch=16,cex=0.5)
  }
}
