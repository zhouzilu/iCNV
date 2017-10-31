#' Plot CNV inference score.
#' @param icnv_res CNV inference result. Result from iCNV_detection() (i.e. iCNV_detection(...))
#' @param h start position of this plot. Default Start of the whole chromosome
#' @param t end position of this plot. Default End of the whole chromosome
#' @param title of this plot. Default "score plot"
#' @return void
#' @examples
#' icnv_res = iCNV_detection(...)
#' pdf(file=paste0(projname,'.pdf'),width=13,height = 10)
#' plotHMMscore(icnv_res,h=100000, t=200000, subj='my favorite subject')
#' dev.off()
#' @export
plotHMMscore=function(icnv_res,h=min(icnv_res[[1]][[1]][[2]]),t=max(icnv_res[[1]][[1]][[2]]),subj="score plot"){
  HMMcalls=icnv_res[[1]]
  CNV=icnv_res[[2]]
  sel=(h<=HMMcalls[[1]][[2]][,1] & HMMcalls[[1]][[2]][,2]<=t)
  scores=lapply(HMMcalls,function(x){x[[7]]})
  scores=t(sapply(scores,function(x)x, simplify = T))[,sel]
  result=lapply(HMMcalls,function(x){x[[1]]})
  result=t(sapply(result,function(x)x, simplify = T))[,sel]
  toplot=scores
  l=1
  image.plot(x=seq(1,ncol(toplot)),y=seq(1,nrow(toplot)),z=t(pmin(pmax(toplot,-l),l)),zlim=c(-l,l),main = subj,ylab='sample',xlab='')
  if (is.null(CNV)){
    addCNVtoplot1(result)
    legend("topright",c("del", "dup"),
    col = c('white','black'),text.col = "green4", pch = c(1,20),cex = 0.75)
  }
  else{
    addCNVtoplot(CNV)
    legend("topright",c("0", "1", "3", "4"),
    col = c('white','grey','magenta','black'),text.col = "green4", pch = c(1,20,20,20),cex = 0.75)
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
  for (i in 1:length(result)){
    hemidel = which(result[[i]]==1)
    homodel = which(result[[i]]==0)
    cp1dup = which(result[[i]]==3)
    cp2dup = which(result[[i]]==4)
    del = which(result[[i]]<2)
    dup = which(result[[i]]>2)
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