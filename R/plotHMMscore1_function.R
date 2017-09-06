#' Plot CNV inference score.
#' @param testres CNV inference result. Output from iCNV_detection()
#' @param h start position of this plot. Default Start of the whole chromosome
#' @param t end position of this plot. Default End of the whole chromosome
#' @param title of this plot. Default "score plot"
#' @return void
#' @examples
#' pdf(file=paste0(projname,'.pdf'),width=13,height = 10)
#' plotHMMscore1(testres,h=100000, t=200000, subj='my favorite subject')
#' dev.off()
#' @export
plotHMMscore1=function(testres,h=min(testres[[1]][[2]]),t=max(testres[[1]][[2]]),subj="score plot"){
  sel=(h<=testres[[1]][[2]][,1] & testres[[1]][[2]][,2]<=t)
  scores=lapply(testres,function(x){x[[7]]})
  scores=t(sapply(scores,function(x)x, simplify = T))[,sel]
  result=lapply(testres,function(x){x[[1]]})
  result=t(sapply(result,function(x)x, simplify = T))[,sel]
  toplot=scores
  l=1
  image.plot(x=seq(1,ncol(toplot)),y=seq(1,nrow(toplot)),z=t(pmin(pmax(toplot,-l),l)),zlim=c(-l,l),main = subj,ylab='sample',xlab='')
  addCNVtoplot1(result)
  legend("topright",c("del", "dup"),
    col = c('white','black'),text.col = "green4", pch = c(1,20),cex = 0.75)
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

