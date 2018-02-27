#' Plot CNV inference score.
#' 
#' Plot out CNV inference score. Each row is a sample, each column is a SNP or, 
#' exon (WES) or bin (WGS). Red color indicate score favor duplication whereas 
#' blue favor deletion.
#' 
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics axis grid legend par plot points
#' @importFrom stats aggregate dnorm dunif kmeans sd
#' @importFrom utils read.table write.table
#' @param icnv_res CNV inference result. Result from iCNV_detection() (i.e. iCNV_detection(...))
#' @param h start position of this plot. Default Start of the whole chromosome
#' @param t end position of this plot. Default End of the whole chromosome
#' @param title of this plot. Default "score plot"
#' @param output generated from output_list_function. If it isn't null, only CNVs in output file will be highlighted. Default NULL
#' @param col Specify if would like to plot in DGV color scheme (red for deletion, blue for duplication and grey for diploid) or default color scheme (blue for deletion, red for duplicatin and and green for diploid) Default NULL
#' @return void
#' @examples
#' plotHMMscore(icnv_res0,h=21000000, t=22000000, title='my favorite subject')
#' plotHMMscore(icnv_res0,h=21000000, t=22000000, title='my favorite subject',col='DGV')
#' @export
plotHMMscore <- function(icnv_res,h=NULL,t=NULL,title="score plot",output=NULL,col=''){
  if(is.null(h)){
    h <- min(icnv_res[[1]][[1]][[2]])
  }
  if(is.null(t)){
    t <- max(icnv_res[[1]][[1]][[2]])
  }
  HMMcalls <- icnv_res[[1]]
  CNV <- icnv_res[[2]]
  sel <- (h<=HMMcalls[[1]][[2]][,1] & HMMcalls[[1]][[2]][,2]<=t)
  Lpos <- HMMcalls[[1]][[2]][sel,]
  scores <- lapply(HMMcalls,function(x){x[[7]]})
  scores <- t(vapply(scores,function(x)x,scores[[1]]))[,sel]
  result <- lapply(HMMcalls,function(x){x[[1]]})
  result <- t(vapply(result,function(x)x,result[[1]]))[,sel]
  toplot <- scores
  l <- 1
  if(col=='DGV'){
    colfunc<-colorRampPalette(c("red","lightgrey","blue"))
    fields::image.plot(x=seq_len(ncol(toplot)),y=seq_len(nrow(toplot)),
      z=t(pmin(pmax(toplot,-l),l)),zlim=c(-l,l),col=colfunc(256),
      main = title,ylab='sample',xlab='')
  }else{
    fields::image.plot(x=seq_len(ncol(toplot)),y=seq_len(nrow(toplot)),
      z=t(pmin(pmax(toplot,-l),l)),zlim=c(-l,l),
      main = title,ylab='sample',xlab='')
  }
  if (!is.null(output)){
    addCNVtoplot2(output,Lpos,col)
  }
  else if (is.null(CNV)){
    addCNVtoplot1(result,col)
    if(col=='DGV'){
      legend("topright",c("del", "dup"),
             col = c('darkred','darkblue'),text.col = "green4", pch = c(20,20),cex = 0.75)
    }else{
      legend("topright",c("del", "dup"),
             col = c('white','black'),text.col = "green4", pch = c(1,20),cex = 0.75)
    }
  }
  else{
    CNV <- t(vapply(CNV,function(x)x,CNV[[1]]))[,sel]
    addCNVtoplot(CNV,col)
    if(col=='DGV'){
      legend("topright",c("0", "1", "3", "4"),
             col = c('darkred','magenta','cyan','blue'),text.col = "green4", 
             pch = c(20,20,20,20),cex = 0.75)
    }else{
      legend("topright",c("0", "1", "3", "4"),
             col = c('white','grey','magenta','black'),text.col = "green4", 
             pch = c(1,20,20,20),cex = 0.75)
    }
  }
}

addCNVtoplot2 <- function(output,Lpos,col){
  for (i in seq_along(output)){
    outputi <- output[[i]]
    if(length(outputi)>0){
      apply(outputi,1,function(x){
        if(x[1]==1){
          del <- which(Lpos[,1]>=x[2]&Lpos[,2]<=x[3])
          if(col=='DGV'){
            points(x=del,y=rep(i,length(del)),col='darkred',pch=16,cex=0.5)
          }else{
            points(x=del,y=rep(i,length(del)),col='white',pch=16,cex=0.5)
          }
        }
        if(x[1]==3){
          dup <- which(Lpos[,1]>=x[2]&Lpos[,2]<=x[3])
          if(col=='DGV'){
            points(x=dup,y=rep(i,length(dup)),col='darkblue',pch=16,cex=0.5)
          }else{
            points(x=dup,y=rep(i,length(dup)),col='black',pch=16,cex=0.5)
          }
        }
      })
      # cat(i,'del:',length(del),' dup:',length(dup),'\n')
    }
  }
}

addCNVtoplot1 <- function(result,col){
  for (i in seq_len(nrow(result))){
    del <- which(result[i,]==1)
    dup <- which(result[i,]==3)
    # cat(i,'del:',length(del),' dup:',length(dup),'\n')
    sel <- del
    if(col=='DGV'){
      points(x=sel,y=rep(i,length(sel)),col='darkred',pch=16,cex=0.5)
    }else{
      points(x=sel,y=rep(i,length(sel)),col='white',pch=16,cex=0.5)
    }
    sel <- dup
    if(col=='DGV'){
      points(x=sel,y=rep(i,length(sel)),col='darkblue',pch=16,cex=0.5)
    }else{
      points(x=sel,y=rep(i,length(sel)),col='black',pch=16,cex=0.5)
    }
  }
}

addCNVtoplot <- function(result,col){
  for (i in seq_len(nrow(result))){
    hemidel <- which(result[i,]==1)
    homodel <- which(result[i,]==0)
    cp1dup <- which(result[i,]==3)
    cp2dup <- which(result[i,]==4)
    # cat(i,'del:',length(del),' dup:',length(dup),'\n')
    if(col=='DGV'){
      sel <- hemidel
      points(x=sel,y=rep(i,length(sel)),col='magenta',pch=16,cex=0.5)
      sel <- homodel
      points(x=sel,y=rep(i,length(sel)),col='darkred',pch=16,cex=0.5)
      sel <- cp1dup
      points(x=sel,y=rep(i,length(sel)),col='cyan',pch=16,cex=0.5)
      sel <- cp2dup
      points(x=sel,y=rep(i,length(sel)),col='darkblue',pch=16,cex=0.5)
    }
    else{
      sel <- hemidel
      points(x=sel,y=rep(i,length(sel)),col='grey',pch=16,cex=0.5)
      sel <- homodel
      points(x=sel,y=rep(i,length(sel)),col='white',pch=16,cex=0.5)
      sel <- cp1dup
      points(x=sel,y=rep(i,length(sel)),col='magenta',pch=16,cex=0.5)
      sel <- cp2dup
      points(x=sel,y=rep(i,length(sel)),col='black',pch=16,cex=0.5)
    }
  }
}
