#' Individual sample plot
#' 
#' Plot relationship between platforms and features for each individual. Only work for muli-platform inference.
#' 
#' @param ngs_plr A list of NGS intensity data. Each entry is an individual. If no NGS data, no need to specify.
#' @param snp_lrr A list of SNP array intensity data. Each entry is an individual. If no SNP array data, no need to specify.
#' @param ngs_baf A list of NGS BAF data. Each entry is an individual. If no NGS data, no need to specify.
#' @param snp_baf A list of SNP array BAF data. Each entry is an individual. If no SNP array data, no need to specify.
#' @param ngs_plr.pos A list of NGS intensity postion data. Each entry is an individual with dimension= (#of bins or exons, 2(start and end position)). If no NGS data, no need to specify.
#' @param snp_lrr.pos A list of SNP array intensity postion data. Each entry is an individual with length=#of SNPs. If no SNP array data, no need to specify.
#' @param ngs_baf.pos A list of NGS BAF postion data. Each entry is an individual with length=#of BAFs. If no NGS data, no need to specify.
#' @param snp_baf.pos A list of SNP array BAF postion data. Each entry is an individual with length=#of BAFs. If no SNP array data, no need to specify.
#' @param icnvres CNV inference result. The output from iCNV_detection()
#' @param I Indicating the position of the individual to plot
#' @param h start position of this plot. Default Start of the whole chromosome
#' @param t end position of this plot. Default End of the whole chromosome
#' @return void
#' @examples
#' pdf(file=paste0(projname,'.pdf'),width=13,height = 10)
#' plotindi(ngs_plr,snp_lrr,ngs_baf,snp_baf,
#'  ngs_plr.pos,snp_lrr.pos,ngs_baf.pos,snp_baf.pos,
#'  icnv_res,I=1,h=20000000, t=30000000)
#' dev.off()
#' @export
plotindi = function(ngs_plr,snp_lrr,ngs_baf,snp_baf,ngs_plr.pos,snp_lrr.pos,ngs_baf.pos,snp_baf.pos,icnvres,I,h=NULL,t=NULL){
  r1L=ngs_plr;r2L=snp_lrr;baf1=ngs_baf;baf2=snp_baf;rpos1=ngs_plr.pos;rpos2=snp_lrr.pos;bpos1=ngs_baf.pos;bpos2=snp_baf.pos
  hmmcalls = icnvres[[1]]
  if (is.null(h)){
    h=min(hmmcalls[[1]][[2]])
  }
  if (is.null(t)){
    t=max(hmmcalls[[1]][[2]])
  }
  r1i=r1L[[I]]
  r2i=r2L[[I]]
  baf1i=baf1[[I]]
  baf2i=baf2[[I]]
  rpos1i=rpos1[[I]]
  rpos2i=rpos2[[I]]
  bpos1i=bpos1[[I]]
  bpos2i=bpos2[[I]]
  res=hmmcalls[[I]]
  result=res[[1]]
  Lposi=res[[2]]
  rt=res[[3]]
  mu=res[[4]]
  sigma=res[[5]]
  score=res[[7]]

  sel=which(h<=Lposi[,1] & Lposi[,2]<=t)
  result=result[sel]
  score=score[sel]
  Lposi=Lposi[sel,]
  sel=which(h<=rpos1i[,1]&rpos1i[,2]<=t)
  r1i=r1i[sel]
  rpos1i=rpos1i[sel,]
  sel=which(h<=rpos2i & rpos2i<=t)
  r2i=r2i[sel]
  rpos2i=rpos2i[sel]
  sel=which(h<=bpos1i & bpos1i<=t)
  baf1i=baf1i[sel]
  bpos1i=bpos1i[sel]
  sel=which(h<=bpos2i & bpos2i<=t)
  baf2i=baf2i[sel]
  bpos2i=bpos2i[sel]
  b=1000
  ttlpos=seq(h-b/2,t+b,by=b)
  n=length(ttlpos)
  # row: z1, baf1, z2, baf2, score
  mat=matrix(NA,ncol=n-1,nrow=5)
  for (i in 1:(n-1)){
    sel=which(rpos2i>=ttlpos[i] & rpos2i<=ttlpos[i+1])
    if (length(sel)>0){
      mat[3,i]=mean(r2i[sel],na.rm=TRUE)
    }
    sel=which(bpos1i>=ttlpos[i] & bpos1i<=ttlpos[i+1])
    if (length(sel)>0){
      mat[2,i]=mean(baf1i[sel],na.rm=TRUE)
    }
    sel=which(bpos2i>=ttlpos[i] & bpos2i<=ttlpos[i+1])
    if (length(sel)>0){
      mat[4,i]=mean(baf2i[sel],na.rm=TRUE)
    }
  }
  for (i in 1:length(r1i)){
    sel=which(ttlpos>=rpos1i[i,1] & ttlpos<=rpos1i[i,2])
    if (length(sel)>0){
      if (all(sel==1)){
        mat[1,sel]=r1i[i]
        if (is.na(mat[2,sel])){mat[2,sel]=1}
      }
      else{
        mat[1,sel-1]=r1i[i]
        for(j in 1:length(sel)){
          if (is.na(mat[2,sel[j]-1])){
            mat[2,sel[j]-1]=1
          }
        }
      }
    }
    else{
      sel=max(which(ttlpos<=rpos1i[i,1]))
      # cat(sel)
      if (!is.infinite(sel)){
        mat[1,sel]=r1i[i]
        if (is.na(mat[2,sel])){
          mat[2,sel]=1
        }
      }
    }
  }
  for (i in 1:length(score)){
    a=Lposi[i,1]
    b=Lposi[i,2]
    if(a==b){
      sel=which(ttlpos<=a)
      mat[5,max(sel)]=score[i]
    }
    else {
      sel=which(ttlpos>=a & ttlpos<=b)
      if (length(sel>0)){
        mat[5,sel-1]=score[i]
      }
      else{
        sel=which(ttlpos<=a)
        mat[5,max(sel)]=score[i]
      }
    }
  }

  sel=which(!is.na(mat[5,]))
  mat=mat[,sel]
  ttlpos=ttlpos[sel]
  n2=ncol(mat)
  rmax=max(abs(mat[c(1,3),]),na.rm=TRUE)
  smax=max(abs(mat[5,]),na.rm=TRUE)
  r1max=max(abs(mat[1,]),na.rm=TRUE)
  r2max=max(abs(mat[3,]),na.rm=TRUE)
  x=1:n2
  par(mfrow=c(1,1))
  l=3
  fields::image.plot((as.matrix((pmin(pmax(mat[5,],-l),l)))),zlim=c(-l,l),axes=FALSE,main='score',ylab=I)
  del = which(result<2)
  dup = which(result>2)
  cat(I,' del:',length(del),' dup:',length(dup),'\n')
  sel=unique(unlist(mapply(function(x,y,pos){which(x<=pos & y>=pos)},Lposi[del,1],Lposi[del,2],MoreArgs = list(pos=ttlpos),SIMPLIFY = FALSE)))-1
  points(x=sel/length(ttlpos),y=rep(0,length(sel)),col='white',pch=20,cex=1)
  sel=unique(unlist(mapply(function(x,y,pos){which(x<=pos & y>=pos)},Lposi[dup,1],Lposi[dup,2],MoreArgs = list(pos=ttlpos),SIMPLIFY = FALSE)))-1
  points(x=sel/length(ttlpos),y=rep(0,length(sel)),col='black',pch=20,cex=1)
  sel=unique(unlist(mapply(function(x,y,pos){which(x==y & abs(x-pos)<1000)},Lposi[del,1],Lposi[del,2],MoreArgs = list(pos=ttlpos),SIMPLIFY = FALSE)))-1
  points(x=sel/length(ttlpos),y=rep(0,length(sel)),col='white',pch=20,cex=1)
  sel=unique(unlist(mapply(function(x,y,pos){which(x==y & abs(x-pos)<1000)},Lposi[dup,1],Lposi[dup,2],MoreArgs = list(pos=ttlpos),SIMPLIFY = FALSE)))-1
  points(x=sel/length(ttlpos),y=rep(0,length(sel)),col='black',pch=20,cex=1)
  legend("topright",c("del", "dup"),
    col = c('white','black'),bg='green2',text.col = c('white','black'), pch = c(20,20),cex = 0.75)

  par(xaxs='i')
  plot(x,y=mat[5,],pch=20,xlab="",type='l', axes=FALSE,col='green',ylab="",ylim=c(-smax,smax),xlim=c(0,n2),lwd=2)
  par(new=TRUE,xaxs='i')
  plot(x,y=mat[1,],pch=20,col='grey',cex=0.5,xlab="", ylab=I,ylim=c(-r1max,r1max),xlim=c(0,n2),main='Sequencing')
  par(new=TRUE,xaxs='i')
  plot(x,y=mat[2,],pch=20,cex=0.5,axes=FALSE,xlab="", ylab="",ylim=c(0,1),xlim=c(0,n2))
  axis(side = 4)
  legend("topright",c("intensity", "BAF","score"),
    col = c('grey','black','green'), pch = c(20,20,20),cex = 0.75)

  par(xaxs='i')
  plot(x,y=mat[5,],pch=20,xlab="",type='l', axes=FALSE,col='green',ylab="",ylim=c(-smax,smax),xlim=c(0,n2),lwd=2)
  par(new=TRUE,xaxs='i')
  plot(x,y=mat[3,],pch=20,col='grey',cex=0.5,xlab="",ylab=I,ylim=c(-r2max,r2max),xlim=c(0,n2),main='SNP')
  par(new=TRUE,xaxs='i')
  plot(x,y=mat[4,],pch=20,cex=0.5,axes=FALSE,xlab="", ylab="",ylim=c(0,1),xlim=c(0,n2))
  axis(side = 4)
  legend("topright",c("intensity", "BAF","score"),
    col = c('grey','black','green'), pch = c(20,20,20),cex = 0.75)
  par(mfrow=c(1,1))
}

