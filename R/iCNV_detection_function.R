#' CNV detection
#' 
#' Copy number variation detection tool for germline data. Able to combine intensity and BAF from SNP array and NGS data.
#' 
#' @param ngs_plr A list of NGS intensity data. Each entry is an individual. If no NGS data, no need to specify.
#' @param snp_lrr A list of SNP array intensity data. Each entry is an individual. If no SNP array data, no need to specify.
#' @param ngs_baf A list of NGS BAF data. Each entry is an individual. If no NGS data, no need to specify.
#' @param snp_baf A list of SNP array BAF data. Each entry is an individual. If no SNP array data, no need to specify.
#' @param ngs_plr.pos A list of NGS intensity postion data. Each entry is an individual with dimension= (#of bins or exons, 2(start and end position)). If no NGS data, no need to specify.
#' @param snp_lrr.pos A list of SNP array intensity postion data. Each entry is an individual with length=#of SNPs. If no SNP array data, no need to specify.
#' @param ngs_baf.pos A list of NGS BAF postion data. Each entry is an individual with length=#of BAFs. If no NGS data, no need to specify.
#' @param snp_baf.pos A list of SNP array BAF postion data. Each entry is an individual with length=#of BAFs. If no SNP array data, no need to specify.
#' @param maxIt An integer number indicate the maximum number of EM iteration if not converged during parameter inference. Default 50.
#' @param visual An indicator variable with value {0,1,2}. 0 indicates no visualization, 1 indicates basic visualization, 2 indicates complete visualization (Note visual 2 only work for single platform and integer CN inferenced)
#' @param projname A string as the name of this project. Default 'iCNV.'
#' @param CN An indicator variable with value {0,1} for whether wants to infer exact copy number. 0 no exact CN, 1 exact CN. Default 0.
#' @param mu A length tree vectur specify means of intensity in mixture normal distribution (Deletion, Diploid, Duplification). Default c(-3,0,2)
#' @param cap A boolean decides whether we cap insane intensity value due to double deletion or mutiple amplification. Default False
#' @keywords CNV, BAF, Platform integration, Intensity
#' @return (1) CNV inference, contains CNV inference, Start and end position for each inference, Conditional probability for each inference, mu for mixture normal, sigma for mixture normal, probability of CNVs, Z score for each inference.
#' @return (2) exact copy number for each CNV inference, if CN=1.
#' @examples
#' # icnv call without genotype (just infer deletion, duplication)
#' projname='icnv.demo.'
#' icnv_res0=iCNV_detection(ngs_plr,snp_lrr,
#'                          ngs_baf,snp_baf,
#'                          ngs_plr.pos,snp_lrr.pos,
#'                          ngs_baf.pos,snp_baf.pos,
#'                          projname=projname,CN=0,mu=c(-3,0,2),cap=TRUE,visual = 1)
#' # icnv call with genotype inference and complete plot
#' projname='icnv.demo.geno.'
#' icnv_res1=iCNV_detection(ngs_plr,snp_lrr,
#'                          ngs_baf,snp_baf,
#'                          ngs_plr.pos,snp_lrr.pos,
#'                          ngs_baf.pos,snp_baf.pos,
#'                          projname=projname,CN=1,mu=c(-3,0,2),cap=TRUE,visual = 2)
#' @export
iCNV_detection = function(ngs_plr=NULL,snp_lrr=NULL,ngs_baf=NULL,snp_baf=NULL,ngs_plr.pos=NULL,snp_lrr.pos=NULL,ngs_baf.pos=NULL,snp_baf.pos=NULL,maxIt=50,visual=0,projname='iCNV.',CN=0,mu=c(-3,0,2),cap=FALSE){
  # Change variable name for easier code writing
  r1L=ngs_plr;r2L=snp_lrr;baf1=ngs_baf;baf2=snp_baf;rpos1=ngs_plr.pos;rpos2=snp_lrr.pos;bpos1=ngs_baf.pos;bpos2=snp_baf.pos
  ptm <- proc.time()
  CNV=NULL
  c=checkdim(r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2)
  if(c[1]&c[2]){
    r2L=NULL;baf2=NULL;rpos2=NULL;bpos2=NULL
    n=length(r1L)
    indivd=seq(1,n)
    HMMcall = mapply(HMMEM,r1L,baf1,rpos1,bpos1,rep(maxIt,n),indivd,rep(list(mu),n),MoreArgs=list(r2i=r2L,baf2i=baf2,rpos2i=rpos2,bpos2i=bpos2,cap=cap),SIMPLIFY=FALSE)
    if (CN!=0){
      print('Start estimating integer copy number...')
      bafzIs=novisualization_ngs(HMMcall,r1L,baf1,rpos1,bpos1)
      CNV=exactCN_ngs(HMMcall,bafzIs[[1]],bafzIs[[2]],bafzIs[[3]],bafzIs[[4]],bafzIs[[5]],bafzIs[[6]],r1L,baf1,rpos1,bpos1)
      if (visual!=0){
        pdf(file=paste0(projname,'visual',visual,'.pdf'),width=13,height = 10)
        plotHMMscore(list(HMMcall,CNV),title=projname)
        dev.off()
      }
    }
    else{
      if (visual!=0){
        pdf(file=paste0(projname,'visual',visual,'.pdf'),width=13,height = 10)
        plotHMMscore(list(HMMcall,CNV),title=projname)
        dev.off()
      }
    }
    print(paste0('Inference time cost:',sum((proc.time() - ptm)[c(1,2)])))
  }
  else if(c[1]&c[3]){
    r1L=NULL;baf1=NULL;rpos1=NULL;bpos1=NULL
    n=length(r2L)
    indivd=seq(1,n)
    HMMcall = mapply(HMMEM,r2L,baf2,rpos2,bpos2,rep(maxIt,n),indivd,rep(list(mu),n),MoreArgs=list(r1i=r1L,baf1i=baf1,rpos1i=rpos1,bpos1i=bpos1,cap=cap),SIMPLIFY=FALSE)
    if (CN!=0){
      print('Start estimating integer copy number...')
      bafzIs=novisualization_snp(HMMcall,r2L,baf2,rpos2,bpos2)
      CNV=exactCN_snp(HMMcall,bafzIs[[1]],bafzIs[[2]],bafzIs[[3]],bafzIs[[4]],bafzIs[[5]],bafzIs[[6]],r2L,baf2,rpos2,bpos2)
      if (visual!=0){
        pdf(file=paste0(projname,'visual',visual,'.pdf'),width=13,height = 10)
        plotHMMscore(list(HMMcall,CNV),title=projname)
        dev.off()
      }
    }
    else{
      if (visual!=0){
        pdf(file=paste0(projname,'visual',visual,'.pdf'),width=13,height = 10)
        plotHMMscore(list(HMMcall,CNV),title=projname)
        dev.off()
      }
    }
    print(paste0('Inference time cost:',sum((proc.time() - ptm)[c(1,2)])))
  }
  else{
    n=length(r1L)
    indivd=seq(1,n)
    HMMcall = mapply(HMMEM,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,rep(maxIt,n),indivd,rep(list(mu),n),rep(list(cap),n),SIMPLIFY=FALSE)
    print(paste0('Inference time cost:',sum((proc.time() - ptm)[c(1,2)])))
    ptm <- proc.time()
    if (visual==2){
      pdf(file=paste0(projname,'visual',visual,'.pdf'),width=13,height = 10)
      if (CN!=0){
        print('Start estimating integer copy number...')
        bafzIs=visualization(HMMcall,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2)
        CNV=exactCN(HMMcall,bafzIs[[4]],bafzIs[[5]],bafzIs[[6]],bafzIs[[7]],bafzIs[[8]],bafzIs[[9]],r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2)
        plotHMMscore(list(HMMcall,CNV),title=projname)
        visualization2(HMMcall,CNV,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2)
      }
      dev.off()
    }
    else if (visual==1){
      pdf(file=paste0(projname,'visual',visual,'.pdf'),width=13,height = 10)
      if (CN!=0){
        print('Start estimating integer copy number...')
        bafzIs=novisualization(HMMcall,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2)
        CNV=exactCN(HMMcall,bafzIs[[1]],bafzIs[[2]],bafzIs[[3]],bafzIs[[4]],bafzIs[[5]],bafzIs[[6]],r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2)       
      }
      plotHMMscore(list(HMMcall,CNV),title=projname)
      dev.off()
    }
    else{
      if (CN!=0){
        print('Start estimating integer copy number...')
        bafzIs=novisualization(HMMcall,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2)
        CNV=exactCN(HMMcall,bafzIs[[1]],bafzIs[[2]],bafzIs[[3]],bafzIs[[4]],bafzIs[[5]],bafzIs[[6]],r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2)
      }
    }
    print(paste0('Visualization time cost:',sum((proc.time() - ptm)[c(1,2)])))
  }
  return(list(HMMcall,CNV))
}

checkdim = function(r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2){
  c1=!all(length(r1L) == c(length(r1L),length(r2L),length(baf1),length(baf2),length(rpos1),length(rpos2),length(bpos1),length(rpos2)))
  c2=all(length(r1L) == c(length(r1L),length(baf1),length(rpos1),length(bpos1)))&length(r1L)!=0
  c3=all(length(r2L) == c(length(r2L),length(baf2),length(rpos2),length(rpos2)))&length(r2L)!=0
  if(c1){
    if(c2 & !c3){
      cat('Dimension of dataset is wrong in array. Use sequencing data only. Exact CN inference down')
      c3=FALSE
    }
    else if(c3 & !c2){
      cat('Dimension of dataset is wrong in sequencing. Use array data only. Exact CN inference down')
      c2=FALSE
    }
    else{
      cat('Dimension of dataset is wrong. Please check.')
      break
    }
  }
  return(c(c1,c2,c3))
}

HMMEM = function(r1i,r2i,baf1i,baf2i,rpos1i,rpos2i,bpos1i,bpos2i,maxIt,ind,mu,cap){
  # Init
  cat('\n\nindividual ',ind,': \n')
  sigma=c(1,1,1)
  p=10^-8
  sigmaM=matrix(NA,ncol=3,nrow=maxIt)
  sigmaM[1,]=sigma
  emissionBLoc = NULL # Use for NGS emission probability
  em1Loc = NULL # Use for combine NGS and ARRAY emission probability
  em2Loc = NULL # Use for combine NGS and ARRAY emission probability
  cat('iteration',1,': p=',p,'; mu=',mu,'; sigma=',sigma,'; sum of difference',Inf,'\n')
  for (i in 2:maxIt){
    res=HMMiEM(r1i,r2i,baf1i,baf2i,rpos1i,rpos2i,bpos1i,bpos2i,pir,pib,mu,sigma,p,cap,emissionBLoc,em1Loc,em2Loc)
    result=res[[1]]
    Lposi=res[[2]]
    rt=res[[3]]
    sigmaM[i,]=sigma=res[[4]]
    score=res[[5]]
    emissionBLoc=res[[6]]
    em1Loc=res[[7]]
    em2Loc=res[[8]]
    dif=sum((sigma-sigmaM[i-1,])^2)
    cat('iteration',i,': p=',p,'; mu=',mu,'; sigma=',sigma,'; sum of difference',dif,'\n')
    if( dif< 0.001){
      break
    }
  }
  return(list(result,Lposi,rt,mu,sigma,p,score))
}

HMMiEM = function(r1i,r2i,baf1i,baf2i,rpos1i,rpos2i,bpos1i,bpos2i,pir,pib,mu,sigma,p,cap,emissionBLoc,em1Loc,em2Loc){
  # 1 for exome and 2 for array
  # logR function
  emissionR1=function(x){dnorm(x,mu[1],sigma[1],log = TRUE)}
  emissionR2=function(x){dnorm(x,mu[2],sigma[2],log = TRUE)}
  emissionR3=function(x){dnorm(x,mu[3],sigma[3],log = TRUE)}
  # BAF function
  # 0 or 1
  bsig1=0.05
  emissionB1=function(x){log(0.5*truncnorm::dtruncnorm(x,0,1,0,bsig1)+0.5*truncnorm::dtruncnorm(x,0,1,1,bsig1))}
  # 0 or 0.5 or 1
  bsig=0.1
  emissionB2=function(x){log(0.25*truncnorm::dtruncnorm(x,0,1,0,bsig)+0.25*truncnorm::dtruncnorm(x,0,1,1,bsig)+0.5*truncnorm::dtruncnorm(x,0,1,0.5,bsig))}
  # 0 or 0.25 or 0.33 or 0.5 or 0.67 or 0.75 or 1
  emissionB3=function(x){log(0.14*truncnorm::dtruncnorm(x,0,1,0,bsig)+0.14*truncnorm::dtruncnorm(x,0,1,1,bsig)+0.16*truncnorm::dtruncnorm(x,0,1,0.5,bsig)+0.14*truncnorm::dtruncnorm(x,0,1,0.25,bsig)+0.14*truncnorm::dtruncnorm(x,0,1,0.33,bsig)+0.14*truncnorm::dtruncnorm(x,0,1,0.67,bsig)+0.14*truncnorm::dtruncnorm(x,0,1,0.75,bsig))}

  # 1
  emission11=NULL
  emission12=NULL
  emission13=NULL
  if(!is.null(r1i)){
    #normalize log ratio
    z1i = (r1i - mean(r1i,na.rm=TRUE))/sd(r1i,na.rm=TRUE);baf1i[is.na(baf1i)]=1;z1i[is.na(z1i)]=0
    if (cap==TRUE){
      # Cap intensity
      z1i=pmax(pmin(z1i,30),-30)
    }
    # emission probability r1
    emissionR11=emissionR1(z1i)
    emissionR12=emissionR2(z1i)
    emissionR13=emissionR3(z1i)
    # emission probability baf1
    rpos1is = rpos1i[,1]
    rpos1ie = rpos1i[,2]

    calexomeLoc = function(x,y,posb){
      return(which(posb>=x & posb<=y))
    }
    if(is.null(emissionBLoc)){
      emissionBLoc = mapply(calexomeLoc,rpos1is,rpos1ie,MoreArgs = list(posb = bpos1i),SIMPLIFY = FALSE)
    }
    calexomeP = function(indb,baf,fx){
      b = baf[indb]
      if(length(b>0)){
        return(sum(c(fx(b)),na.rm=TRUE))
      }
      else{
        return(fx(0))
      }
    }
    emissionB11 = (mapply(calexomeP,emissionBLoc,MoreArgs = list(baf = baf1i, fx=emissionB1)))
    emissionB12 = (mapply(calexomeP,emissionBLoc,MoreArgs = list(baf = baf1i, fx=emissionB2)))
    emissionB13 = (mapply(calexomeP,emissionBLoc,MoreArgs = list(baf = baf1i, fx=emissionB3)))

    # combine emission from r and baf
    emission11=emissionR11+emissionB11
    emission12=emissionR12+emissionB12
    emission13=emissionR13+emissionB13
  }
  # 2
  emission21=NULL
  emission22=NULL
  emission23=NULL
  if(!is.null(r2i)){
    #normalize log ratio
    z2i = (r2i-mean(r2i,na.rm=TRUE))/sd(r2i,na.rm=TRUE);baf2i[is.na(baf2i)]=1;z2i[is.na(z2i)]=0
    if (cap==TRUE){
      # Cap intensity
      z2i=pmax(pmin(z2i,30),-30)
    }
    # emission probability r2
    emissionR21=emissionR1(z2i)
    emissionR22=emissionR2(z2i)
    emissionR23=emissionR3(z2i)
    emissionB21=(emissionB1(baf2i))
    emissionB22=(emissionB2(baf2i))
    emissionB23=(emissionB3(baf2i))
    # combine emission from r and baf
    emission21=emissionR21+emissionB21
    emission22=emissionR22+emissionB22
    emission23=emissionR23+emissionB23
  }

  if(!is.null(emission11)&!is.null(emission21)){
    # match position
    sel= unique(unlist(mapply(function(x,y,pos){which(x<=pos & y>=pos)},rpos1is,rpos1ie,MoreArgs = list(pos=rpos2i),SIMPLIFY = FALSE)))
    rpos2iUexon = rpos2i[-sel]
    rpos2iUexon = matrix(rep(rpos2iUexon,2),ncol=2,byrow=FALSE)
    Lposi = rbind(rpos2iUexon,rpos1i)
    Lposi = Lposi[order(Lposi[,1]),]
    pos_exom=floor((rpos1is+rpos1ie)/2)

    calem1Loc = function(x,y,pos1){
      return(which(pos1>=x & pos1<=y))
    }
    calem2Loc = function(x,y,pos2){
      return(which(pos2>=x & pos2<=y))
    }
    if(is.null(em1Loc)){
      em1Loc = mapply(calem1Loc,Lposi[,1],Lposi[,2],MoreArgs=list(pos1=pos_exom),SIMPLIFY = FALSE)
    }
    if(is.null(em2Loc)){
      em2Loc = mapply(calem2Loc,Lposi[,1],Lposi[,2],MoreArgs=list(pos2=rpos2i),SIMPLIFY = FALSE)
    }
    
    calem12 = function(Loc1,Loc2,em1,em2){
      r=sum(c(sum(em1[Loc1]),sum(em2[Loc2])),na.rm=TRUE)
      return(r)
    }
    # total emission probability
    emission1=mapply(calem12,em1Loc,em2Loc,MoreArgs = list(em1=emission11,em2=emission21))
    emission2=mapply(calem12,em1Loc,em2Loc,MoreArgs = list(em1=emission12,em2=emission22))
    emission3=mapply(calem12,em1Loc,em2Loc,MoreArgs = list(em1=emission13,em2=emission23))

  }else if(!is.null(emission11)&is.null(emission21)){
    Lposi = rpos1i
    emission1=emission11
    emission2=emission12
    emission3=emission13
  }else if (is.null(emission11)&!is.null(emission21)){
    Lposi = matrix(rep(rpos2i,2),ncol=2,byrow=FALSE)
    emission1=emission21
    emission2=emission22
    emission3=emission23
  }else{
    cat('Intensity and BAF emission probability are null.')
  }

  # cat('emission 1 range:', range(emission1),'emission 2 range:', range(emission2),'emission 3 range:', range(emission3),'\n')

  n = nrow(Lposi)
  D=70000
  d=pmax(rowMeans(Lposi)[2:n]-rowMeans(Lposi)[1:(n-1)],1)
  f=exp(-d/D)
  # p=10^-8
  q=1/6
  transitionarray = array(NA,dim=c(3,3,n-1))
  transitionarray[1,1,]=transition11=log(f*(1-q)+(1-f)*p)
  transitionarray[1,2,]=transition12=log(f*q+(1-f)*(1-2*p))
  transitionarray[1,3,]=transition13=log((1-f)*p)
  transitionarray[2,1,]=transition21=log(rep(p,length(d)))
  transitionarray[2,2,]=transition22=log(rep(1-2*p,length(d)))
  transitionarray[2,3,]=transition23=log(rep(p,length(d)))
  transitionarray[3,1,]=transition31=log((1-f)*p)
  transitionarray[3,2,]=transition32=log(f*q+(1-f)*(1-2*p))
  transitionarray[3,3,]=transition33=log(f*(1-q)+(1-f)*p)
  emissionmatrix = matrix(NA,nrow=3,ncol=n)
  emissionmatrix[1,]=emission1
  emissionmatrix[2,]=emission2
  emissionmatrix[3,]=emission3
  # viterbi algorithm
  v=matrix(NA,nrow=n,ncol=3)
  pointer=matrix(NA,nrow=(n+1),ncol=3)
  v[1,1]=emission1[1]+log(p)
  v[1,2]=emission2[1]+log(1-2*p)
  v[1,3]=emission3[1]+log(p)
  pointer[1,1]=0
  pointer[1,2]=0
  pointer[1,3]=0

  for (i in 2:n){
    # max1=max((v[i-1,1]+(transition11[i-1])),(v[i-1,2]+(transition21[i-1])),(v[i-1,3]+(transition31[i-1])))
    max1=max(v[i-1,]+transitionarray[,1,i-1])
    v[i,1]=emission1[i]+max1
    if ((v[i-1,1]+(transition11[i-1]))==max1){
      pointer[i,1]=1
    } else if ((v[i-1,2]+(transition21[i-1]))==max1){
      pointer[i,1]=2
    } else{
      pointer[i,1]=3
    }
    # max2=max((v[i-1,1]+(transition12[i-1])),(v[i-1,2]+(transition22[i-1])),(v[i-1,3]+(transition32[i-1])))
    max2=max(v[i-1,]+transitionarray[,2,i-1])
    v[i,2]=emission2[i]+max2
    if ((v[i-1,1]+(transition12[i-1]))==max2){
      pointer[i,2]=1
    } else if ((v[i-1,2]+(transition22[i-1]))==max2){
      pointer[i,2]=2
    } else {
      pointer[i,2]=3
    }
    # max3=max((v[i-1,1]+(transition13[i-1])),(v[i-1,2]+(transition23[i-1])),(v[i-1,3]+(transition33[i-1])))
    max3=max(v[i-1,]+transitionarray[,3,i-1])
    v[i,3]=emission3[i]+max3
    if ((v[i-1,1]+(transition13[i-1]))==max3){
      pointer[i,3]=1
    } else if ((v[i-1,2]+(transition23[i-1]))==max3){
      pointer[i,3]=2
    } else {
      pointer[i,3]=3
    }
  }

  # termination
  max4=max(v[n,1],v[n,2],v[n,3])
  if (v[n,1]==max4){
    pointer[n+1,1]=1
  } else if (v[n,2]==max4) {
    pointer[n+1,1]=2
  } else {
    pointer[n+1,1]=3
  }
  # traceback
  traceback=pointer[n+1,1]
  result=rep(NA,n)
  for (i in seq(n+1,2,-1)){
    if (traceback==1){
      result[i-1]=1
      traceback=pointer[i-1,1]
    } else if (traceback==2) {
      result[i-1]=2
      traceback=pointer[i-1,2]
    } else {
      result[i-1]=3
      traceback=pointer[i-1,3]
    }
  }

  # index matrix
  I=matrix(0,ncol=3,nrow=n)
  I[which(result==1),1]=1
  I[which(result==2),2]=1
  I[which(result==3),3]=1

  # posterior probability
  # fast pmin and pmax function
  pmin. <- function(k,x) (x+k - abs(x-k))/2 
  pmax. <- function(k,x) (x+k + abs(x-k))/2 
  # Or E step
  # Forward Algorithm
  at=matrix(NA,nrow=n,ncol=3)
  cta=rep(NA,n)
  at[1,1]=(emission1[1])+log(p)
  at[1,2]=(emission2[1])+log(1-2*p)
  at[1,3]=(emission3[1])+log(p)
  for (i in 2:n){
    wt=sum(range(at[i-1,]))/2
    #s1=wt+log(sum(exp(pmin(pmax(c(at[i-1,1]+transition11[i-1]-wt,at[i-1,2]+transition21[i-1]-wt,at[i-1,3]+transition31[i-1]-wt),-745),709))))
    s1=wt+log(sum(exp(pmin.(pmax.(at[i-1,]+transitionarray[,1,i-1]-wt,-745),709))))
    
    at[i,1]=emission1[i]+s1
    #s2=wt+log(sum(exp(pmin(pmax(c(at[i-1,1]+transition12[i-1]-wt,at[i-1,2]+transition22[i-1]-wt,at[i-1,3]+transition32[i-1]-wt),-745),709))))
    s2=wt+log(sum(exp(pmin.(pmax.(at[i-1,]+transitionarray[,2,i-1]-wt,-745),709))))
    at[i,2]=emission2[i]+s2
    #s3=wt+log(sum(exp(pmin(pmax(c(at[i-1,1]+transition13[i-1]-wt,at[i-1,2]+transition23[i-1]-wt,at[i-1,3]+transition33[i-1]-wt),-745),709))))
    s3=wt+log(sum(exp(pmin.(pmax.(at[i-1,]+transitionarray[,3,i-1]-wt,-745),709))))
    at[i,3]=emission3[i]+s3
    wt=sum(range(at[i,]))/2
    cta[i]=-(wt+log(sum(exp(pmin.(pmax.(at[i,]-wt,-745),709)))))
    at[i,]=at[i,]+cta[i]
  }

  # Backward Algorithm
  bt=matrix(NA,nrow=n,ncol=3)
  ctb=rep(NA,n)
  bt[n,1]=0
  bt[n,2]=0
  bt[n,3]=0
  for (i in seq(n-1,1,-1)){
    wt=mean(range(c(emission1[i+1],emission2[i+1],emission3[i+1])))
    s1=wt+log(sum(exp(pmin.(pmax.(transitionarray[1,,i]+emissionmatrix[,i+1]+bt[i+1,]-wt,-745),709))))
    s2=wt+log(sum(exp(pmin.(pmax.(transitionarray[2,,i]+emissionmatrix[,i+1]+bt[i+1,]-wt,-745),709))))
    s3=wt+log(sum(exp(pmin.(pmax.(transitionarray[3,,i]+emissionmatrix[,i+1]+bt[i+1,]-wt,-745),709))))
    wt=mean(range(c(s1,s2,s3)))
    ctb[i]=-(wt+log(sum(exp(pmin.(pmax.(c(s1-wt,s2-wt,s3-wt),-745),709)))))
    bt[i,]=ctb[i]+(c(s1,s2,s3))
  }
  abt=pmin.(pmax.(at+bt,-745),709)
  rt=(exp(abt)+1e-323)/(rowSums(exp(abt))+1e-323)
  rt=rt/rowSums(rt)
  ct=cta+ctb

  #calculate score
  calscore=function(x){
    sel=which.max(x[c(1,3)])
    p=sqrt(-log(x[2]))
    if (length(sel)==0){
      cat('Probability out of bound! Try cap your intensity data.')
    }
    if (sel==1){return (-1*p)}
    else {return (p)}
  }
  score=apply(rt,1,calscore)

  # M step
  getZs = function(Loc1,Loc2,z1,z2){
    return(mean(c((z1[Loc1]),(z2[Loc2])),na.rm=TRUE))
  }
  if(!is.null(emission11)&!is.null(emission21)){
    zs=mapply(getZs, em1Loc,em2Loc, MoreArgs = list(z1=z1i,z2=z2i))
  }else if(!is.null(emission11)&is.null(emission21)){
    zs=z1i
  }else if(is.null(emission11)&!is.null(emission21)){
    zs=z2i
  }else{cat('Intensity and BAF emission probability are NULL')}

  sigmahat= rep(sqrt(sum(rt*(cbind(zs-mu[1],zs-mu[2],zs-mu[3])^2))/sum(rt)),3)
  return (list(result,Lposi,rt,sigmahat,score,emissionBLoc,em1Loc,em2Loc))
}

visualization = function(testres,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2){
  cat('\n','Generating visualization plot. This may take a while...','\n')
  result=lapply(testres,function(x){x[[1]]})
  Lpos=lapply(testres,function(x){x[[2]]})
  ttldipposlist=mapply(function(x){which(x==2)},result,SIMPLIFY=FALSE)
  ttldelposlist=mapply(function(x){which(x==1)},result,SIMPLIFY=FALSE)
  ttldupposlist=mapply(function(x){which(x==3)},result,SIMPLIFY=FALSE)
  nzvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldipposlist,Lpos,SIMPLIFY = FALSE)
  nbafs=unlist(mapply(function(x){x[,2]},nzvsbaf))
  nzs=unlist(mapply(function(x){x[,1]},nzvsbaf))
  nI=unlist(mapply(function(x){x[,3]},nzvsbaf))
  dzvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldelposlist,Lpos,SIMPLIFY = FALSE)
  dbafs=unlist(mapply(function(x){x[,2]},dzvsbaf))
  dbafs[is.na(dbafs)]=0
  dzs=unlist(mapply(function(x){x[,1]},dzvsbaf))
  dI=unlist(mapply(function(x){x[,3]},dzvsbaf))
  pzvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldupposlist,Lpos,SIMPLIFY = FALSE)
  pbafs=unlist(mapply(function(x){x[,2]},pzvsbaf))
  pbafs[is.na(pbafs)]=0
  pzs=unlist(mapply(function(x){x[,1]},pzvsbaf))
  pI=unlist(mapply(function(x){x[,3]},pzvsbaf))
  x.lim=range(c(range(pzs,na.rm=TRUE),range(dzs,na.rm=TRUE),range(nzs,na.rm=TRUE)))
  plot(x=nzs,y=nbafs,pch=20,cex=0.5,xlim=x.lim,ylim=c(0,0.5),col='grey',xlab='log ratio zscore',ylab='mBAF')
  grid()
  points(x=dzs[dI==2],y=dbafs[dI==2],pch=20,cex=0.8,col='blue')
  points(x=dzs[dI==2.5],y=dbafs[dI==2.5],pch=20,cex=0.8,col='cyan')
  points(x=dzs[dI==1],y=dbafs[dI==1],pch=20,cex=0.8,col='cornflowerblue')
  points(x=dzs[dI==1.5],y=dbafs[dI==1.5],pch=20,cex=0.8,col='deepskyblue')
  points(x=pzs[pI==2],y=pbafs[pI==2],pch=20,cex=0.8,col='red')
  points(x=pzs[pI==2.5],y=pbafs[pI==2.5],pch=20,cex=0.8,col='pink')
  points(x=pzs[pI==1],y=pbafs[pI==1],pch=20,cex=0.8,col='magenta')
  points(x=pzs[pI==1.5],y=pbafs[pI==1.5],pch=20,cex=0.8,col='purple')
  legend(x.lim[1],0.5, c("Del SNPs only", "Del SNPs in Exon", "Del Exon w/ BAFs", "Del Exon w/ no BAF","Dup SNPs only","Dup SNPs in Exon","Dup Exon w/BAFs","Dup Exon w/ no BAF"),
    col = c('blue','cyan','cornflowerblue','deepskyblue','red','pink','magenta','purple'),text.col = "green4", pch = c(20,20,20,20),cex = 0.75)
  return(list(nbafs,nzs,nI,dbafs,dzs,dI,pbafs,pzs,pI,x.lim))
}

# Pair intensity and BAF
BAFvsZ=function(r1i,r2i,baf1i,baf2i,rpos1i,rpos2i,bpos1i,bpos2i,CNVs,Lposi){
  bafs=c()
  rs=c()
  I=c()
  if (length(CNVs)>0){
    z1i=(r1i-mean(r1i,na.rm=TRUE))/sd(r1i,na.rm=TRUE)
    z2i=(r2i-mean(r2i,na.rm=TRUE))/sd(r2i,na.rm=TRUE)
    mbaf1i=0.5-abs(baf1i-0.5)
    mbaf2i=0.5-abs(baf2i-0.5)
    for (j in 1:length(CNVs)){
      s=Lposi[CNVs[j],1]
      t=Lposi[CNVs[j],2]
      if (s==t){
        sel1=which(rpos2i==s)
        sel2=which(bpos2i==s)
        bafs=c(bafs,mbaf2i[sel2])
        rs=c(rs,z2i[sel1])
        I=c(I,2)
      }
      else{
        sel1=which(rpos2i>=s&rpos2i<=t)
        sel2=which(bpos2i>=s&bpos2i<=t)
        if (length(sel1)>0){
          bafs=c(bafs,mbaf2i[sel2])
          rs=c(rs,z2i[sel1])
          I=c(I,rep(2.5,length(sel1)))
        }
        sel1=which(rpos1i[,1]>=s&rpos1i[,2]<=t)
        sel2=which(bpos1i>=s&bpos1i<=t)
        if (length(sel2)>0){
          bafs=c(bafs,mbaf1i[sel2])
          rs=c(rs,rep(z1i[sel1],length(sel2)))
          I=c(I,rep(1,length(sel2)))
        }
        else{
          bafs=c(bafs,0)
          rs=c(rs,z1i[sel1])
          I=c(I,1.5)
        }
      }
    }
  }
  return(cbind(rs,bafs,I))
}

# Infer exact Copy Number
exactCN = function(testres,dbafs,dzs,dI,pbafs,pzs,pI,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2){
  result=lapply(testres,function(x){x[[1]]})
  Lpos=lapply(testres,function(x){x[[2]]})
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
  mresult=lapply(res,function(x){x[[1]]})
  mLpos=lapply(res,function(x){x[[2]]})
  # Distribution mean inference using K means
  dz1=dzs[dI==1|dI==1.5]
  if(length(dz1)<3){
    dzprob11=function(x){log(0.99)}
    dzprob12=function(x){log(0.01)}
  }
  else{
    dz1fit <- kmeans(dz1, 2)
    dmu1=aggregate(dz1,by=list(dz1fit$cluster),FUN=mean)$x
    dzprob11=function(x){dnorm(x,max(dmu1),1,log=TRUE)}
    dzprob12=function(x){dnorm(x,min(dmu1),1,log=TRUE)}
  }
  dz2=dzs[dI==2|dI==2.5]
  if(length(dz2)<3){
    dzprob21=function(x){log(0.99)}
    dzprob22=function(x){log(0.01)}
  }
  else{
    dz2fit <- kmeans(dz2, 2)
    dmu2=aggregate(dz2,by=list(dz2fit$cluster),FUN=mean)$x
    dzprob21=function(x){dnorm(x,max(dmu2),1,log=TRUE)}
    dzprob22=function(x){dnorm(x,min(dmu2),1,log=TRUE)}
  }
  dbprob1=function(x){log(0.5*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.5*truncnorm::dtruncnorm(x,0,1,1,0.1))}
  dbprob2=function(x){log(0.45*truncnorm::dtruncnorm(x,0,1,0,0.01)+0.45*truncnorm::dtruncnorm(x,0,1,1,0.01)+0.1*dunif(x,0,1,log=TRUE))}
  
  pz1=pzs[pI==1|pI==1.5]
  if(length(pz1)<3){
    pzprob11=function(x){log(0.99)}
    pzprob12=function(x){log(0.01)}
  }
  else{
    pz1fit <- kmeans(pz1, 2)
    pmu1=aggregate(pz1,by=list(pz1fit$cluster),FUN=mean)$x
    pzprob11=function(x){dnorm(x,min(pmu1),1,log=TRUE)}
    pzprob12=function(x){dnorm(x,max(pmu1),1,log=TRUE)}
  }
  pz2=pzs[pI==2|pI==2.5]
  if(length(pz2)<3){
    pzprob21=function(x){log(0.99)}
    pzprob22=function(x){log(0.01)}
  }
  else{
    pz2fit <- kmeans(pz2, 2)
    pmu2=aggregate(pz2,by=list(pz2fit$cluster),FUN=mean)$x
    pzprob21=function(x){dnorm(x,min(pmu2),1,log=TRUE)}
    pzprob22=function(x){dnorm(x,max(pmu2),1,log=TRUE)}
  }
  pbprob1=function(x){log(0.25*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,1,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,0.33,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,0.67,0.1))}
  pbprob2=function(x){log(0.2*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,1,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.5,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.25,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.75,0.1))}
  CN1=mapply(exactCNi,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,mresult,mLpos,MoreArgs=list(dzprob11=dzprob11,dzprob12=dzprob12,dzprob21=dzprob21,dzprob22=dzprob22,dbprob1=dbprob1,dbprob2=dbprob2,pzprob11=pzprob11,pzprob12=pzprob12,pzprob21=pzprob21,pzprob22=pzprob22,pbprob1=pbprob1,pbprob2=pbprob2),SIMPLIFY = FALSE)
  CN=mapply(function(CNs,mLpos,result,Lpos){
    CNi=result
      for (i in 1:length(CNs)){
        if (CNs[i]!=2){
          s=mLpos[i,1]
          t=mLpos[i,2]
          sel=which(Lpos[,1]>=s & Lpos[,2]<=t)
          CNi[sel]=CNs[i]
        }
      }
    return(CNi)
    },CN1,mLpos,result,Lpos,SIMPLIFY = FALSE)
  return(CN)
}

# Infer exact Copy Number single platform
exactCN_ngs = function(testres,dbafs,dzs,dI,pbafs,pzs,pI,r1L,baf1,rpos1,bpos1){
  result=lapply(testres,function(x){x[[1]]})
  Lpos=lapply(testres,function(x){x[[2]]})
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
  mresult=lapply(res,function(x){x[[1]]})
  mLpos=lapply(res,function(x){x[[2]]})
  # Distribution mean inference using K means
  dz1=dzs[dI==1|dI==1.5]
  if(length(dz1)<3){
    dzprob11=function(x){log(0.99)}
    dzprob12=function(x){log(0.01)}
  }else{
    dz1fit <- kmeans(dz1, 2)
    dmu1=aggregate(dz1,by=list(dz1fit$cluster),FUN=mean)$x
    dzprob11=function(x){dnorm(x,max(dmu1),1,log=TRUE)}
    dzprob12=function(x){dnorm(x,min(dmu1),1,log=TRUE)}
  }
  dbprob1=function(x){log(0.5*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.5*truncnorm::dtruncnorm(x,0,1,1,0.1))}
  dbprob2=function(x){log(0.45*truncnorm::dtruncnorm(x,0,1,0,0.01)+0.45*truncnorm::dtruncnorm(x,0,1,1,0.01)+0.1*dunif(x,0,1,log=TRUE))}
  
  pz1=pzs[pI==1|pI==1.5]
  if(length(pz1)<3){
    pzprob11=function(x){log(0.99)}
    pzprob12=function(x){log(0.01)}
  }else{
    pz1fit <- kmeans(pz1, 2)
    pmu1=aggregate(pz1,by=list(pz1fit$cluster),FUN=mean)$x
    pzprob11=function(x){dnorm(x,min(pmu1),1,log=TRUE)}
    pzprob12=function(x){dnorm(x,max(pmu1),1,log=TRUE)}    
  }
  pbprob1=function(x){log(0.25*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,1,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,0.33,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,0.67,0.1))}
  pbprob2=function(x){log(0.2*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,1,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.5,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.25,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.75,0.1))}
  CN1=mapply(exactCNi_ngs,r1L,baf1,rpos1,bpos1,mresult,mLpos,MoreArgs=list(dzprob11=dzprob11,dzprob12=dzprob12,dbprob1=dbprob1,dbprob2=dbprob2,pzprob11=pzprob11,pzprob12=pzprob12,pbprob1=pbprob1,pbprob2=pbprob2),SIMPLIFY = FALSE)
  CN=mapply(function(CNs,mLpos,result,Lpos){
    CNi=result
      for (i in 1:length(CNs)){
        if (CNs[i]!=2){
          s=mLpos[i,1]
          t=mLpos[i,2]
          sel=which(Lpos[,1]>=s & Lpos[,2]<=t)
          CNi[sel]=CNs[i]
        }
      }
    return(CNi)
    },CN1,mLpos,result,Lpos,SIMPLIFY = FALSE)
  return(CN)
}

# Likelihood estimator
exactCNi_ngs = function(r1i,baf1i,rpos1i,bpos1i,resulti,Lposi,dzprob11,dzprob12,dbprob1,dbprob2,pzprob11,pzprob12,pbprob1,pbprob2){
  CNi=resulti
  z1i=r1i#(r1i-mean(r1i,na.rm=TRUE))/sd(r1i,na.rm=TRUE)
  baf1i[is.na(baf1i)]=1
  z1i[is.na(z1i)]=0
  for (j in 1:length(resulti)){
    if (CNi[j]<2){
      s=Lposi[j,1]
      t=Lposi[j,2]
      del1=0
      del2=0
      sel1=which(rpos1i[,1]>=s&rpos1i[,2]<=t)
      sel2=which(bpos1i>=s&bpos1i<=t)
      if (length(sel2)>0){
        del1=del1+sum(dzprob11(z1i[sel1]))+sum(dbprob1(baf1i[sel2]))
        del2=del2+sum(dzprob12(z1i[sel1]))+sum(dbprob2(baf1i[sel2]))
      }
      else{
        del1=del1+sum(dzprob11(z1i[sel1]))+sum(dbprob1(1))
        del2=del2+sum(dzprob12(z1i[sel1]))+sum(dbprob2(1))
      }
      if (del2>del1){
        CNi[j]=0
      }
    }
    else if(CNi[j]>2){
      s=Lposi[j,1]
      t=Lposi[j,2]
      dup1=0
      dup2=0
      sel1=which(rpos1i[,1]>=s&rpos1i[,2]<=t)
      sel2=which(bpos1i>=s&bpos1i<=t)
      if (length(sel2)>0){
        dup1=dup1+sum(pzprob11(z1i[sel1]))+sum(pbprob1(baf1i[sel2]))
        dup2=dup2+sum(pzprob12(z1i[sel1]))+sum(pbprob2(baf1i[sel2]))
      }
      else{
        dup1=dup1+sum(pzprob11(z1i[sel1]))+sum(pbprob1(1))
        dup2=dup2+sum(pzprob12(z1i[sel1]))+sum(pbprob2(1))
      }
      if (dup2>dup1){
        CNi[j]=4
      }
    }
  }
  return (CNi)
}


# Likelihood estimator
exactCNi = function(r1i,r2i,baf1i,baf2i,rpos1i,rpos2i,bpos1i,bpos2i,resulti,Lposi,dzprob11,dzprob12,dzprob21,dzprob22,dbprob1,dbprob2,pzprob11,pzprob12,pzprob21,pzprob22,pbprob1,pbprob2){
  CNi=resulti
  z1i=r1i#(r1i-mean(r1i,na.rm=TRUE))/sd(r1i,na.rm=TRUE)
  z2i=r2i#(r2i-mean(r2i,na.rm=TRUE))/sd(r2i,na.rm=TRUE)
  baf1i[is.na(baf1i)]=1
  baf2i[is.na(baf2i)]=1
  z1i[is.na(z1i)]=0
  z2i[is.na(z2i)]=0
  for (j in 1:length(resulti)){
    if (CNi[j]<2){
      s=Lposi[j,1]
      t=Lposi[j,2]
      if (s==t){
        sel1=which(rpos2i==s)
        sel2=which(bpos2i==s)
        if ((dzprob21(z2i[sel1])+dbprob1(baf2i[sel2]))<(dzprob22(z2i[sel1])+dbprob2(baf2i[sel2]))){
          CNi[j]=0
        }
      }
      else{
        del1=0
        del2=0
        sel1=which(rpos2i>=s&rpos2i<=t)
        sel2=which(bpos2i>=s&bpos2i<=t)
        if (length(sel1)>0){
          del1=sum(dzprob21(z2i[sel1]))+sum(dbprob1(baf2i[sel2]))
          del2=sum(dzprob22(z2i[sel1]))+sum(dbprob2(baf2i[sel2]))
        }
        sel1=which(rpos1i[,1]>=s&rpos1i[,2]<=t)
        sel2=which(bpos1i>=s&bpos1i<=t)
        if (length(sel2)>0){
          del1=del1+sum(dzprob11(z1i[sel1]))+sum(dbprob1(baf1i[sel2]))
          del2=del2+sum(dzprob12(z1i[sel1]))+sum(dbprob2(baf1i[sel2]))
        }
        else{
          del1=del1+sum(dzprob11(z1i[sel1]))+sum(dbprob1(1))
          del2=del2+sum(dzprob12(z1i[sel1]))+sum(dbprob2(1))
        }
        if (del2>del1){
          CNi[j]=0
        }
      }
    }
    else if(CNi[j]>2){
      s=Lposi[j,1]
      t=Lposi[j,2]
      if (s==t){
        sel1=which(rpos2i==s)
        sel2=which(bpos2i==s)
        if ((pzprob21(z2i[sel1])+pbprob1(baf2i[sel2]))<(pzprob22(z2i[sel1])+pbprob2(baf2i[sel2]))){
          CNi[j]=4
        }
      }
      else{
        dup1=0
        dup2=0
        sel1=which(rpos2i>=s&rpos2i<=t)
        sel2=which(bpos2i>=s&bpos2i<=t)
        if (length(sel1)>0){
          dup1=sum(pzprob21(z2i[sel1]))+sum(pbprob1(baf2i[sel2]))
          dup2=sum(pzprob22(z2i[sel1]))+sum(pbprob2(baf2i[sel2]))
        }
        sel1=which(rpos1i[,1]>=s&rpos1i[,2]<=t)
        sel2=which(bpos1i>=s&bpos1i<=t)
        if (length(sel2)>0){
          dup1=dup1+sum(pzprob11(z1i[sel1]))+sum(pbprob1(baf1i[sel2]))
          dup2=dup2+sum(pzprob12(z1i[sel1]))+sum(pbprob2(baf1i[sel2]))
        }
        else{
          dup1=dup1+sum(pzprob11(z1i[sel1]))+sum(pbprob1(1))
          dup2=dup2+sum(pzprob12(z1i[sel1]))+sum(pbprob2(1))
        }
        if (dup2>dup1){
          CNi[j]=4
        }
      }
    }
  }
  return (CNi)
}


visualization2 = function(testres,CNV,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2){
  cat('\n','Generating visualization2. This may take a while...','\n')
  Lpos=lapply(testres,function(x){x[[2]]})
  ttldipposlist=mapply(function(x){which(x==2)},CNV,SIMPLIFY=FALSE)
  ttldel1poslist=mapply(function(x){which(x==1)},CNV,SIMPLIFY=FALSE)
  ttldel2poslist=mapply(function(x){which(x==0)},CNV,SIMPLIFY=FALSE)
  ttldup1poslist=mapply(function(x){which(x==3)},CNV,SIMPLIFY=FALSE)
  ttldup2poslist=mapply(function(x){which(x==4)},CNV,SIMPLIFY=FALSE)
  nzvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldipposlist,Lpos,SIMPLIFY = FALSE)
  nbafs=unlist(mapply(function(x){x[,2]},nzvsbaf))
  nbafs[is.na(nbafs)]=0
  nzs=unlist(mapply(function(x){x[,1]},nzvsbaf))
  nI=unlist(mapply(function(x){x[,3]},nzvsbaf))
  cat('step1 of 5','\n')
  d1zvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldel1poslist,Lpos,SIMPLIFY = FALSE)
  d1bafs=unlist(mapply(function(x){x[,2]},d1zvsbaf))
  d1bafs[is.na(d1bafs)]=0
  d1zs=unlist(mapply(function(x){x[,1]},d1zvsbaf))
  d1I=unlist(mapply(function(x){x[,3]},d1zvsbaf))
  cat('step2 of 5','\n')
  d2zvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldel2poslist,Lpos,SIMPLIFY = FALSE)
  d2bafs=unlist(mapply(function(x){x[,2]},d2zvsbaf))
  d2bafs[is.na(d2bafs)]=0
  d2zs=unlist(mapply(function(x){x[,1]},d2zvsbaf))
  d2I=unlist(mapply(function(x){x[,3]},d2zvsbaf))
  cat('step3 of 5','\n')
  p1zvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldup1poslist,Lpos,SIMPLIFY = FALSE)
  p1bafs=unlist(mapply(function(x){x[,2]},p1zvsbaf))
  p1bafs[is.na(p1bafs)]=0
  p1zs=unlist(mapply(function(x){x[,1]},p1zvsbaf))
  p1I=unlist(mapply(function(x){x[,3]},p1zvsbaf))
  cat('step4 of 5','\n')
  p2zvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldup2poslist,Lpos,SIMPLIFY = FALSE)
  p2bafs=unlist(mapply(function(x){x[,2]},p2zvsbaf))
  p2bafs[is.na(p2bafs)]=0
  p2zs=unlist(mapply(function(x){x[,1]},p2zvsbaf))
  p2I=unlist(mapply(function(x){x[,3]},p2zvsbaf))
  cat('step5 of 5','\n')
  Is=c(nI,d1I,d2I,p1I,p2I)
  CNV=c(rep(2,length(nI)),rep(1,length(d1I)),rep(0,length(d2I)),rep(3,length(p1I)),rep(4,length(p2I)))
  d.f=data.frame(zs=c(nzs,d1zs,d2zs,p1zs,p2zs),bafs=c(nbafs,d1bafs,d2bafs,p1bafs,p2bafs),Is=factor(c(nI,d1I,d2I,p1I,p2I),labels=c("SNPs only", "SNPs in Exon", "Exon w/ BAFs", "Exon w/ no BAF")),CNV=c(rep(2,length(nI)),rep(1,length(d1I)),rep(0,length(d2I)),rep(3,length(p1I)),rep(4,length(p2I))))
  # save(d.f,file='ggplot_df.rda')
  p=ggplot2::ggplot(d.f,ggplot2::aes(x=zs,y=bafs))+ggplot2::geom_point(ggplot2::aes(colour=factor(CNV),shape=Is),alpha=0.75)+
    ggplot2::scale_colour_manual(name="", values = c("0"="blue", "1"="cyan", "2"="grey","3"="orange", "4"="red"))+ggplot2::stat_density2d(ggplot2::aes(colour=factor(CNV),linetype=Is),h=c(8,0.3))
  print(p)
  cat('\n','Finish generating visualization2 plot...','\n')
}

novisualization = function(testres,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2){
  cat('\n','Calculating novisualization...','\n')
  result=lapply(testres,function(x){x[[1]]})
  Lpos=lapply(testres,function(x){x[[2]]})
  ttldelposlist=mapply(function(x){which(x==1)},result,SIMPLIFY=FALSE)
  ttldupposlist=mapply(function(x){which(x==3)},result,SIMPLIFY=FALSE)
  dzvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldelposlist,Lpos,SIMPLIFY = FALSE)
  dbafs=unlist(mapply(function(x){x[,2]},dzvsbaf))
  dbafs[is.na(dbafs)]=0
  dzs=unlist(mapply(function(x){x[,1]},dzvsbaf))
  dI=unlist(mapply(function(x){x[,3]},dzvsbaf))
  pzvsbaf=mapply(BAFvsZ,r1L,r2L,baf1,baf2,rpos1,rpos2,bpos1,bpos2,ttldupposlist,Lpos,SIMPLIFY = FALSE)
  pbafs=unlist(mapply(function(x){x[,2]},pzvsbaf))
  pbafs[is.na(pbafs)]=0
  pzs=unlist(mapply(function(x){x[,1]},pzvsbaf))
  pI=unlist(mapply(function(x){x[,3]},pzvsbaf))
  x.lim=range(c(range(pzs,na.rm=TRUE),range(dzs,na.rm=TRUE)))
  return(list(dbafs,dzs,dI,pbafs,pzs,pI,x.lim))
}

novisualization_ngs = function(testres,r1L,baf1,rpos1,bpos1){
  cat('\n','Calculating novisualization...','\n')
  result=lapply(testres,function(x){x[[1]]})
  Lpos=lapply(testres,function(x){x[[2]]})
  ttldelposlist=mapply(function(x){which(x==1)},result,SIMPLIFY=FALSE)
  ttldupposlist=mapply(function(x){which(x==3)},result,SIMPLIFY=FALSE)
  dzvsbaf=mapply(BAFvsZ_ngs,r1L,baf1,rpos1,bpos1,ttldelposlist,Lpos,SIMPLIFY = FALSE)
  dbafs=unlist(mapply(function(x){x[,2]},dzvsbaf))
  dbafs[is.na(dbafs)]=0
  dzs=unlist(mapply(function(x){x[,1]},dzvsbaf))
  dI=unlist(mapply(function(x){x[,3]},dzvsbaf))
  pzvsbaf=mapply(BAFvsZ_ngs,r1L,baf1,rpos1,bpos1,ttldupposlist,Lpos,SIMPLIFY = FALSE)
  pbafs=unlist(mapply(function(x){x[,2]},pzvsbaf))
  pbafs[is.na(pbafs)]=0
  pzs=unlist(mapply(function(x){x[,1]},pzvsbaf))
  pI=unlist(mapply(function(x){x[,3]},pzvsbaf))
  x.lim=range(c(range(pzs,na.rm=TRUE),range(dzs,na.rm=TRUE)))
  return(list(dbafs,dzs,dI,pbafs,pzs,pI,x.lim))
}

# Pair intensity and BAF, NGS
BAFvsZ_ngs=function(r1i,baf1i,rpos1i,bpos1i,CNVs,Lposi){
  bafs=c()
  rs=c()
  I=c()
  if (length(CNVs)>0){
    z1i=(r1i-mean(r1i,na.rm=TRUE))/sd(r1i,na.rm=TRUE)
    mbaf1i=0.5-abs(baf1i-0.5)
    for (j in 1:length(CNVs)){
      s=Lposi[CNVs[j],1]
      t=Lposi[CNVs[j],2]
      sel1=which(rpos1i[,1]>=s&rpos1i[,2]<=t)
      sel2=which(bpos1i>=s&bpos1i<=t)
      if (length(sel2)>0){
        bafs=c(bafs,mbaf1i[sel2])
        rs=c(rs,rep(z1i[sel1],length(sel2)))
        I=c(I,rep(1,length(sel2)))
      }
      else{
        bafs=c(bafs,0)
        rs=c(rs,z1i[sel1])
        I=c(I,1.5)
      }

    }
  }
  return(cbind(rs,bafs,I))
}

# Pair intensity and BAF
BAFvsZ_snp=function(r2i,baf2i,rpos2i,bpos2i,CNVs,Lposi){
  bafs=c()
  rs=c()
  I=c()
  if (length(CNVs)>0){
    z2i=(r2i-mean(r2i,na.rm=TRUE))/sd(r2i,na.rm=TRUE)
    mbaf2i=0.5-abs(baf2i-0.5)
    for (j in 1:length(CNVs)){
      s=Lposi[CNVs[j],1]
      t=Lposi[CNVs[j],2]
      sel1=which(rpos2i==s)
      sel2=which(bpos2i==s)
      bafs=c(bafs,mbaf2i[sel2])
      rs=c(rs,z2i[sel1])
      I=c(I,2)
    }
  }
  return(cbind(rs,bafs,I))
}



novisualization_snp = function(testres,r2L,baf2,rpos2,bpos2){
  cat('\n','Calculating novisualization...','\n')
  result=lapply(testres,function(x){x[[1]]})
  Lpos=lapply(testres,function(x){x[[2]]})
  ttldelposlist=mapply(function(x){which(x==1)},result,SIMPLIFY=FALSE)
  ttldupposlist=mapply(function(x){which(x==3)},result,SIMPLIFY=FALSE)
  dzvsbaf=mapply(BAFvsZ_snp,r2L,baf2,rpos2,bpos2,ttldelposlist,Lpos,SIMPLIFY = FALSE)
  dbafs=unlist(mapply(function(x){x[,2]},dzvsbaf))
  dbafs[is.na(dbafs)]=0
  dzs=unlist(mapply(function(x){x[,1]},dzvsbaf))
  dI=unlist(mapply(function(x){x[,3]},dzvsbaf))
  pzvsbaf=mapply(BAFvsZ_snp,r2L,baf2,rpos2,bpos2,ttldupposlist,Lpos,SIMPLIFY = FALSE)
  pbafs=unlist(mapply(function(x){x[,2]},pzvsbaf))
  pbafs[is.na(pbafs)]=0
  pzs=unlist(mapply(function(x){x[,1]},pzvsbaf))
  pI=unlist(mapply(function(x){x[,3]},pzvsbaf))
  x.lim=range(c(range(pzs,na.rm=TRUE),range(dzs,na.rm=TRUE)))
  return(list(dbafs,dzs,dI,pbafs,pzs,pI,x.lim))
}


# Infer exact Copy Number SNP
exactCN_snp = function(testres,dbafs,dzs,dI,pbafs,pzs,pI,r2L,baf2,rpos2,bpos2){
  result=lapply(testres,function(x){x[[1]]})
  Lpos=lapply(testres,function(x){x[[2]]})
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
  mresult=lapply(res,function(x){x[[1]]})
  mLpos=lapply(res,function(x){x[[2]]})
  # Distribution mean inference using K means
  dz2=dzs[dI==2|dI==2.5]
  if(length(dz2)<3){
    dzprob21=function(x){log(0.99)}
    dzprob22=function(x){log(0.01)}
  }else{
    dz2fit <- kmeans(dz2, 2)
    dmu2=aggregate(dz2,by=list(dz2fit$cluster),FUN=mean)$x
    dzprob21=function(x){dnorm(x,max(dmu2),1,log=TRUE)}
    dzprob22=function(x){dnorm(x,min(dmu2),1,log=TRUE)}
  }
  dbprob1=function(x){log(0.5*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.5*truncnorm::dtruncnorm(x,0,1,1,0.1))}
  dbprob2=function(x){log(0.45*truncnorm::dtruncnorm(x,0,1,0,0.01)+0.45*truncnorm::dtruncnorm(x,0,1,1,0.01)+0.1*dunif(x,0,1,log=TRUE))}
  pz2=pzs[pI==2|pI==2.5]
  if(length(pz2)<3){
    pzprob21=function(x){log(0.99)}
    pzprob22=function(x){log(0.01)}
  }else{
    pz2fit <- kmeans(pz2, 2)
    pmu2=aggregate(pz2,by=list(pz2fit$cluster),FUN=mean)$x
    pzprob21=function(x){dnorm(x,min(pmu2),1,log=TRUE)}
    pzprob22=function(x){dnorm(x,max(pmu2),1,log=TRUE)}
  }
  pbprob1=function(x){log(0.25*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,1,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,0.33,0.1)+0.25*truncnorm::dtruncnorm(x,0,1,0.67,0.1))}
  pbprob2=function(x){log(0.2*truncnorm::dtruncnorm(x,0,1,0,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,1,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.5,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.25,0.1)+0.2*truncnorm::dtruncnorm(x,0,1,0.75,0.1))}
  CN1=mapply(exactCNi_snp,r2L,baf2,rpos2,bpos2,mresult,mLpos,MoreArgs=list(dzprob21=dzprob21,dzprob22=dzprob22,dbprob1=dbprob1,dbprob2=dbprob2,pzprob21=pzprob21,pzprob22=pzprob22,pbprob1=pbprob1,pbprob2=pbprob2),SIMPLIFY = FALSE)
  CN=mapply(function(CNs,mLpos,result,Lpos){
    CNi=result
      for (i in 1:length(CNs)){
        if (CNs[i]!=2){
          s=mLpos[i,1]
          t=mLpos[i,2]
          sel=which(Lpos[,1]>=s & Lpos[,2]<=t)
          CNi[sel]=CNs[i]
        }
      }
    return(CNi)
    },CN1,mLpos,result,Lpos,SIMPLIFY = FALSE)
  return(CN)
}


# Likelihood estimator
exactCNi_snp = function(r2i,baf2i,rpos2i,bpos2i,resulti,Lposi,dzprob21,dzprob22,dbprob1,dbprob2,pzprob21,pzprob22,pbprob1,pbprob2){
  CNi=resulti
  z2i=r2i#(r2i-mean(r2i,na.rm=TRUE))/sd(r2i,na.rm=TRUE)
  baf2i[is.na(baf2i)]=1
  z2i[is.na(z2i)]=0
  for (j in 1:length(resulti)){
    if (CNi[j]<2){
      s=Lposi[j,1]
      t=Lposi[j,2]
      sel1=which(rpos2i==s)
      sel2=which(bpos2i==s)
      if ((dzprob21(z2i[sel1])+dbprob1(baf2i[sel2]))<(dzprob22(z2i[sel1])+dbprob2(baf2i[sel2]))){
        CNi[j]=0
      }
    }
    else if(CNi[j]>2){
      s=Lposi[j,1]
      t=Lposi[j,2]
      sel1=which(rpos2i==s)
      sel2=which(bpos2i==s)
      if ((pzprob21(z2i[sel1])+pbprob1(baf2i[sel2]))<(pzprob22(z2i[sel1])+pbprob2(baf2i[sel2]))){
        CNi[j]=4
      }
    }
  }
  return (CNi)
}




