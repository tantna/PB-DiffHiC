exactTestPoisson <- function(dataMatrix, meanMatrix, group1Ind, group2Ind, verbose=TRUE) {
  
  dataMatrix = round(dataMatrix)
  y1 <- rowSums(dataMatrix[,group1Ind])
  y2 <- rowSums(dataMatrix[,group2Ind])
  m1 <- rowSums(meanMatrix[,group1Ind])
  m2 <- rowSums(meanMatrix[,group2Ind])
  
  N <- rowSums( dataMatrix[,c(group1Ind,group2Ind)] )
  
  pvals <- rep(NA, nrow(dataMatrix))
  
  for (i in 1:length(pvals)) {
    v <- 0:N[i]
    p.top <- dpois(v, lambda=m1[i]) * dpois((N[i]-v), lambda=m2[i])
    p.obs <- dpois(y1[i], lambda=m1[i]) * dpois(y2[i], lambda=m2[i])
    p.bot <- dpois(N[i], lambda=m1[i]+m2[i])
    keep <- p.top <= p.obs
    pvals[i] <- sum(p.top[keep]/p.bot)
    if(N[i]>10000) {	
      pvals[i]<- chisq.test(matrix(c(y1[i],y2[i],(sum(y1)-y1[i]),(sum(y2)-y2[i])),2,2))$p.value
    } 
    if (verbose)
      if (i%%1000 == 0)
        cat(".")
  }
  if (verbose)
    cat("\n")
  # pvals[N==0]<-0  
  pvals
}
#------scbnm--------------
calcNormFactors <- function(dataMatrix, refColumn=1, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {
  if( !is.matrix(dataMatrix) )
    stop("'dataMatrix' needs to be a matrix")
  if( refColumn > ncol(dataMatrix) )
    stop("Invalid 'refColumn' argument")
  apply(dataMatrix,2,.calcFactorWeighted,ref=dataMatrix[,refColumn], logratioTrim=logratioTrim, 
        sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
}
.calcFactorWeighted <- function(obs, ref, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {
  
  if( all(obs==ref) )
    return(1)
  
  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
  
  # remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  if (doWeighting) 
    2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  else
    2^( mean(logR[keep], na.rm=TRUE) )
}
sage.test2 <- function(x, y, n1=sum(x), n2=sum(y))
  #	Binomial probabilities for comparing SAGE libraries
  #	Gordon Smyth
  #	15 Nov 2003.  Last modified 21 Jan 2004.
{
  if(any(is.na(x)) || any(is.na(y))) stop("missing values not allowed")
  nn<-sum(x)*sqrt(n1/n2)
  mm<-sum(y)*sqrt(n2/n1)
  x <- round(x)
  y <- round(y)
  if(any(x<0) || any(y<0)) stop("x and y must be non-negative")
  if(length(x) != length(y)) stop("x and y must have same length")
  n1 <- round(n1)
  n2 <- round(n2)
  #if(!missing(n1) && any(x>n1)) stop("x cannot be greater than n1")
  #if(!missing(n2) && any(y>n2)) stop("y cannot be greater than n2")
  size <- x+y
  p.value <- rep(1,length(x))
  if(n1==n2) {
    i <- (size>0)
    if(any(i)) {
      x <- pmin(x[i],y[i])
      size <- size[i]
      p.value[i] <- pbinom(x,size=size,prob=0.5)+pbinom(size-x+0.5,size=size,prob=0.5,lower.tail=FALSE)
    }
    return(p.value)
  }
  prob <- n1/(n1+n2)
  nn <- round(nn)
  mm <- round(mm)
  
  if(any(big <- size>10000)) {
    ibig <- (1:length(x))[big]
    for (i in ibig) p.value[i] <- chisq.test(matrix(c(x[i],y[i],(nn-x[i]),(mm-y[i])),2,2))$p.value
  }
  size0 <- size[size>0 & !big]
  if(length(size0)) for (isize in unique(size0)) {
    i <- (size==isize)
    p <- dbinom(0:isize,p=prob,size=isize)
    o <- order(p)
    cumsump <- cumsum(p[o])[order(o)]
    p.value[i] <- cumsump[x[i]+1]
  }
  p.value
}


#(1) SBNM method to choose normalization factor.
SCBNM<-function(xx,housekeep,a=0.05)
{
  library(edgeR)
  library(statmod)
  fW <- calcNormFactors(xx,logratioTrim=0.3, sumTrim=0.05)[2]
  fW4<-seq(fW-0.5,fW+0.5,0.1)
  q4<-rep(0,length(xx[,1]))
  n<-length(fW4)
  fdr<-rep(0,n)
  for(j in 1:n){
    sMm4<-sage.test2(xx[,1], xx[,2], n1=sum(xx[,1])/sqrt(fW4[j]), n2=sum(xx[,2])*sqrt(fW4[j]))
    #q4<-sMm4*length(xx[,2])/rank(sMm4)
    fdr[j]<-sum(sMm4[housekeep]<a)
  }
  fw5<-fW4[which.min(abs(fdr-a))]
  fW41<-seq(fw5-0.25,fw5+0.25,0.005)
  q41<-rep(0,length(xx[,1]))
  n1<-length(fW41)
  fdr1<-rep(0,n1)
  for(j in 1:n1){
    sMm41<-sage.test2(xx[,1], xx[,2], n1=sum(xx[,1])/sqrt(fW41[j]), n2=sum(xx[,2])*sqrt(fW41[j]))
    #q41<-sMm41*length(xx[,2])/rank(sMm41)
    fdr1[j]<-sum(sMm41[housekeep]<a)
    if(j %% 10==0) cat(".")
  }
  fw51<-fW41[which.min(abs(fdr1-a))]
  return(fw51)
}

cal_scale=function(Datalist,Hkind){
  #Create the input data to calculate the scaling factor
  vec_hk=DiagCount=list()
  for (i in seq_along(Datalist)) {
    DiagCount[[i]]=Datalist[[i]][1:Hkind]
    vec_hk[[i]] = c(DiagCount[[i]],Datalist[[i]])
  }
  ## Calculate scaling factor for data
  factor_vec = c(1)
  n <- length(vec_hk[[1]])
  
  for (m in 2:4){
    fx = cbind(vec_hk[[1]], vec_hk[[m]])
    
    # starttime = Sys.time()
    factor_scbn <- SCBNM(fx, housekeep = 1:Hkind, a=0.05) 
    # endtime=Sys.time()
    # t=endtime-starttime
    factor_vec = c(factor_vec,factor_scbn)
  }
  print(factor_vec)
  return(factor_vec)
}

PB_two <- function(Datalist,Scale){
  vec <- Datalist
  factor_vec <- Scale
  
  ## Calculate p-values and select significants
  # scale_scbn = scale_scbn
  xx = cbind(vec[[1]],vec[[2]],vec[[3]],vec[[4]])
  
  lambda <- rowSums(xx)/sum(xx)
  cs <- colSums(xx)
  effM <- cs*factor_vec
  expMeanAdj <- outer(lambda, effM)
  exactPadj <-exactTestPoisson(dataMatrix = xx, group1Ind=1:2, group2Ind=3:4,  
                               meanMatrix = expMeanAdj, verbose = TRUE)

  p_adjust <- p.adjust(exactPadj, method = "BH")
  return(list(pv = exactPadj,qv = p_adjust,scale = factor_vec))
}

