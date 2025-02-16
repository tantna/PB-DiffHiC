library(Matrix)
library(mvtnorm)

# gussfilter(mm)
gaussfilter = function(mu,sigma,ksize,rawmat,gauss_dis){
  mean = c(mu,mu)
  s = diag(sigma,2,2)
  c0=ceiling(ksize/2)
  #gM ：gauss matrix
  gM = matrix(0,ksize,ksize)
  
  for (i in 1:ksize){
    for (j in i:ksize){
      gM[i,j]=dmvnorm(c(i-c0,j-c0),mean,s)
    }
  }
  
  gM[lower.tri(gM)] <- t(gM)[lower.tri(gM)]
  gM = gM*(1/dmvnorm(mean,mean,s))
  
  #paddding row matrix
  M1 = as.matrix(rawmat)
  MT=t(M1)
  diag(MT)=0
  M=M1+MT
  #pm = matrix(0,dim(M)+2,dim(M)+2)
  pad_size=as.integer((ksize-1)/2)
  pm = matrix(0,nrow(M)+2*pad_size,nrow(M)+2*pad_size)
  
  pm[(pad_size+1):(nrow(pm)-pad_size),(pad_size+1):(ncol(pm)-pad_size)] = M
  
  pm[(pad_size+1):(nrow(pm)-pad_size),1:pad_size] = M[,1:pad_size]#padding left
  
  pm[(pad_size+1):(nrow(pm)-pad_size),(ncol(pm)-pad_size+1):ncol(pm)] = M[,(ncol(M)-pad_size+1):ncol(M)]#padding right
  
  pm[1:pad_size,(pad_size+1):(ncol(pm)-pad_size)] = M[1:pad_size,]#padding upper
  
  pm[(ncol(pm)-pad_size+1):ncol(pm),(pad_size+1):(ncol(pm)-pad_size)] = M[(ncol(M)-pad_size+1):ncol(M),]#padding down
  
  pm[1:pad_size,1:pad_size] = M[1:pad_size,1:pad_size]
  
  pm[1:pad_size,(ncol(pm)-pad_size+1):ncol(pm)] = M[1:pad_size,(ncol(M)-pad_size+1):ncol(M)]
  
  pm[(ncol(pm)-pad_size+1):ncol(pm),1:pad_size] = M[(ncol(M)-pad_size+1):ncol(M),1:pad_size]
  
  pm[(ncol(pm)-pad_size+1):ncol(pm),(ncol(pm)-pad_size+1):ncol(pm)] = M[(ncol(M)-pad_size+1):ncol(M),(ncol(M)-pad_size+1):ncol(M)]
  
  #gaussian filter
  r = matrix(0,nrow(M),nrow(M))
  for (i in 1:nrow(M)) {
    for (j in i:min(i+gauss_dis,nrow(M))) {
      r[i,j] = sum(pm[(i):(i+2*pad_size),(j):(j+2*pad_size)] * gM)
    }
  }
  return(r)
}


#data process
gauss2vec=function(data,test_dis,msize,keepdis,ksize,mu=0,sigma=1,data2mat=TRUE,gauss=TRUE){
  if (data2mat){
    A=spMatrix(msize,msize,i=data$x1,j=data$y1,x=data$ifs)
  }else{
    A=as.matrix(data)
  }
  if (gauss){
    A=gaussfilter(mu,sigma,ksize,A,test_dis)data
  }
  h_keep=0
  vec=vector();from=vector();to=vector()
  for (k in 1:test_dis) {
    if (k!=msize) {
      Diag <- diag(A[1:(msize-k+1),k:msize])
      if (k<=keepdis){
        h_keep=h_keep+length(Diag)  }
    }
    else {
      Diag <- A[1,k]
    }
    vec <- append(vec,Diag)
    from=append(from,c(1:(msize-k+1)))
    to=append(to,c(k:msize))
  }
  return(list(vec=vec,h_keep=h_keep,from=from,to=to))
}
