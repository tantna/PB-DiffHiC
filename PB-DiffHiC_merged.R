library(SCBN)

cal_scale_merged=function(Datalist,Hkind){
  ##combine two matrix as a matrix,to input SCBN
  n <- length(Datalist[[1]])
  Hic_mat <- matrix(0,n,4)
  Hic_mat[,2] <- Datalist[[1]]
  Hic_mat[,4] <- Datalist[[2]]
  Hic_mat[,1] = Hic_mat[,3] = 10000  #gene length
  #print(Hic_mat)
  factor_scbn <- SCBN(orth_gene=Hic_mat, hkind=1:Hkind, a=0.05)
  # factor_scbn
  scale_scbn <- factor_scbn$scbn_val
  return(c(1,scale_scbn))
}

PB_merged <- function(Datalist,Hkind,scale_factor,Scale=TRUE){
  ##combine two matrix as a matrix,to input SCBN
  n <- length(Datalist[[1]])
  Hic_mat <- matrix(0,n,4)
  Hic_mat[,2] <- Datalist[[1]]
  Hic_mat[,4] <- Datalist[[2]]
  Hic_mat[,1] = Hic_mat[,3] = 10000  #gene length
  #print(Hic_mat)
  
  ## Calculate scaling factor for data
  if (Scale){
    # factor_scbn
    scale_scbn <- cal_scale_merged(Datalist,Hkind)[2]
  }else{
    scale_scbn=scale_factor[2]
  }
  print(scale_scbn)
  
  ## Calculate p-values and select significants
  # scale_scbn = scale_scbn
  orth_gene <- Hic_mat
  x <- orth_gene[, 2]
  y <- orth_gene[, 4]
  lengthx <- orth_gene[, 1]
  lengthy <- orth_gene[, 3]
  n1 <- sum(x)
  n2 <- sum(y)
  p_value <- sageTestNew(x, y, lengthx, lengthy, n1, n2, scale_scbn)
  #head(p_value)
  
  p_adjust <- p.adjust(p_value, method = "BH")
  #head(p_adjust)
  return(list(pv = p_value,qv=p_adjust,scale=scale_scbn))
}
