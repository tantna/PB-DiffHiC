library(SCBN)

PB_merged <- function(Data1,Data2,Hkind){
  ##combine two matrix as a matrix,to input SCBN
  n <- length(Data1)
  Hic_mat <- matrix(0,n,4)
  Hic_mat[,2] <- Data1
  Hic_mat[,4] <- Data2
  Hic_mat[,1] = Hic_mat[,3] = 10000  #gene length
  #print(Hic_mat)
  
  ## Calculate scaling factor for data
  factor_scbn <- SCBN(orth_gene=Hic_mat, hkind=1:Hkind, a=0.05)
  # factor_scbn
  scale_scbn <- factor_scbn$scbn_val
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
  return(list(hicmat=Hic_mat,pv = p_value,qv=p_adjust,scale_factor=scale_scbn))
}

