#########
# qn.R
# Wenyi Wang
#########

qn <- function(cels){
    #1. given n arrays of length p, form X of dimension
    # p × n where each array is a column;
    # 2. sort each column of X to give Xsort;
    # 3. take the means across rows of Xsort and assign this
    # mean to each element in the row to get X sort;
    # 4. get Xnormalized by rearranging each column of
    #  X sort to have the same ordering as original X

  order.X <- apply(cels,2,order)
  Xsort <- apply(cels,2,sort)
  mean.Xsort <- apply(Xsort,1,median)

  Xnorm <- matrix(0,nrow(cels),ncol(cels))
  
  for(i in 1:ncol(cels)){
    Xnorm[order.X[,i],i] <- mean.Xsort
  }

  ##to double check
#  order.Xnorm <- apply(Xnorm,2,order)
#  identical(order.Xnorm,order.X)
  return(Xnorm)
}
