ssum <- function(x1,x2){
  a <- as.vector(x1)
  b <- as.vector(x2)

  tmp1 <- is.na(a)&!is.na(b) # if a[i] has missing value and b[i] has no missing value,tmp1[i]=TRUE, otherwise tmp1[i]=FALSE.(i=1,...,length(a))
  tmp2 <- !is.na(a)&is.na(b)
  tmp3 <- !is.na(a)&!is.na(b)

  x3 <- rep(NA,length(a))
  x3[tmp1] <- b[tmp1]
  x3[tmp2] <- a[tmp2]
  x3[tmp3] <- (a[tmp3]+b[tmp3])/2

  dim(x3) <- dim(x1)
  return(x3)
}

