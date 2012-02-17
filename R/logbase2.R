logbase2 <- function(x){
  y <- rep(0,length(x))
  y[x>0] <- log2(x[x>0])
  return(y)
}
