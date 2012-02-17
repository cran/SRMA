calc.post <- function(d,prior,f1i,f3i,error){
                      ##strand.spec=FALSE){
  ##if(strand.spec==FALSE){
  post <- matrix(NA,nrow(d),3)
  for(i in 1:nrow(d)){
    ll1 <- ll2 <- ll3 <- 0
    for(j in 1:2){
      ll1 <- ll1+calc.ll(d[i,j],f1i[i,j],error[i,j])
      ll2 <- ll2+calc.ll(d[i,j],0,error[i,j])
      ll3 <- ll3+calc.ll(d[i,j],f3i[i,j],error[i,j])
    }
    log.denom <- c(ll1,ll2,ll3)+log(prior)
    post[i,] <- exp(log.denom)/sum(exp(log.denom), na.rm=T)
  }
  ##}else{
  ##  post <- matrix(NA,nrow(d),3)
  ##  for(i in 1:nrow(d)){
  ## ll1 <- ll2 <- ll3 <- NULL
  ## for(j in 1:2){
  ##   ll1 <- c(ll1,calc.ll(d[i,j],f1i[i,j],error[i,j]))
  ##   ll2 <- c(ll2,calc.ll(d[i,j],0,error[i,j]))
  ##   ll3 <- c(ll3,calc.ll(d[i,j],f3i[i,j],error[i,j]))
  ## }
  ## log.denom <- c(sum(ll1),sum(ll2),sum(ll3))+log(prior)
  ## post[i,] <- exp(log.denom)/sum(exp(log.denom))
  ##}
  return(post)
}
