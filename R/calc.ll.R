calc.ll <- function(x,center,error){
  ##only calculates log likelihood for univariate normal distribution
  ll <- -1/2*log(2*pi)-log(error)-1/2*(x-center)^2/(error^2)
  
  ##ll <- -k/2*log(2*pi)-1/2*log(det(variance))-1/2*t(x-center)%*%(solve(variance))%*%(x-center)
  return(ll)
}
