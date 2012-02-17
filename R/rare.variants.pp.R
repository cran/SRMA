rare.variants.pp <- function(d,ai,bad,center,cov,prior,dbsnpID){

  ##assign a smaller covariance to the het cluster
  posterior <- matrix(NA,nrow(d),3)
  for(k in 1:sum(!bad)){
    ll1 <- c(calc.ll(d[!bad,][k,1],center[1],error=sqrt(cov[1,1])),
             calc.ll(d[!bad,][k,2],center[2],error=sqrt(cov[2,2])))
    ll2 <- c(calc.ll(d[!bad,][k,1],0,error=sqrt(cov[1,1])),
             calc.ll(d[!bad,][k,2],0,error=sqrt(cov[2,2])))
   ## ll2 <- c(calc.ll(d[!bad,][k,1],0,error=0.1),
   ##          calc.ll(d[!bad,][k,2],0,error=0.1))
    ll3 <- c(calc.ll(d[!bad,][k,1],-1.7,error=sqrt(cov[1,1])),
             calc.ll(d[!bad,][k,2],-1.7,error=sqrt(cov[2,2])))
    s <- c(ll1[1],ll2[1],ll3[1])
    as <- c(ll1[2],ll2[2],ll3[2])
  ##  if(!identical(order(s),order(as))){

    ####begin, change on 051311
    if(order(s)[3]<order(as)[3]&order(s)[3]==1){
      ##comb.s <- s
      if(length(dbsnpID)>0){
        comb.s <- as
      }else{
        comb.s <- s
      }
    }
    if(order(s)[3]>order(as)[3]&order(as)[3]==1){
    ##  comb.s <- as
      if(length(dbsnpID)>0){
        comb.s <- s
      }else{
        comb.s <- as
      }
    }
    
    if((order(s)[3]>1&order(as)[3]>1)|(order(s)[3]==order(as)[3])){
    ####end, change on 051311
      
      ## browser()
      ##only use one strand, if the other strand has very low average intensities
      ai1 <- sum(ai[!bad,1]<9)
      ai2 <- sum(ai[!bad,2]<9)
      ai3 <- sum(!bad)/3
      if(ai1>ai3&ai1>ai2){
        comb.s <- as
      }else{
        if(ai2>ai3&ai1<ai2){
          comb.s <- s
        }else{
          comb.s <- s+as                
        }
      }
    }
  
    ##log.denom <- comb.s+log(prior)
    ##posterior[!bad,][k,] <- exp(log.denom)/sum(exp(log.denom))
    num <- exp(comb.s+log(prior))
    if(sum(num==0)==length(num)){
      num <- rep(1,length(num))
    }
    posterior[!bad,][k,] <- num/sum(num)
  }
  return(posterior)
}
