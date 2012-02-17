position.post <- function(indx1,ref,het,hom,bp.indx,alldelta,f1,f3,
                          errors1,errors3,round){
  ##  prior <- c((1-af)^2,2*af*(1-af),af^2)
  nsample<-dim(alldelta)[3]        
  prior <- c(0.999,0.0009,0.0001)
  d <- cbind(alldelta[indx1,bp.indx,,1],alldelta[indx1,4-bp.indx,,2])
  if(round==1){
    errors <- errors1
    f1i <- cbind(f1[indx1,bp.indx,,1],f1[indx1,4-bp.indx,,2])
    f3i <- cbind(f3[indx1,bp.indx,,1],f3[indx1,4-bp.indx,,2])
    error <- cbind(errors[indx1,bp.indx,,1],errors[indx1,4-bp.indx,,2])
  }
  if(round==2){
    f1i <- f1[indx1,,]
    f3i <- f3[indx1,,]
    error1 <- errors1[indx1,,]
    error3 <- errors3[indx1,,]
  }
  ##if(round==1){
  ## rf1i <- as.matrix(d)-as.matrix(f1i)+cbind(mu1,mu1)
  ## rf2i <- as.matrix(d)
  ## rf3i <- as.matrix(d)-as.matrix(f3i)+cbind(mu3,mu3)
  ##
  ##if(round==2){
  ## rf1i <- as.matrix(d)-as.matrix(f1i)+as.matrix(mu1)
  ## rf2i <- as.matrix(d)
  ## rf3i <- as.matrix(d)-as.matrix(f3i)+as.matrix(mu3)
  ##
    
  if(round==1){
    post <- calc.post(d,prior,f1i,f3i,error)
    ##note a posterior prob<0.9985 indicates non-ref
  }
  if(round==2){
    post <- calc.post.2(d,prior,f1i,f3i,error1,error3)
  }
  genotypes <- rep(NA,nsample)
  genotypes[post[,1]>=prior[1]] <- ref
  genotypes[post[,3]>0.1&post[,2]<0.8] <- hom
  genotypes[is.na(genotypes)] <- het

  output <- list(post,genotypes)
  names(output) <- c("post","genotypes")
  return(output)
}
