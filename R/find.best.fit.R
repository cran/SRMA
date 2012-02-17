find.best.fit <- function(bic,allfit){
  ##for unadjusted fits only
  ## browser()
  while(sum(!is.na(bic))>0){
    x <- order(bic)[1]
    if(x==3){
      break
    }
    fit <- allfit[[x]]
    m <- t(fit$parameters$mean)
    m <- m[order(m[,1],decreasing=TRUE),]
    if(x==1){
      a1 <- sum(m[1,]>0)==2
      a2 <- sum(m[1,]>m[2,])==2&sum(m[2,]>m[3,])==2
      a3 <- sum(m[3,]<0)==2
      a4 <- sqrt(sum((m[1,]-m[2,])^2))>sqrt(sum((m[2,]-c(0,0))^2))|fit$parameters$pro[3]>=0.1
      a5 <- sqrt(sum((m[3,]-m[2,])^2))>sqrt(sum((m[2,]-c(0,0))^2))|fit$parameters$pro[3]>=0.1
      a6 <- sqrt(sum((m[1,]-c(0,0))^2))>sqrt(sum((m[2,]-c(0,0))^2))|sqrt(sum((m[3,]-c(0,0))^2))>sqrt(sum((m[2,]-c(0,0))^2))
      if(!(a1&a2&a3&a4&a5&a6)){
        bic[x] <- NA
      }else{
        break
      }
    }
    if(x==2){
      a1 <- sum(m[1,]>m[2,])==2
      d1 <- sqrt(sum((m[1,]-c(0,0))^2))
      d2 <- sqrt(sum((m[2,]-c(0,0))^2))
      if(d1==d2){
        break
      }
      ####begin, change on 051311
      if(d1>d2){
        a2 <- sum(m[1,]>0)==2
        a3 <- sqrt(sum((m[1,]-m[2,])^2))>d2
        ##a4 <- sqrt(sum((m[1,]-c(0,0))^2))>d2
      }
      if(d1<d2){
        a2 <- sum(m[2,]<0)==2
        a3 <- sqrt(sum((m[2,]-m[1,])^2))>d1
        ##a4 <- sqrt(sum((m[2,]-c(0,0))^2))>d1
      }
      ##if(!(a1&a2&a3&a4)){
      if(!(a1&a2&a3)){
        bic[x] <- NA
      }else{
        break
      }
      ####end, change on 051311
    }  
  }
  return(x)
}
  
