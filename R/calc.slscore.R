calc.slscore <- function(d,genotypes){
  a1 <- which(genotypes==1)
  a2 <- which(genotypes==2)
  a3 <- which(genotypes==3)

 ## browser()
  slscore.i <- matrix(NA,length(genotypes),2)
  if(length(a1)==length(genotypes)|length(a2)==length(genotypes)|length(a3)==length(genotypes)){
    for(i in 1:length(genotypes)){
      al <- apply((apply(d[-i,],1,function(a,b){
        abs(a-b)
      },b=d[i,])),1,median)
      bl <- abs(d[i,])
      for(k in 1:2){
        slscore.i[i,k] <- (bl[k]-al[k])/max(bl[k],al[k])
      }
    }
  }else{
    if(length(a1)>0){
      for(i in 1:length(a1)){
        if(length(c(a1[i],a2,a3))==length(genotypes)){
          al <- rep(0,2)
        }else{
          al <- apply(apply(d[-c(a1[i],a2,a3),,drop=FALSE],1,function(a,b){
            abs(a-b)
          },b=d[a1[i],]),1,median)
        }
        bl <- apply(apply(d[-a1,,drop=FALSE],1,function(a,b){
          abs(a-b)
        },b=d[a1[i],]),1,median)
        for(k in 1:2){
          slscore.i[a1[i],k] <- (bl[k]-al[k])/max(bl[k],al[k])
        }
      }
    }
    if(length(a2)>0){
      for(i in 1:length(a2)){
        if(length(c(a2[i],a1,a3))==length(genotypes)){
          al <- rep(0,2)
        }else{
          al <- apply(apply(d[-c(a2[i],a1,a3),,drop=FALSE],1,function(a,b){
            abs(a-b)
          },b=d[a2[i],]),1,median)
        }
        bl <- apply(apply(d[-a2,,drop=FALSE],1,function(a,b){
          abs(a-b)
        },b=d[a2[i],]),1,median)
        for(k in 1:2){
          slscore.i[a2[i],k] <- (bl[k]-al[k])/max(bl[k],al[k])
        }
      }
    }
    if(length(a3)>0){
      for(i in 1:length(a3)){
        if(length(c(a3[i],a1,a2))==length(genotypes)){
          al <- rep(0,2)
        }else{
          al <- apply(apply(d[-c(a3[i],a1,a2),,drop=FALSE],1,function(a,b){
            abs(a-b)
          },b=d[a3[i],]),1,median)
        }
        bl <- apply(apply(d[-a3,,drop=FALSE],1,function(a,b){
          abs(a-b)
        },b=d[a3[i],]),1,median)
        for(k in 1:2){
          slscore.i[a3[i],k] <- (bl[k]-al[k])/max(bl[k],al[k])
        }
      }
    }
  }
  return(slscore.i)
}

