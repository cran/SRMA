#######
## find.badamps.R
## Wenyi Wang
#######
###
## criterion 1: [R>=0.9]
###

find.badamps <- function(allqm=allqm,mode.amplevel="average"){
  badamp <- matrix(0,dim(allqm)[1],dim(allqm)[4])

  for(i in 1:dim(allqm)[1]){
    if(mode.amplevel=="average"){
      r <- (allqm[i,3,1,]+allqm[i,3,2,])/2
    }
    if(mode.amplevel=="maximum"){
      r <- mapply(max,allqm[i,3,1,],allqm[i,3,2,])
    }
##    d <- (allqm[i,4,1,]+allqm[i,4,2,])/2
##    t <- (allqm[i,6,1,]+allqm[i,6,2,])/2
    badamp[i,] <- ifelse(r>=0.9,0,1)
  }
  ##overall failure rate:
  ##sum(badamp)/(nrow(badamp)*ncol(badamp))
  rownames(badamp) <- dimnames(allqm)[[1]]
  colnames(badamp) <- dimnames(allqm)[[4]]
  return(badamp)
}

