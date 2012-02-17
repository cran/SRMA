single.array.lm.all.alleles<-function(reorg.Data, verbose=FALSE){
  alldelta.qn<-reorg.Data$alldelta.qn
  allsigma.qn<-reorg.Data$allsigma.qn
  fragl<-rep(reorg.Data$exonlength, reorg.Data$exonlength)
  allpairs.qn<-reorg.Data$allpairs.qn
  allgc<-rep(reorg.Data$gc, reorg.Data$exonlength)
  sample.names<-dimnames(alldelta.qn)[[3]]
  stopifnot(dim(alldelta.qn)[1]==length(allgc))
  n.samp<-dim(alldelta.qn)[3]
  mu1 <- mu3 <- rep(NA,n.samp)
  f1 <- f3 <- errors <- array(NA,dim(alldelta.qn))
  for(i in 1:n.samp){
    if(verbose) cat('Processing Sample', i, ':', sample.names[i], '\n')
    tmp <- !is.na(alldelta.qn[,1,i,1])
    output <- correct.lr(i=i,alldelta.qn=alldelta.qn,allsigma.qn=allsigma.qn,allpairs.qn=allpairs.qn,allgc=allgc,
                         verbose=verbose)
    f1[tmp,,i,] <- array(output$f1,c(sum(tmp),3,2))
    f3[tmp,,i,] <- array(output$f3,c(sum(tmp),3,2))
    errors[tmp,,i,] <- array(output$error,c(sum(tmp),3,2))
    mu1[i] <- output$mu1
    mu3[i] <- output$mu3
    gc()
  }
  chipadj.r1.qn <- list(f1,f3,mu1,mu3,errors)
  names(chipadj.r1.qn) <- c("f1","f3","mu1","mu3","errors")
  return(chipadj.r1.qn)
}
