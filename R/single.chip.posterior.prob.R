single.chip.posterior.prob <- function(reorg.Data, verbose=FALSE){
  alldelta.qn<-reorg.Data$alldelta.qn
  allsigma.qn<-reorg.Data$allsigma.qn
  fragl<-rep(reorg.Data$exonlength, reorg.Data$exonlength)
  allpairs.qn<-reorg.Data$allpairs.qn
  allgc<-rep(reorg.Data$gc, reorg.Data$exonlength)
  sample.names<-dimnames(alldelta.qn)[[3]]
  stopifnot(dim(alldelta.qn)[1]==length(allgc))
  n.samp<-dim(alldelta.qn)[3]
  
  cat('Correcting logratios for all alleles in single array...\n')
  chipadj.r1.qn <- single.array.lm.all.alleles(reorg.Data, verbose=verbose)
  ##calculate posterior probability for every base and every position then choose principal alternative allele
  cat('calculating posterior probabilities for all bases and choosing principal alternative allele...\n')
  ref.info<-reorg.Data$ref.info
  snps <- ref.info
  names(snps)[1:2] <- c("ExonID","BasePos")
  bpfrompp1.qn <- calc.pp.bp(snps,ref.info,alldelta.qn,allpairs.qn,f1=chipadj.r1.qn$f1,f3=chipadj.r1.qn$f3,
                             errors=chipadj.r1.qn$errors,
                             verbose=verbose)
  ##do second round logratio adjustment
  cat('adjusting logratios for selected principal alternative alleles... \n')
  fragl.mat<-cbind(fragl, fragl)
  allgc.mat<-cbind(allgc, allgc)
  allbp.indx <- bpfrompp1.qn$all.bp.indx
  allpost <- bpfrompp1.qn$all.post
  ##sample.names<-dimnames(alldelta.qn)[[3]]
  nsample<-length(sample.names)
  nsite<-nrow(alldelta.qn)
  f1 <- f3 <- errors1 <- errors3 <- array(NA,dim=c(nsite,nsample,2))
  mu1 <- mu3 <- NULL
  for(i in 1:nsample){
    if(verbose) cat('Adjusting Sample',i,':', sample.names[i] ,'\n')
    good <- !is.na(alldelta.qn[,1,i,1])
    output <- correct.lr.singlebp(i=i,alldelta=alldelta.qn,allsigma=allsigma.qn,
                                  fragl=fragl.mat,allpairs=allpairs.qn,
                                  allgc=allgc.mat,allbp.indx=allbp.indx,allpost=allpost,
                                  verbose=F) #add verbose argument
    f1[good,i,] <- matrix(output$f1,sum(good),2)
    f3[good,i,] <- matrix(output$f3,sum(good),2)
    errors1[good,i,] <- matrix(output$error1,sum(good),2)
    errors3[good,i,] <- matrix(output$error3,sum(good),2)

    mu1 <- rbind(mu1,apply(f1[,i,],2,median,na.rm=T))
    mu3 <- rbind(mu3,apply(f3[,i,],2,median,na.rm=T))
  }
  chipadj.r2.qn <- list(f1,f3,mu1,mu3,errors1,errors3)
  names(chipadj.r2.qn) <- c("f1","f3","mu1","mu3","errors1","errors3")
  ##calculate single base posterior probability
  cat('Calculating single-base posterior probability...\n')
  ampnames<-unlist(unique(ref.info$Fragment))
  post.prob<-list()
  for(i in 1:length(ampnames)){
    if(verbose) print(paste(i, ':', ampnames[i]))
    snps <- ref.info[ref.info$Fragment==ampnames[i],]
    names(snps)[1:2] <- c("ExonID","BasePos")
    pp2 <- calc.pp.bp.2(snps=snps,refCalls=ref.info,alldelta=alldelta.qn,allpairs=allpairs.qn,
                      f1=f1,f3=f3,errors1=errors1,errors3=errors3,
                       allbp.indx=allbp.indx)
    allpost<-pp2$all.post
    post.prob<- c(post.prob,list(allpost))
  }
  sc.posterior.prob<-list(post.prob=post.prob, allbp.indx=allbp.indx, chipadj.qn=chipadj.r2.qn)
  return(sc.posterior.prob)
}
