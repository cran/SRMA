call.genotype<-function(reorg.Dat, singlechip.pp, dbsnp.info=NULL){
  refCalls<-reorg.Dat$ref.info
  ##proofing dbsnp.info dataset
  if(!is.null(dbsnp.info)){
    if(!all(c("rs","Ref","Alt.Allele","ExonID","BasePos") %in% colnames(dbsnp.info) )) {
      stop('dbSNP data has to contain columns rs, Ref, Alt.Allele, ExonID, BasePos')} else {
        keep.indx1 <- which(paste(dbsnp.info$ExonID, dbsnp.info$BasePos) %in% paste(refCalls$Fragment, refCalls$FragPos))
        dbsnp.info<-ifelse(length(keep.indx1)==0, NULL, dbsnp.info[keep.indx1,])
      }
  }
  ##initialize parameters
  ampnames<-as.vector(unique(reorg.Dat$ref.info$Fragment))
  alldelta.qn<-reorg.Dat$alldelta.qn
  allsigma.qn<-reorg.Dat$allsigma.qn
  allpairs.qn<-reorg.Dat$allpairs.qn
  sample.names<-dimnames(alldelta.qn)[[3]]
  nsample<-length(sample.names)
  ###
  allbp.indx <- singlechip.pp$allbp.indx
  post.prob<-singlechip.pp$post.prob
  f1<-singlechip.pp$chipadj.qn$f1
  f3<-singlechip.pp$chipadj.qn$f3
  mu1<-singlechip.pp$chipadj.qn$mu1
  mu3<-singlechip.pp$chipadj.qn$mu3
  alloutput <- list()
  ##j is the number of run for paralelle.
  for(i in 1:length(ampnames)){
    cat('Working on amplicon',i, ':', ampnames[i], as.character(Sys.time()), '\n')
    snps <- refCalls[refCalls$Fragment==ampnames[i],]
    names(snps)[1:2] <- c("ExonID","BasePos")
    allpost<-post.prob[[i]]
    output <- srma(snps=snps,allpost=allpost,alldelta=alldelta.qn,
                   allsigma=allsigma.qn,
                   allpairs=allpairs.qn,
                   allbp.indx=allbp.indx,dbsnp.info=dbsnp.info,
                   f1=f1,f3=f3,mu1=mu1,mu3=mu3,refCalls=refCalls,plot=FALSE)
    alloutput <- c(alloutput,list(output))
  }
  names(alloutput)<-ampnames
  srma.output<-list(output=alloutput, sample.names=sample.names, refCalls=refCalls)
  return(srma.output)
}
