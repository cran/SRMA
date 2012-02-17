get.delta.sigma.qn<-function(qc.out, ref.info){
  if(!all(c("Fragment","FragPos","Ref") %in% colnames(ref.info))) stop('ref.info need to contain columns with fragment name (Fragment), position (FragPos) and reference allele (Ref) information!')
  chipType<-qc.out$par$chipType
  dataDir<-qc.out$par$dataDir
  chipOrder<-qc.out$par$chipOrder
  ampmap<-qc.out$par$ampmap
  badamp<-qc.out$badamp
  celint.qn<-get.intensity(chipType=chipType,
                           dataDir=dataDir,
                           chipOrder=chipOrder,
                           ampmap=ampmap,
                           badamp=badamp,
                           chip.normalize=T,
                           amp.normalize=T
                           )
# This take about 7 hours   
  All.qn<-get.delta.sigma(chipType=chipType,
                          dataDir=dataDir,
                          chipOrder=chipOrder,
                          celint=celint.qn
                          )
  unitnames <- names(All.qn$all.delta.qn) 
  exonnames <- unique(substr(unitnames,1,nchar(unitnames)-5)) 
  exonnames<-exonnames[regexpr("AFFX",exonnames)>0|regexpr("ENSE",exonnames)>0]
  n.bases <- 0
  for(i in 1:length(exonnames)){
    unitindx <- which(regexpr(exonnames[i],unitnames)>0)[1]
    n.bases <- n.bases+nrow(All.qn$all.delta.qn[[unitindx]])
  }
  reorg.data<-reorganize.datafiles(delta=All.qn$all.delta.qn,
                                   sigma=All.qn$all.sigma.qn,
                                   pairs=All.qn$all.pairs.qn,
                                   n.bases=n.bases,
                                   n.samples=length(All.qn$celnames),
                                   sample.names=All.qn$celnames,
                                   exon.names=exonnames,
                                   ref.info=ref.info)
  return(reorg.data)
}
      
