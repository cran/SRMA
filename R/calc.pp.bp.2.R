calc.pp.bp.2 <- function(snps,refCalls,alldelta,allpairs,f1,f3,errors1,errors3,
                         allbp.indx,...){

  nsample<-dim(alldelta)[3]      
  all.post <- array(NA,c(nrow(snps),nsample,3)) #snpID,sampleID,basepair,genotype
  all.geno <- matrix(NA,nrow(snps),nsample)
  
  for(i in 1:nrow(snps)){
    frag <- substr(as.character(unlist(snps$ExonID[i])),1,15)
    fragpos <- snps$BasePos[i]
    indx1 <- which(refCalls$Fragment==frag&refCalls$FragPos==fragpos)
    if(length(indx1)>0){
      bad <- is.na(alldelta[indx1,1,,1])
      basepair <- allpairs[indx1,4-allbp.indx[indx1],2]
      ref.a <- substr(basepair,1,1)
      alt.a <- substr(basepair,2,2)
      ref <- paste(ref.a,ref.a,sep="")
      het <- basepair
      hom <- paste(alt.a,alt.a,sep="")
      output <- position.post(indx1=indx1,ref=ref,het=het,hom=hom,
                              bp.indx=allbp.indx[indx1],
                              alldelta=alldelta,f1=f1,f3=f3,
                              errors1=errors1,errors3=errors3,round=2)
      all.post[i,!bad,] <- output$post[!bad,]
      all.geno[i,!bad] <- output$genotypes[!bad]
    }
  }
    
  output <- list(all.post,all.geno)
  names(output) <- c("all.post","all.geno")
  return(output)
}
  
