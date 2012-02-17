###################
#This function calls position.post function in helpfunctionsv2.R ---NZ Jan 18
#position.post calls calc.post, calc.post.2, calc.ll
#
######################
calc.pp.bp <- function(snps,refCalls,alldelta.qn,allpairs.qn,f1,f3,errors,verbose=F, ...){
  stopifnot(all(c(dim(alldelta.qn)[1], dim(alldelta.qn)[1], dim(refCalls)[1])==nrow(snps)))
  ##added to ensure the dimentions
  nsample<-dim(alldelta.qn)[3] ### added to replace 40 
  ## nrow(snps) are used to replace site number

  all.post <- array(NA,c(nrow(snps),nsample,3,3)) #snpID,sampleID,basepair,genotype
  all.geno <- array(NA,c(nrow(snps),nsample,3))
  all.bp.indx <- rep(NA,nrow(snps))
  
  for(i in 1:nrow(snps)){
    if(verbose)  cat('Calculating posterior probability for Base', i, '\n')
    frag <- substr(as.character(unlist(snps$ExonID[i])),1,15)
    fragpos <- snps$BasePos[i]
    indx1 <- which(refCalls$Fragment==frag&refCalls$FragPos==fragpos)
    if(length(indx1)>0){
      bad <- is.na(alldelta.qn[indx1,1,,1])
      het <- hom <- NULL
      for(j in 1:3){
        basepair <- allpairs.qn[indx1,4-j,2]
        ref.a <- substr(basepair,1,1)
        alt.a <- substr(basepair,2,2)
        ref <- paste(ref.a,ref.a,sep="")
        het <- c(het,basepair)
        hom <- c(hom,paste(alt.a,alt.a,sep=""))
        output <- position.post(indx1=indx1,ref=ref,het=het[j],hom=hom[j],bp.indx=j,
                                alldelta=alldelta.qn,f1=f1,f3=f3,
                                errors1=errors,errors3=errors,round=1)
        all.post[i,!bad,j,] <- output$post[!bad,]
        all.geno[i,!bad,j] <- output$genotypes[!bad]
      }
    
      post <- all.post[i,!bad,,,drop=FALSE]
      a <- matrix(NA,3,3)
      for(j in 1:3){
        for(k in 1:3){
          a[j,k] <- sum(post[1,,j,k])
        }
      }
      b <- which(a[,3]==max(a[,3]))
      all.bp.indx[i] <- b[1]
    }
  }
  
  output <- list(all.post,all.geno,all.bp.indx)
  names(output) <- c("all.post","all.geno","all.bp.indx")
  return(output)
}

