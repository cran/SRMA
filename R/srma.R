srma <- function(snps,allpost,alldelta,allsigma,allpairs,allbp.indx,dbsnp.info,
                  f1,f3,mu1,mu3,refCalls, plot=FALSE, ...){
  nsample<-dim(alldelta)[3]
  cat <- NULL
  refsites <- NULL
  ##rs <- NULL
  calls <- posteriors <- matrix(NA,nrow(snps),nsample)
  ##rscore <- matrix(NA,nrow(snps),2)
  ##notes <- rep(NA,nrow(snps))
  slscore <- array(NA,c(nrow(snps),nsample,2))
    
  for(i in 1:nrow(snps)){
    ##print(i)
    frag <- substr(as.character(unlist(snps$ExonID[i])),1,15)
    fragpos <- snps$BasePos[i]
    indx1 <- which(refCalls$Fragment==frag&refCalls$FragPos==fragpos)
    bad <- is.na(alldelta[indx1,1,,1])
    basepair <- allpairs[indx1,4-allbp.indx[indx1],2]
    ref.a <- substr(basepair,1,1)
    alt.a <- substr(basepair,2,2)
    ref <- paste(ref.a,ref.a,sep="")
    het <- basepair
    hom <- paste(alt.a,alt.a,sep="")
    
    post <- allpost[i,!bad,]
    if(sum(!bad)>5){
      if(plot==TRUE){
        ##postscript(paste("figs/snps/snp",i,frag,"-",fragpos,".ps",sep=""))
        postscript(paste("snp",i,frag,"-",fragpos,".eps",sep=""), ...)
        par(mfrow=c(1,1),mar=c(4,5,5,5),pty="s")
      }
        
      fit <- cross.chip(indx1=indx1,post=post,ref=ref,het=het,hom=hom,
                        ref.a=ref.a,
                        alt.a=alt.a,frag=frag,fragpos=fragpos,
                        alldelta=alldelta,allsigma=allsigma,f1=f1,f3=f3,
                        mu1=mu1,mu3=mu3,
                        dbsnp.info=dbsnp.info,
                        bad=bad,final.bp=allbp.indx[indx1],plot=plot) 
      if(plot==TRUE){
        mtext(paste(frag,"-",fragpos," ",fit$BasePair,sep=""),outer=TRUE,line=-3)
        dev.off()
      }
        
      ##rscore[i,] <- fit$Reliability
      ##notes[i] <- fit$note
      slscore[i,!bad,] <- fit$slscore
      ##  if(max(abs(fit$Reliability))<10&min(abs(fit$Reliability))<3|
      ##     fit$Reliability[1]*fit$Reliability[2]<0|!is.na(fit$note))
      ##  if(max(abs(fit$Reliability))<10&(min(abs(fit$Reliability))<3|fit$Reliability[1]*fit$Reliability[2]<0)&sum(fit$genotypes==3,na.rm=T)<4){
      ##    rs <- c(rs,FALSE)
      ##  }else{
      ##    rs <- c(rs,TRUE)
      ##  }
    }else{
      genotypes <- rep(NA,length(bad))
      genotypes[!bad] <- 1
      posterior <- matrix(NA,length(bad),3)
      posterior[!bad] <- 1
      ##rs <- c(rs,FALSE)
      fit <- list(genotypes=genotypes,fit=list(posterior=posterior))
    }
    refsites <- c(refsites,sum(fit$genotypes[!bad]==1)==sum(!bad))
    cat <- c(cat,fit$category) 
    calls[i,!bad] <- c(ref,het,hom)[fit$genotypes[!bad]]
    if(sum(!bad)==1) posteriors[i,!bad] <- max(fit$fit$posterior[!bad,]) else { #add this to avoid error when there is only one good sample
      posteriors[i,!bad] <- apply(fit$fit$posterior[!bad,],1,max)
    }
  }
  ##output <- list(calls,posteriors,cat,refsites,rs,rscore,notes,slscore)
  ##names(output) <- c("calls","posteriors","category","refsites",
  ##                   "reliability","rscore",
  ##                   "notes","slscore")
  output <- list(calls,posteriors,cat,refsites,slscore)
  names(output) <- c("calls","posteriors","category","refsites",
                     "slscore")
  return(output)
}
