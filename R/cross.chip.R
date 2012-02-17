cross.chip <- function(indx1,post,ref,het,hom,ref.a,alt.a,frag,fragpos,
                       alldelta,allsigma,f1,f3,mu1,mu3,dbsnp.info,bad,
                       final.bp,plot=FALSE){

  ## 2. Now only look at the chosen allele, perform EM if common SNPs
  ## There are three categories: common SNPs w/ EM, rare SNPs w/o EM, refs
  d <- cbind(alldelta[indx1,final.bp,,1],alldelta[indx1,4-final.bp,,2])
  ai <- cbind(allsigma[indx1,final.bp,,1],allsigma[indx1,4-final.bp,,2])
  
  f1i <- f1[indx1,,]
  f3i <- f3[indx1,,]
  
  ##rf1i <- as.matrix(f1i)-mu1
  ##rf3i <- as.matrix(f3i)-mu3

  dbsnpID <-ifelse( is.null(dbsnp.info)|nrow(dbsnp.info)==0, numeric(0), which(dbsnp.info$ExonID==frag&dbsnp.info$BasePos==fragpos))
                                        
    #####Begin, change on May 13,2011
  ref.cond <- max(post[,2],na.rm=T)<0.001&max(post[,3],na.rm=T)<0.0001
  hom.var.cond <- max(post[,2],na.rm=T)<0.001&max(post[,1],na.rm=T)<0.0001
  ##browser()
  ##if(max(post[,2], na.rm=T)<0.001&max(post[,3], na.rm=T)<0.0001){ #add na.rm arguments
  if(ref.cond|hom.var.cond){
    ##print("all ref")
    genotypes <- rep(NA,nrow(d))
    genotypes[!bad] <- ifelse(ref.cond,1,3)
    indx <- bic <- NULL
    posterior <- matrix(NA,nrow(d),3)
    if(ref.cond){
      posterior[!bad,1] <- 1
      posterior[!bad,2:3] <- 0
    }
    if(hom.var.cond){
      posterior[!bad,1:2] <- 0
      posterior[!bad,3] <- 1
    }
  #######End, change on May 13, 2011
    fit <- list(posterior=posterior)
    
    if(plot==TRUE){
      plot(d[!bad,2],d[!bad,1],xlim=c(-4,4),ylim=c(-4,4),
           xlab=paste("log(",het,")-",sep=""),
           ylab=paste("log(",het,")+",sep=""))
      abline(v=0,col="grey",lwd=2,lty=2)
      abline(h=0,col="grey",lwd=2,lty=2)
    }
    cat <- "ref only"
  }else{
    tmp2 <- kmeans(na.omit(post[,2]),2,nstart=25)#add na.omit by NZ)
    tmp3 <- kmeans(na.omit(post[,3]),2,nstart=25)
    
    if(max(tmp2$center)/min(tmp2$center)>2&max(post[,2])>0.001){
      heti <- (tmp2$cluster==which(tmp2$center==max(tmp2$center)))&post[,2]>0.001
    }else{
      heti <- rep(FALSE,dim(post)[1])
    }
    
    if(max(tmp3$center)/min(tmp3$center)>2&max(post[,3])>0.0001){
      ##   homi <- post[,3]> min(tmp3$center)*5
      homi <- tmp3$cluster==which(tmp3$center==max(tmp3$center))
    }else{
      homi <- rep(FALSE,dim(post)[1])
    }
    
    if(sum(heti&homi)>0){
      tmp4 <- post[heti&homi,2:3]
      if(which(tmp4==max(tmp4))==1){
        homi[heti&homi]=FALSE
      }else{
        heti[heti&homi] <- FALSE
      }
    }
    
    
    ##indicate <- sum(homi)>2&max(tmp2$center)>0.01&max(tmp2$center)/min(tmp2$center)>8
    refi <- !heti&!homi
    indicate <- sum(refi)>2 &(sum(heti)>4|sum(homi)>2)&max(tmp2$center)>0.01&max(tmp2$center)/min(tmp2$center)>8  #modified on 4/4/11 replaced indicate <- (sum(heti)>4|sum(homi)>2)&max(tmp2$center)>0.01&max(tmp2$center)/min(tmp2$center)>8
    ##browser()
    
    if(indicate){
      ##browser()
      ## check if this is a dbSNP
      if(length(dbsnpID)>0&&dbsnp.info$Alt.Allele[dbsnpID]!="-"){
        stopifnot(as.character(dbsnp.info$Ref[dbsnpID])==ref.a)
      }
      ###
      ## perform EM
      ###
      ##browser()
      em3 <- mclustModel(data=d[!bad,],BICvalues=mclustBIC(d[!bad,]),G=3)
      em2 <- mclustModel(data=d[!bad,],BICvalues=mclustBIC(d[!bad,]),G=2)
      em1 <- mclustModel(data=d[!bad,],BICvalues=mclustBIC(d[!bad,]),G=1)

      allfit <- list(fit1=em3,fit2=em2,NULL)
      bic <- c(min(-em3$bic),min(-em2$bic),min(-em1$bic))

      ####begin, change on 051311
      ##indx <- find.best.fit(bic=bic,allfit=allfit)
      indx <- ifelse(length(dbsnpID)>0,order(bic)[1],
                     find.best.fit(bic=bic,allfit=allfit))
      ####end, change on 051311
      fit <- allfit[[indx]]      
      fit$posterior <- matrix(NA,nrow(d),3) # used to be ifelse(indx==1,3,2)) this gives error when indx==2/5 6--NZ
      if(indx!=3){
        m <- t(fit$parameter$mean)
      }
      
      if(indx<3){
        cat <- "common SNPs w/ EM"
      }else{
        cat <- "ref only"
      }
      genotypes <- rep(NA,nrow(d))
      if(indx==1){
        fit$posterior[!bad,] <- fit$z[,order(m[,1],decreasing=TRUE)]
        genotypes[!bad] <- map(fit$posterior[!bad,])
      }
      if(indx==2){
        ##fit$posterior[!bad,] <- fit$z[,order(m[,1],decreasing=TRUE)]
        order.m <- order(m[,1],decreasing=TRUE)
        m <- m[order.m,]
        if(sqrt(sum(m[1,]^2))>=sqrt(sum(m[2,]^2))){
          if(sqrt(sum(m[2,]^2))>0.5&sum(m[2,]<0)==2){
            indx <- 5 ##genotype 1, 3
          }else{
            indx <- 4 ##genotype 1, 2
          }
        }else{
          if(sqrt(sum(m[1,]^2))>0.5&sum(m[1,]>0)==2){
            indx <- 5
          }else{
            indx <- 6 ##genotype 2, 3
          }
        }       
      }
#################This chunk is from the old version  replaced on 3/12/11    
#      if(indx==3){ fit$posterior[!bad,] <- matrix(rep(c(1,0),sum(!bad)), sum(!bad),2,byrow=T)
#		genotypes[!bad] <- 1      }
#      if(indx==5){  tmp <- map(fit$posterior[!bad,])
#        genotypes[!bad] <- ifelse(tmp==1,1,3)      }
#      if(indx==6){       genotypes[!bad] <- map(fit$posterior[!bad,])+1      }
#   if(indx==4){        genotypes[!bad] <- map(fit$posterior[!bad,])      }
##################################

      if(indx==3){
        if(em1$parameter$mean[1]<0&em1$parameter$mean[2]<0){
          fit$posterior[!bad,] <- matrix(rep(c(0,0,1),sum(!bad)),
                                         sum(!bad),3,byrow=T)
          genotypes[!bad] <- 3
        }else{
          fit$posterior[!bad,] <- matrix(rep(c(1,0,0),sum(!bad)),
                                         sum(!bad),3,byrow=T)
          genotypes[!bad] <- 1
        }
      }
      if(indx==5){
        fit$posterior[!bad,c(1,3)] <- fit$z[,order.m]
        tmp <- map(fit$posterior[!bad,c(1,3)])
        genotypes[!bad] <- ifelse(tmp==1,1,3)
      }
      
      if(indx==6){
        fit$posterior[!bad,c(2,3)] <- fit$z[,order.m]
        genotypes[!bad] <- map(fit$posterior[!bad,c(2,3)])+1
      }
      
      if(indx==4){
        fit$posterior[!bad,c(1,2)] <- fit$z[,order.m]
        genotypes[!bad] <- map(fit$posterior[!bad,c(1,2)])
      }
      ##browser()
      ##Plot the calls
      if(plot==TRUE){
        if(indx!=3){
          mclust2Dplot(data=d[!bad,],parameters=fit$parameter,z=fit$z,xlim=c(-4,4),ylim=c(-4,4))
        }
        if(indx==3){
          plot(d[!bad,1],d[!bad,2],xlim=c(-4,4),ylim=c(-4,4))
        }
        abline(v=0,col="grey",lwd=2,lty=2)
        abline(h=0,col="grey",lwd=2,lty=2)
      }
    }
    if(!indicate){
      if(sum(heti)>0|sum(homi)>0){
        cat <- "rare variants"
        ## check if this is a dbSNP
        if(length(dbsnpID)>0){
          if(dbsnp.info$Alt.Allele[dbsnpID]!="-"){
            stopifnot(as.character(dbsnp.info$Ref[dbsnpID])==ref.a)
          }
          prior <- c(0.9752,0.0247,0.0001)
        }else{
          prior <- c(0.999^2,2*0.001*(1-0.001),0.001^2)
        }
        
        genotypes <- rep(NA,nrow(d))
        centers <- apply(d[!bad,][!heti,,drop=FALSE],2,median)
        #######
        ##05-03-2011
        if(centers[1]<0&centers[2]<0){
          centers <- c(1.7,1.7)
        }	
        ##05-03-2011
        #######			        
        if(sum(!heti&!homi)>sum(heti|homi)){
          cov.points <- !heti&!homi
        }else{
          if(sum(heti)>sum(homi)){
            cov.points <- heti
          }else{
            cov.points <- homi
          }
        }
        covs <- cov(d[!bad,][cov.points,])
        ##browser()
        ######
        ##begin, change on 051311
        posterior <- rare.variants.pp(d=d,
                                      ai=ai,bad=bad,center=centers,
                                      cov=covs,
                                      prior=prior,dbsnpID=dbsnpID)
        ##end, change on 051311
        ######
       
        a <- posterior[,2]
        a1 <- max(a,na.rm=T)
        b <- sort(a[!is.na(a)])
        b1 <- b[length(b)-1]
        if(b1!=0){
          if(a1/b1>10^20){
            if(length(dbsnpID)>0){
              posterior[!is.na(posterior[,2])&posterior[,2]==a1,] <- c(0,1,0)
            }
          }
        }
        fit <- list(posterior=posterior)
        genotypes[!bad] <- apply(posterior[!bad,],1,function(a){
          which(a==max(a))[1]
        })
        
        ##check the location of the heterozygotes
        if(sum(genotypes[!bad]!=1)>0){
          a1 <- which(genotypes[!bad]==1)
          a2 <- which(genotypes[!bad]==2)
          a3 <- which(genotypes[!bad]==3)
          if(length(a1)>0&length(a2)>0&length(a3)==0){
            m1 <- apply(d[!bad,][a1,,drop=FALSE],2,mean)
            m2 <- apply(d[!bad,][a2,,drop=FALSE],2,mean)
            if(sum(m2<m1)<2|sum(m2^2)>sum((m2-m1)^2)|(sum(m1^2)<sum(m2^2)&sum(m1<0)==1)){
              posterior[!bad,][a2,] <- c(1,0,0)
              genotypes[!bad][a2] <- 1
            }
          }
          if(length(a1)==0&length(a2)>0&length(a3)>0){
            m2 <- apply(d[!bad,][a2,,drop=FALSE],2,mean)
            m3 <- apply(d[!bad,][a3,,drop=FALSE],2,mean)
            if(sum(m3<m2)<2|sum(m2^2)>sum((m3-m2)^2)|(sum(m3^2)<sum(m2^2)&sum(m2>0)==1)){
              posterior[!bad,][a2,] <- c(0,0,1)
              genotypes[!bad][a2] <- 3
            }
          }
        }
        
        if(plot==TRUE){
          plot(d[!bad,2],d[!bad,1],xlim=c(-4,4),ylim=c(-4,4),
               xlab=paste("log(",het,")-",sep=""),
               ylab=paste("log(",het,")+",sep=""))
               ##main="Raw data")
          abline(v=0,col="grey",lwd=2,lty=2)
          abline(h=0,col="grey",lwd=2,lty=2,)
          if(sum(genotypes==2,na.rm=T)>0){
            points(d[!bad&genotypes==2,2],d[!bad&genotypes==2,1],col="blue")
          }
          if(sum(genotypes==3,na.rm=T)>0){
            points(d[!bad&genotypes==3,2],d[!bad&genotypes==3,1],col="red")
          }
          center <- rbind(centers[2:1],c(0,0),c(-1.7,-1.7))
          scale <- matrix(c(covs[2,2],covs[2,1],covs[1,2],covs[1,1]),2,2)
          for(w in 1:length(unique(genotypes))){
            lines(ellipse(x=scale,centre=center[unique(genotypes)[w],],level=0.90),
                  col="darkgrey")
          }
          points(center[unique(genotypes),1],center[unique(genotypes),2],col="darkgrey",pch=19,cex=.5)
        }
        indx <- bic <- NULL
      }else{
        cat <- "ref only"
        if(plot==TRUE){
          plot(d[!bad,2],d[!bad,1],xlim=c(-4,4),ylim=c(-4,4),
               xlab=paste("log(",het,")-",sep=""),
               ylab=paste("log(",het,")+",sep=""))
          abline(v=0,col="grey",lwd=2,lty=2)
          abline(h=0,col="grey",lwd=2,lty=2)
        }
        genotypes <- rep(NA,nrow(d))
        genotypes[!bad] <- 1
        indx <- bic <- NULL
        posterior <- matrix(NA,nrow(d),3)
        posterior[!bad,1] <- 1
        posterior[!bad,2:3] <- 0
        fit <- list(posterior=posterior)
      }
    }
  }
  ##when all samples are called as het, change them back to ref
  if(sum(genotypes[!bad]==2)==sum(!bad)){
    genotypes[!bad] <- 1
  }
  if(sum(genotypes[!bad]==2)==sum(!bad)-1&sum(genotypes[!bad]==1)==1){
    genotypes[!bad] <- 1
  }
  
###here if there are nonref calls, confirm with the dbSNPs if possible
  note <- NA
  if(length(dbsnpID)>0&&sum(genotypes[!bad]!=1)>0&&dbsnp.info$Alt.Allele[dbsnpID]!="-"){
    aa <- unlist(strsplit(as.character(dbsnp.info$Alt.Allele[dbsnpID]),split="/"))
    if(!dbsnp.info$rs[dbsnpID]%in%c("rs3219466","rs1008681")){
      if(aa!=alt.a){
        note <- "wrong alt"
      }
    }else{
      stopifnot(alt.a=="A"|alt.a=="T")
    }
  }
  
  ##rs <- check.probe.quality(d[!bad,],genotypes[!bad])
  slscore <- calc.slscore(d[!bad,],genotypes[!bad])
  ##output <- list(d,final.bp,het,fit,genotypes,bic,indx,cat,rs,note,slscore)
  ##names(output) <- c("data","AltBaseIndx","BasePair","fit","genotypes","BIC","Best model",
  ##                   "category","Reliability","note","slscore")
  output <- list(d,final.bp,het,fit,genotypes,bic,indx,cat,slscore)
  names(output) <- c("data","AltBaseIndx","BasePair","fit","genotypes","BIC","Best model",
                     "category","slscore")
  return(output)
}
