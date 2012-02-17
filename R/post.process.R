post.process<-function(srma.output, dbsnp.info=NULL, sls.cut=0.67, ampmap=NULL){
  alloutput<-srma.output$output
  ampnames<-names(alloutput)
  refCalls<-srma.output$refCalls
  sample.names<-srma.output$sample.names
  n.sample<-length(sample.names)
  n.base<-nrow(refCalls)
  cat('Organizing the output...\n')
  
  srmacalls.m <- srmaqs.m <- matrix(NA,n.base,n.sample)
  ##srmars.m <- array(NA,c(n.base,2))
  ##srmanotes <- rep(NA,n.base)
  srmaslscore <- array(NA,c(n.base,n.sample,2))
  ##calculate bad quality matrix
  qc.m <- matrix(NA,nrow(refCalls),n.sample)
  for(i in 1:length(ampnames)){
    exoni<- alloutput[[i]]
    x1 <- ifelse(is.na(exoni$calls),FALSE,TRUE)
    qc.m[refCalls$Fragment==ampnames[i],] <- x1
  }
  cat('Checking footprint...\n')
  n.footprint <- 0
  m.footprint <- 0
  n.snps <- 0
  for(i in 1:length(ampnames)){
    cat('Checking', ampnames[i], '... \n')
    exoni <- alloutput[[i]]
    snpsi <- which(!exoni$refsites)
    n.snps <- n.snps+length(snpsi)
    if(length(snpsi)>1){
      snpsi1 <- check.footprint(snpsi=snpsi,exoni=exoni,exonname=ampnames[i], refCalls=refCalls,dbsnp.info=dbsnp.info)
      if(length(snpsi1)<length(snpsi)){
        footprint.indx <- snpsi[!snpsi%in%snpsi1]
        m.footprint <- m.footprint+length(footprint.indx)
        ref <- toupper(as.character(unlist(refCalls$Ref[refCalls$Fragment==ampnames[i]][footprint.indx])))
        n.footprint <- n.footprint+sum(apply(cbind(ref,exoni$calls[footprint.indx,]),1,function(a){
          ref.allele <- a[1]
          ref <- paste(ref.allele,ref.allele,sep="")
          call <- a[2:(n.sample+1)]
          return(sum(!is.na(call)&call!=ref))
        }))
        ##change het calls to ref calls in the footprint positions
        ref.fp <- toupper(as.character(refCalls$Ref[refCalls$Fragment==ampnames[i]][footprint.indx]))
        missing <- is.na(exoni$calls[footprint.indx,])
        exoni$calls[footprint.indx,] <- matrix(rep(paste(ref.fp,ref.fp,sep=""),n.sample),
                                               length(footprint.indx),n.sample)
        exoni$calls[footprint.indx,][missing] <- NA
        ## change posteriors to .99 if they are <.99
        pp <- exoni$posteriors[footprint.indx,]
        exoni$posterior[footprint.indx,][!missing&pp<0.99] <- 0.99
      }
    }
    ##make matrices with footprint corrected
    srmacalls.m[refCalls$Fragment==ampnames[i],] <- exoni$calls
    srmaqs.m[refCalls$Fragment==ampnames[i],] <- exoni$posteriors
    ##srmars.m[refCalls$Fragment==ampnames[i],] <- exoni$rscore
    ##srmanotes[refCalls$Fragment==ampnames[i]] <- exoni$notes
    srmaslscore[refCalls$Fragment==ampnames[i],,] <- exoni$slscore
  }
  

  cat('Identifying SNPs...\n')
  temp.ref<-toupper(as.character(refCalls$Ref))
  temp.ref.m<-matrix(paste(temp.ref, temp.ref, sep=""), nrow=nrow(refCalls), ncol=ncol(srmacalls.m))
  snps.indx <- apply(srmacalls.m!=temp.ref.m, 1, any, na.rm=T) 
  allsnps <- data.frame(refCalls[snps.indx,],srmacalls.m[snps.indx,])
  colnames(allsnps)[4:ncol(allsnps)] <- sample.names
  exons <- as.character(unique(allsnps[,1]))
  
  cat('Identifying insertion/deletions...\n')
  if(is.null(ampmap)) {
    exon.leng<-table(refCalls[,1])
    temp.map<-data.frame(names(exon.leng), as.vector(exon.leng), 1, as.vector(exon.leng)) 
    arrayampmap<-cbind(temp.map, temp.map)
    colnames(arrayampmap)<-c("exon", "exon_len",   "exon_start", "exon_end",   "amplicon"  , "amp_len"  ,  "amp_start" , "amp_end")                                  
  } else {arrayampmap<-ampmap }
  indelsnps <- NULL
  allsnps2 <- NULL
  for(i in 1:length(exons)){
    ampmap <- arrayampmap[which(regexpr(exons[i],arrayampmap$amplicon)>0),,drop=FALSE]
    if(nrow(ampmap)==0) next
    for(k in 1:nrow(ampmap)){
      if(nchar(as.character(ampmap[k,'amplicon']))==17){
        subsnps <- allsnps[allsnps[,1]==exons[i]&allsnps[,2]>=ampmap$'amp_start'[k]&allsnps[,2]<=ampmap$'amp_end'[k],,drop=FALSE]
      }else{
        subsnps <- allsnps[allsnps[,1]==exons[i],,drop=FALSE]
      }
      out <- find.indels(subsnps=subsnps,code=1)
      indelsnps <- rbind(indelsnps,out$indelsnps) 
      if(k>1){
        x1 <- paste(out$singlesnps$Fragment,out$singlesnps$FragPos)%in%paste(indelsnps$Fragment,
                                                                             indelsnps$FragPos)
        out$singlesnps <- out$singlesnps[!x1,]
      }
      allsnps2 <- rbind(allsnps2,out$singlesnps)
    }
    id1 <- paste(allsnps2[,1],"-",allsnps2[,2],sep="") #335 unique
    id2 <- paste(indelsnps[,1],"-",indelsnps[,2],sep="") #296 unique
    
    allsnps2 <- allsnps2[match(unique(id1),id1),]
    indelsnps <- indelsnps[match(unique(id2),id2),]
  }
  
  allsnps2 <- allsnps2[!is.na(allsnps2[,1]),]
  
  cat('Applying Silhouette score cutoff... \n')
  cutoff<-sls.cut
  
  slscore <- t(sapply(1:n.base,function(a,srmaslscore){
    apply(srmaslscore[a,,],1,min)
  },srmaslscore=srmaslscore))
  
  slscore2 <- t(sapply(1:n.base,function(a,srmaslscore){
    apply(srmaslscore[a,,],2,median,na.rm=T)
  },srmaslscore=srmaslscore))
  
  tmp.finalcall <- t(apply(allsnps2,1,function(a,slscore,refCalls,cutoff,dbsnp.info){
    frag <- as.character(unlist(a[1]))
    fragpos <- as.numeric(as.character(unlist(a[2])))
    ref.allele <- toupper(as.character(unlist(a[3])))
    ref <- paste(ref.allele,ref.allele,sep="")
    score <- slscore[refCalls$Fragment==frag&refCalls$FragPos==fragpos,]
    sc <- median(score[!is.na(score)])
    call <- as.character(unlist(a[4:length(a)]))
    if(sc<=cutoff){
      call[!is.na(call)] <- "X"
    }else{
      if(sum(dbsnp.info[,'ExonID']%in%frag&dbsnp.info[,'BasePos']%in%fragpos)==0){
        call[!is.na(call)&call!=ref&!is.na(score)&score<=cutoff] <- "N"
      }
    }
    if(sum(!is.na(call)&call!="N"&call!="X"&call!=ref)==0){
      goodsnps <- "Bad"
    }else{
      goodsnps <- "Good"
    }
    return(c(goodsnps,call))
  },slscore=slscore,refCalls=refCalls,cutoff=cutoff,dbsnp.info))
  
  
  allsnps3 <- allsnps2
  allsnps3[,4:ncol(allsnps3)] <- tmp.finalcall[,-1]
  tmp2 <- tmp.finalcall[,1]=="Good" 
  allsnps3 <- allsnps3[tmp2,] 
 
  ##srma.metric<-list(srma.qc=qc.m, srma.calls=srmacalls.m, srma.qs=srmaqs.m, srma.rs=srmars.m, srma.slscore=srmaslscore)
  if(is.null(indelsnps))indelsnps<-data.frame(matrix(numeric(0), nrow=0, ncol=3+n.sample))
  if(is.null(allsnps3)) allsnps3<-data.frame(matrix(numeric(0), nrow=0, ncol=3+n.sample))
  colnames(indelsnps)[1:3]<-colnames(allsnps3)[1:3]<-c( 'Fragment',  'FragPos',  'Ref')   

  srma.metric<-list(srma.qc=qc.m, srma.calls=srmacalls.m, srma.qs=srmaqs.m, srma.slscore=srmaslscore)
  footprint.par<-list(base.n=n.footprint, sample.n=m.footprint)
  final.output<-new('SRMAResults', footprint=footprint.par, srmaMat=srma.metric, srmaIndelsnps=indelsnps, srmaSnps=allsnps3, sampleNames=sample.names, slscoreCutoff=sls.cut)
  return(final.output)
}


