check.footprint <- function(snpsi,exoni,exonname,refCalls,dbsnp.info){
  ##define 2 functions
  num.homs <- function(x,calls,ref){
    call <- calls[x,]
    all.alleles <- unique(unlist(strsplit(call,split="")))
    alt.allele <- all.alleles[!is.na(all.alleles)&all.alleles!=ref[x]]
    return(call==paste(alt.allele,alt.allele,sep=""))
  }
  num.hets <- function(x,calls,ref){
    call <- calls[x,]
    all.alleles <- unique(unlist(strsplit(call,split="")))
    alt.allele <- all.alleles[!is.na(all.alleles)&all.alleles!=ref[x]]
    return(call==paste(ref[x],alt.allele,sep=""))
  }
  calls <- exoni$calls[snpsi,,drop=FALSE]
  refcallsi<-refCalls[refCalls$Fragment==exonname,]
  ref <- toupper(as.character(refcallsi[snpsi,'Ref']))
  goodsnp <- rep(TRUE,length(ref))
  
  hom <- t(sapply(1:length(ref),num.homs,calls=calls,ref=ref))
  het <- t(sapply(1:length(ref),num.hets,calls=calls,ref=ref))
  
  id <- 1
  while(id<=length(snpsi)){
    nb <- which(snpsi<=snpsi[id]+12&snpsi>=snpsi[id]-12&snpsi!=snpsi[id])
    if(length(nb)>0){
      hom.sample <- hom[id,]
      for(i in nb){
        ##scenario 1, hom decreases delta in the neighborhood
        tmp1 <- sum(het[i,!hom.sample],na.rm=T)<3
        ##tmp1 <- sum(het[i,!hom.sample],na.rm=T)==0
        tmp2 <- sum(hom[i,hom.sample],na.rm=T)<sum(hom.sample,na.rm=T)
        tmp3 <- sum(het[i,hom.sample],na.rm=T)>0
        
        ##scenario 2, hom increases delta in the neighborhood
        tmp4 <- sum(het[i,!hom.sample],na.rm=T)>(sum(!hom.sample,na.rm=T)-3)
        tmp5 <- sum(het[i,hom.sample],na.rm=T)<3
        
        ##if(tmp1&&tmp2&&tmp3)
        if((tmp1&&tmp2&&tmp3)||(tmp4&&tmp2&&tmp5)){
          goodsnp[i] <- FALSE
        }
      }
    }
    id <- id+1
  }
  ##add back reported dbSNPs
  ##dbsnps <- dbsnp.info[dbsnp.info$ExonID==exonname,]
  ##if(nrow(dbsnps)>0)   goodsnp[which(refcallsi$FragPos %in% dbsnps$BasePos)] <- TRUE
  ##changed on 10/13/11
  ####
  if(!is.null(dbsnp.info) && nrow(dbsnp.info)!=0 && sum(dbsnp.info$ExonID==exonname)>0) {
    dbsnps <- dbsnp.info[dbsnp.info$ExonID==exonname,]
    goodsnp[which(refcallsi$FragPos %in% dbsnps$BasePos)] <- TRUE
  }
  return(snpsi[goodsnp])
}  
