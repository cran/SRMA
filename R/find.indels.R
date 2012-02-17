find.indels <- function(subsnps,code){
n.sample<-ncol(subsnps)-3
  if(nrow(subsnps)<=2){
    indelsnps <- NULL
    singlesnps <- subsnps
  }else{
    ref.allele <- toupper(as.character(subsnps$Ref))
    ref <- paste(ref.allele,ref.allele,sep="")
    cl.matrix <- matrix(NA,length(ref.allele),n.sample)
    for(i in 1:nrow(subsnps)){
      call <- subsnps[i,4:ncol(subsnps)]
      if(code==1){
        snpid <- which(call!=ref[i])
        refid <- which(call==ref[i])
      }
      if(code==2){
        snpid <- which(call!=ref.allele[i])
        refid <- which(call==ref.allele[i])
      }
      cl <- rep(NA,n.sample)
      if(code==1){
        cl[snpid] <- sapply(as.character(unlist(call[snpid])),function(a,ref.allele){
          b <- unlist(strsplit(a,split=""))
          if(sum(b%in%ref.allele)>0){
            return(1)
          }else{
            return(2)
          }
        },ref.allele=ref.allele[i])
        cl[refid] <- 0
      }
      if(code==2){
        cl[snpid] <- sapply(as.character(unlist(call[snpid])),function(a,ref.allele){
          if(sum(a%in%c("A","C","T","G"))>0){
            return(2)
          }else{
            return(1)
          }
        },ref.allele=ref.allele[i])
        cl[refid] <- 0
      }
      cl.matrix[i,] <- unlist(cl)
      ##print(i)
      ##print(cl.matrix)
    }
    homs <- apply(cl.matrix,1,function(a){
      sum(!is.na(a)&a==2)>0
    })
    hets <- apply(cl.matrix,1,function(a){
      sum(!is.na(a)&a==1)>3&sum(!is.na(a)&a==1)<sum(!is.na(a))-1
    })
    goodsnpsid <- homs|hets
    
    tmp1 <- cl.matrix[!goodsnpsid,,drop=FALSE]
    if(nrow(tmp1)<=2){
      indelsnps <- NULL
      singlesnps <- subsnps
    }else{
      tmp2 <- apply(tmp1,2,sum,na.rm=T)
      if(sum(tmp2>2)>0){
        indelsnps <- NULL
        sampleid <- which(tmp2>2)
        singlesnps <- subsnps
        rm.snps <- rep(FALSE,nrow(tmp1))
        for(k in 1:length(sampleid)){
          x <- tmp1[,sampleid[k]]
          string <- which(x==1)          
          ##indelsnps <- rbind(indelsnps,
          ##                   subsnps[!goodsnpsid,][string[1]:string[length(string)],])
          indelsnps <- rbind(indelsnps,subsnps[!goodsnpsid,][string,])
          ## tmp1[string[1]:string[length(string)],sampleid[k]] <- 99
          tmp1[string,sampleid[k]] <- 99
          ##singlesnps[!goodsnpsid,4:ncol(subsnps)][string[1]:string[length(string)],
          ##                             sampleid[k]] <- NA
          singlesnps[!goodsnpsid,4:ncol(subsnps)][string,sampleid[k]] <- NA
        }
        rm.snps <- apply(tmp1,1,function(a){
          sum(!is.na(a)&a==99)>0
        })
        singlesnps <- singlesnps[-which(!goodsnpsid)[rm.snps],]
      }else{
        indelsnps <- NULL
        singlesnps <- subsnps
      }
    }
  }

  id1 <- paste(singlesnps[,1],"-",singlesnps[,2],sep="") #335 unique
  id2 <- paste(indelsnps[,1],"-",indelsnps[,2],sep="") #296 unique
  
  singlesnps <- singlesnps[match(unique(id1),id1),]
  indelsnps <- indelsnps[match(unique(id2),id2),]
  
  out <- list(indelsnps,singlesnps)
  names(out) <- c("indelsnps","singlesnps")
  return(out)
}
