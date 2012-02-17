#####
## reorganize delta, sigma and pairs into an array format
#####
reorganize.datafiles <- function(delta,sigma,pairs,n.bases,n.samples,sample.names,
                                 exon.names,ref.info){
  
  alldelta.qn <- allsigma.qn <- array(NA,c(n.bases,3,n.samples,2))
  dimnames(alldelta.qn) <- dimnames(allsigma.qn) <- list(1:n.bases,c("bp1","bp2","bp3"),sample.names,
                                                         c("sense","antisense"))
  allpairs.qn <- array(NA,c(n.bases,3,2))
  dimnames(allpairs.qn) <- list(1:n.bases,c("bp1","bp2","bp3"),
                                c("sense","antisense"))

  unitnames <- names(delta)
  start <- 1
  new.ref.info <- NULL
  gc <- exonlength <- NULL
  for(i in exon.names){
    print(i)
    indx <- which(regexpr(i,unitnames)>0)
    y <- ref.info[ref.info$Fragment==i,c("Fragment","FragPos","Ref")]
    new.ref.info <- rbind(new.ref.info,y)
    exonlength <- c(exonlength,nrow(y))
    gc <- c(gc,sum(y$Ref%in%c("c","g"))/nrow(y))
    
##    end <- start+nrow(delta[[indx[1]]])-1
    end <- start+nrow(y)-1
    for(j in 1:2){
      alldelta.qn[start:end,,,j] <- delta[[indx[j]]]
      allsigma.qn[start:end,,,j] <- sigma[[indx[j]]]
      allpairs.qn[start:end,,j] <- pairs[[indx[j]]]
    }
    start <- end+1
  }

  output <- list(alldelta.qn,allsigma.qn,allpairs.qn,new.ref.info,gc,exonlength)
  names(output) <- c("alldelta.qn","allsigma.qn","allpairs.qn","ref.info","gc","exonlength")
  return(output)
}

