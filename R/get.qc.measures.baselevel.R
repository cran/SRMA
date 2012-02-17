######
## get.qc.measures.baselevel.R
## Wenyi Wang
## prepare the measures used for quality assessment and normalization
#####

get.qc.measures.baselevel <- function(celint=celint,
                                      chipType=chipType,
                                      exonList=NA,
                                      celnames=celnames,
                                      mode.baselevel="refcallrate"){

  all.logPM <- all.refc <- all.dm <-
    all.am <- all.rle <- vector("list",length(names(celint)))

  unitnames <- names(celint)

  exons <- c("AFFX-TagIQ-EX",unique(substr(unitnames[regexpr("ENSE",unitnames)>0],1,15)))
  if(sum(is.na(exonList))==0){
    subexons <- exons[exonList]
  }else{
    subexons <- exons
  }

  cdf <- AffymetrixCdfFile$byChipType(chipType,tags="byExon");
  if(chipType=="MitoP-2r520651"){
    units <- 10:nbrOfUnits(cdf);
    cellQuartets <- getCellQuartets(cdf,units=units);
  }else{
    cellQuartets <- getCellQuartets(cdf);
  }

  for(i in 1:length(subexons)){
    ##print(i)
    ## read in probe intensities
    unitindx <- which(regexpr(subexons[i],unitnames)>0)
    for(j in 1:2){
      cells <- t(cellQuartets[[unitindx[j]]][[1]])
      x <- celint[[unitindx[j]]]
      dim(x) <- c(nrow(x)/4,4,length(celnames))
      dimnames(x) <- list(rownames(cells),c("T","G","C","A"),celnames)

      if(mode.baselevel=="refcallrate"){
        if(j==1){
          x.st <- x
        }
        if(j==2){
          x.at <- x
        }
      }
      ## calculate RLE, reference call (1-reference, 0-nonreference), log2(PM/maxMM)
      PM.indx <- sapply(dimnames(x)[[1]],function(a,b){
        which(b==a)
      },b=dimnames(x)[[2]])

      logPM <- refc <- dm <- am <- matrix(0,dim(x)[1],length(celnames))

      for(h in 1:dim(x)[3]){
        for(k in 1:dim(x)[1]){
          logPM[k,h] <- log2(x[k,PM.indx[k],h])
          if(mode.baselevel=="average"){
            if(which(x[k,,h]==max(x[k,,h]))[1]==PM.indx[k]){
              refc[k,h] <- 1
            }else{
              refc[k,h] <- 0
            }
          }
          dm[k,h] <- logPM[k,h]-max(log2(x[k,-PM.indx[k],h]))
          am[k,h] <- 1/2*(logPM[k,h]+max(log2(x[k,-PM.indx[k],h])))
        }
      }
      rle <- logPM-apply(logPM,1,median)

      all.rle[[unitindx[j]]] <- rle
      all.am[[unitindx[j]]] <- am
      all.logPM[[unitindx[j]]] <- logPM
      if(mode.baselevel=="average"){
        all.refc[[unitindx[j]]] <- refc
      }
      all.dm[[unitindx[j]]] <- dm
    }
    if(mode.baselevel=="refcallrate"){
      ref <- dimnames(x.st)[[1]]
  ##calcalate calls using a simple base caller as below for each sample
      refc <- NULL
      for(k in 1:dim(x.st)[3]){
        call <- call.ref.base(st=x.st[,,k],at=x.at[,,k])
        a <- mapply(identical,ref,call)
        refc <- cbind(refc,ifelse(a==TRUE,1,0))
    }
     all.refc[[unitindx[1]]] <- all.refc[[unitindx[2]]] <- refc
    }
  }
  names(all.logPM) <- names(all.refc) <- names(all.dm) <- names(all.am) <-
    names(all.rle) <- unitnames
  output <- list(all.logPM=all.logPM,all.refc=all.refc,all.dm=all.dm,
                 all.am=all.am,all.rle=all.rle)
  return(output)
}





