######
# get.delta.sigma: get quantile normalized delta and sigma for each base
# Wenyi Wang
#####

get.delta.sigma <- function(chipType,dataDir,chipOrder,celint,exonList=NA){

  ##load R packages if not yet loaded
  require(R.utils)
  require(R.filesets)
  require(aroma.core)
  require(aroma.affymetrix)

  ##read in annotation file
  cdf <- AffymetrixCdfFile$byChipType(chipType,tags="byExon");
  if(chipType=="MitoP-2r520651"){
    units <- 10:nbrOfUnits(cdf);
    cellQuartets <- getCellQuartets(cdf,units=units);
  }else{
    cellQuartets <- getCellQuartets(cdf);
  }
  unitnames <- names(cellQuartets)

  exons <- unique(c("AFFX-TagIQ-EX",unique(substr(unitnames[regexpr("ENSE",unitnames)>0],1,15))))

  if(sum(is.na(exonList))==0){
    subexons <- exons[exonList]
  }else{
    subexons <- exons
  }


  all.logPM.qn <- all.delta.qn <- all.sigma.qn <- all.pairs.qn <- vector("list",length(unitnames))
  for(i in 1:length(subexons)){
    ##print(i)
    unitindx <- which(regexpr(subexons[i],unitnames)>0)
    for(j in 1:2){
      cells <- t(cellQuartets[[unitindx[j]]][[1]])
      x <- celint[[unitindx[j]]]
      celnames <- colnames(x)
      dim(x) <- c(nrow(x)/4,4,ncol(x))
      dimnames(x) <- list(rownames(cells),c("T","G","C","A"),celnames)

      ## calculate RLE, reference call (1-reference, 0-nonreference),log2(PM/maxMM)
      PM.indx <- sapply(dimnames(x)[[1]],function(a,b){
        which(b==a)
      },b=dimnames(x)[[2]])

      logPM <- matrix(NA,dim(x)[1],dim(x)[3])
      delta <- sigma <- array(NA,c(dim(x)[1],3,dim(x)[3]))
      pairs <- matrix(NA,dim(x)[1],3)

      refbase <- rownames(x)
      for(k in 1:dim(x)[1]){
        for(h in 1:dim(x)[3]){
          logPM[k,h] <- log2(x[k,PM.indx[k],h])
          delta[k,,h] <- logPM[k,h]-log2(x[k,-PM.indx[k],h])
          sigma[k,,h] <- 1/2*(logPM[k,h]+log2(x[k,-PM.indx[k],h]))
        }
        
       ## if(regexpr(pattern="st",text=unitnames[unitindx[j]])>0){                                                                                                                                                                                                                       
          pairs[k,] <- unlist(lapply(paste(refbase[k],dimnames(x)[[2]][-PM.indx[k]],sep=""),complement))
       ##}else{                                                                                                                                                                                                                                                                         
       ##pairs[k,] <- paste(refbase[k],dimnames(x)[[2]][-PM.indx[k]],sep="")                                                                                                                                                                                                            
       ## }    
      }

      all.logPM.qn[[unitindx[j]]] <- logPM
      all.delta.qn[[unitindx[j]]] <- delta
      all.sigma.qn[[unitindx[j]]] <- sigma
      all.pairs.qn[[unitindx[j]]] <- pairs
    }
  }
  names(all.logPM.qn) <- names(all.delta.qn) <- names(all.sigma.qn) <- names(all.pairs.qn) <- unitnames
  output <- list(celnames,all.logPM.qn,all.delta.qn,all.sigma.qn,all.pairs.qn)
  names(output) <- c("celnames","all.logPM.qn","all.delta.qn","all.sigma.qn","all.pairs.qn")
  return(output)
}





