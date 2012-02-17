######
# get.intensity: get raw or quantile normalized intensities for each base
# Wenyi Wang
# For each base, each strand: M_1, M_2, M_3, S_1, S_2, S_3, pm, base pairs
#####

get.intensity <- function(chipType,dataDir,chipOrder,ampmap=NA,badamp=NA,
                          exonList=NA,chip.normalize=FALSE,amp.normalize=FALSE){
  ##object "ampmap" provides mapping between amplicons and exons
  ##object "badamp" provides vectors of 0 and 1 for all amplicons in each sample
  
  ##load R packages if not yet loaded
  ## require(R.utils)
  ##require(R.filesets)
  ##require(aroma.core)
  ##require(aroma.affymetrix)
  
  ##read in raw cel files
  cdf <- AffymetrixCdfFile$byChipType(chipType,tags="byExon");
  acs <- AromaCellSequenceFile$byChipType(chipType);
  csR <- AffymetrixCelSet$byName(dataDir, cdf=cdf);
  log <- verbose <- Arguments$getVerbose(-10, timestamp=TRUE)
  if(chip.normalize==TRUE){
    qn1 <- QuantileNormalization(csR,subsetToUpdate=NULL,typesToUpdate=NULL,
                                 targetDistribution=NULL,subsetToAvg=NULL,
                                 typesToAvg=NULL)
    csR <- process(qn1,verbose=log)
  }
  
  if(chipType=="MitoP-2r520651"){
    units <- 10:nbrOfUnits(cdf);
    cellQuartets <- getCellQuartets(cdf,units=units);
  }else{
    cellQuartets <- getCellQuartets(cdf);
  }
  unitnames <- names(cellQuartets)
  ##  chipindx <- match(chipOrder,getFullNames(csR))
  
  ##Perform quantile.normalization at amplicon level
  celint <- vector("list",length(unitnames))
  names(celint) <- unitnames
  
  exons <- unique(c("AFFX-TagIQ-EX",substr(unitnames[regexpr("ENSE",unitnames)>0],1,15)))
  if(sum(is.na(exonList))==0){
    subexons <- exons[exonList]
  }else{
    subexons <- exons
  }
  for(i in 1:length(subexons)){
    ##print(i)
    ## read in probe intensities
    unitindx <- which(regexpr(subexons[i],unitnames)>0)
    
    ## get both sense and antisense data
    cells.st <- t(cellQuartets[[unitindx[1]]][[1]])
    x.st <- extractMatrix(csR,cells=cells.st)
    ## the layout is 1:n for probe1,n+1:2n for probe2,2n+1:3n for probe3, 3n+1:4n for probe4
    cells.at <- t(cellQuartets[[unitindx[2]]][[1]])
    x.at <- extractMatrix(csR,cells=cells.at)
    
    ## identify bad chips, do not normalize using these data
    if(amp.normalize==TRUE){
      x <- rbind(x.st,x.at)[,chipOrder]
      new.x <- amp.spec.qn(subexons[i],x,badamp,ampmap)
      colnames(new.x) <- colnames(x.st)
      a <- nrow(x.st)
      celint[[unitindx[1]]] <- new.x[1:a,]
      celint[[unitindx[2]]] <- new.x[(a+1):(2*a),]
    }else{
      celint[[unitindx[1]]] <- x.st
      celint[[unitindx[2]]] <- x.at
    }
  }
  return(celint)
}



