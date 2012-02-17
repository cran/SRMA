######
# quality.control: identify bad amplicons with the criterion: R<0.9
#####

quality.control <- function(chipType,dataDir,chipOrder,
                            output.measures=c("R","D","T"),
                            mode.baselevel="refcallrate",
                            mode.amplevel="average",
                            ampmap=NA,exonList=NA){

  celint <-get.intensity(chipType=chipType,dataDir=dataDir,chipOrder=chipOrder,
                         exonList=exonList)
  
  celnames <- colnames(celint[[1]])
  ##calculate processing time
  time <- round(12*length(celnames)/10)
  print(paste("Processing ",length(celnames),
              " samples... will take ~",time," hours"))
  
  baselevel <- get.qc.measures.baselevel(celint=celint,
                                         chipType=chipType,
                                         exonList=exonList,
                                         celnames=celnames,
                                         mode.baselevel=mode.baselevel)
  
  all.rle <- baselevel$all.rle
  all.refc <- baselevel$all.refc
  all.dm <-  baselevel$all.dm
  all.am <-  baselevel$all.am
  
  ##identify corresponding amplicons for the exonList
  if(sum(is.na(exonList))==0){
    unitnames <- names(celint)
    exons <- c("AFFX-TagIQ-EX",
               unique(substr(unitnames[regexpr("ENSE",unitnames)>0],1,15)))
    subexons <- exons[exonList]
    ampList <- c(0,which(ampmap$exon%in%subexons))
  }else{
    ampList <- NA
  }
  
  allqm <- get.qc.measures.amplevel(ampmap=ampmap,celnames=celnames,
                                    all.rle=all.rle,
                                    all.refc=all.refc,
                                    all.dm=all.dm,
                                    all.am=all.am,
                                    ampList=ampList)
  
  output <- list(NULL)
  badamp <- find.badamps(allqm,mode.amplevel)
  output <- c(output,list(badamp=badamp))
  if("R"%in%output.measures){
    R <- allqm[,3,,]
    output <- c(output,list(R=R))
  }
  if("D"%in%output.measures){
    D <- allqm[,4,,]
    output <- c(output,list(D=D))
  }
  if ("T"%in%output.measures){
    T <- allqm[,6,,]
    output <- c(output,list(T=T))
  }
  output$par<-list(chipType=chipType,dataDir=dataDir,chipOrder=chipOrder, ampmap=ampmap)
  return(output)
}



