####
## get.qc.measures.amplevel.R: calculate amplicon specific QC measures
## Wenyi Wang
####

get.qc.measures.amplevel <- function(ampmap=ampmap,celnames=celnames,
                                     all.rle=all.rle,
                                     all.refc=all.refc,
                                     all.dm=all.dm,
                                     all.am=all.am,
                                     ampList=NA){
  amp <- as.character(ampmap$amplicon)
  tmp <- nchar(amp)

  allqm <- array(0,c(length(amp)+1,7,2,length(celnames)))
  dimnames(allqm) <- list(c("AFFX-TagIQ-EX",amp),
                          c("med.rle","iqr.rle","ref.cr","med.dm","mad.dm",
                            "med.am","mad.am"),
                          c("sense","antisense"),
                          celnames)
  
  if(sum(is.na(ampList))>0){
    allqm[1,,,] <- get.qm(id=0,rle=all.rle,refc=all.refc,dm=all.dm,am=all.am,ampmap=ampmap)
    ampOrder <- 1:length(amp)
  }else{
    if(sum(ampList==0)>0){
      allqm[1,,,] <- get.qm(id=0,rle=all.rle,refc=all.refc,dm=all.dm,am=all.am,ampmap=ampmap)
    }
    ampOrder <- ampList[ampList!=0]
  }
  for(i in ampOrder){
    ##print(i)
    if(tmp[i]>17){
      print("Warning: this amplicon name has more letters than accepted.")
      next
    }else{
      allqm[1+i,,,] <- get.qm(id=i,
                              rle=all.rle,
                              refc=all.refc,
                              dm=all.dm,
                              am=all.am,
                              ampmap=ampmap,mode.amplevel="average")
    }
  }
  return(allqm)
}

