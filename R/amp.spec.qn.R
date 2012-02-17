amp.spec.qn <- function(exon.i,x,badamp,ampmap){
  ##1. when an exon corresponds to one amplicon (1 to 1, or 2 to 1)
  ##   apply quantile normalization on the exon, with bad exons removed
  ##2. when an exon correponds to multiple amplicons (1 to 2, 3, 4...)
  ##   translate the exon into the amplicon structure,
  ##   normalize on the amplicon, then put the data back to exon,
  ##   for the overlapping positions, take average when both are good
  ##3. when two exons corresponds to 3 amplicons (only one case)
  ##   the ENSE00001269987 show be normalized directly with
  ##   ENSE00001269987-ENSE00001270109_1 qc data
  ##   the ENSE00001270109 will be transformed into the amplicon structure,
  ##   normalized and transformed back

  amps <- as.character(ampmap$amplicon[regexpr(exon.i,ampmap$amplicon)>0])
  ##there are two types of amplicons, 1. full exon, 2, part of an exon
  category <- ifelse(nchar(amps[1])==15|exon.i=="AFFX-TagIQ-EX",1,2)

  if(category==1){
    if(exon.i=="AFFX-TagIQ-EX"){
      bad <- rep(0,ncol(x))
    }else{
      bad <- badamp[which(rownames(badamp)==amps[1]),]
    }
    x.qn <- matrix(NA,nrow(x),ncol(x))
    if(sum(bad==0)>0){
      x.qn[,bad==0] <- qn(x[,bad==0,drop=FALSE])
    }
  }

  if(category==2){
    x.qn <- matrix(NA,nrow(x),ncol(x))
    rows <- which(ampmap$exon==exon.i)
    startpos <- ampmap$exon_start[rows]
    endpos <- ampmap$exon_end[rows]
    for(j in length(amps)){
      bad <- badamp[which(rownames(badamp)==amps[j]),]
      l <- nrow(x)/8
      seg <- startpos[j]:endpos[j]
      st.pos <- c(startpos[j]:endpos[j],l+startpos[j]:endpos[j],l*2+startpos[j]:endpos[j],
                  l*3+startpos[j]:endpos[j])
      at.pos <- l*4+st.pos
      subx <- x[c(st.pos,at.pos),]
      stopifnot(nrow(subx)==(endpos[j]-startpos[j]+1)*8)

      subx.qn <- matrix(NA,nrow(subx),ncol(subx))
      if(sum(bad==0)>0){
        subx.qn[,bad==0] <- qn(subx[,bad==0,drop=FALSE])
      }
      x.qn[c(st.pos,at.pos),] <- ssum(x.qn[c(st.pos,at.pos),],subx.qn)
    }
  }
  return(x.qn)
}
