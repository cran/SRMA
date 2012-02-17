get.qm <- function(id,rle,refc,dm,am,ampmap,mode.amplevel="average"){
  unitnames <- names(rle)
  qm <- array(NA,c(7,2,ncol(rle[[1]])))
  if(id==0){
    for(j in 1:2){
      qm[1,j,] <- apply(rle[[j]],2,median)
      qm[2,j,] <- apply(rle[[j]],2,IQR)
      qm[4,j,] <- apply(dm[[j]],2,median)
      qm[5,j,] <- apply(dm[[j]],2,mad)
      qm[6,j,] <- apply(am[[j]],2,median)
      qm[7,j,] <- apply(am[[j]],2,mad)
    }
    if(mode.amplevel=="average"){
      for(j in 1:2){
        qm[3,j,] <- apply(refc[[j]],2,sum)/nrow(refc[[j]])
      }
    }
    if(mode.amplevel=="maximum"){
      tmp <- ifelse(refc[[1]]+refc[[2]]>0,1,0)
      qm[3,1,] <- qm[3,2,] <- apply(tmp,2,sum)/nrow(tmp)
    }
  }else{
    frag <- as.character(unlist(ampmap$amplicon[id]))

    if(nchar(frag)==15){
      indx <- which(regexpr(pattern=ampmap$exon[id],unitnames)>0)
      startpos <- ampmap$exon_start[id]
      endpos <- ampmap$exon_end[id]
    }

    if(nchar(frag)==17){ ##only take the positions that are non-overlapping between amps
      exon <- as.character(ampmap$exon[id])
      ampid <- which(regexpr(pattern=exon,ampmap$amplicon)>0)
      o.ampid <- ampid[order(ampmap$exon_start[ampid])]
      a <- which(o.ampid==id)
      if(a==1){
        startpos <- ampmap$exon_start[o.ampid[a]]
        endpos <- ampmap$exon_start[o.ampid[a+1]]
      }else{
        if(a==length(o.ampid)){
          startpos <- ampmap$exon_end[o.ampid[a-1]]
          endpos <- ampmap$exon_end[o.ampid[a]]
        }else{
          startpos <- ampmap$exon_end[o.ampid[a-1]]
          endpos <- ampmap$exon_start[o.ampid[a+1]]
        }
      }
      indx <-  which(regexpr(pattern=exon,unitnames)>0)
    }

    for(j in 1:2){
      qm[1,j,] <- apply(rle[[indx[j]]][startpos:endpos,],2,median)
      qm[2,j,] <- apply(rle[[indx[j]]][startpos:endpos,],2,IQR)
      qm[4,j,] <- apply(dm[[indx[j]]][startpos:endpos,],2,median)
      qm[5,j,] <- apply(dm[[indx[j]]][startpos:endpos,],2,mad)
      qm[6,j,] <- apply(am[[indx[j]]][startpos:endpos,],2,median)
      qm[7,j,] <- apply(am[[indx[j]]][startpos:endpos,],2,mad)
    }
    if(mode.amplevel=="average"){
      for(j in 1:2){
        qm[3,j,] <- apply(refc[[indx[j]]][startpos:endpos,],2,sum)/(endpos-startpos+1)
      }
    }
    if(mode.amplevel=="maximum"){
      tmp <- ifelse(refc[[indx[1]]][startpos:endpos,]+refc[[indx[2]]][startpos:endpos,]>0,1,0)
      qm[3,1,] <- qm[3,2,] <- apply(tmp,2,sum)/nrow(tmp)
    }
  }

  return(qm)
}
