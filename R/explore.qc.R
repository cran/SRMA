######
## explore.qc.R: plot quality control measures for a sample
## Wenyi Wang
######

explore.qc <- function(qc,sample.ID,fig.dir,fig.name){

  bad <- qc$badamp[,sample.ID]
  d <- (qc$D[,1,sample.ID]+qc$D[,2,sample.ID])/2
  t <- (qc$T[,1,sample.ID]+qc$T[,2,sample.ID])/2

  postscript(paste(fig.dir,fig.name,sep=""))
  par(mar=c(4,5,5,5),pty="s")
  plot(d~t,xlab="Average intensity",ylab="Log ratio",
       xlim=c(6,13),ylim=c(-1,3),cex=0.7)
  points(d[bad==1]~t[bad==1],pch=23,cex=0.7,col="grey")
  dev.off()
}
