\name{get.qc.measures.amplevel}
\alias{get.qc.measures.amplevel}
\title{Function for calculatig amplicon specific QC measures}
\description{
A function is used to calculate seven metrics that can be used for
quality assessment for all amplicons and both strands: the median and
IQR of relative log expression, the reference call rate, the median and
MAD of log intensity ratios and the median and MAD of average intensities.
}
\usage{
get.qc.measures.amplevel(ampmap,celnames,all.rle,all.refc,
                         all.dm,all.am,ampList)}
\arguments{
  \item{ampmap}{a data frame, mapping between amplicons and
    exons. Structure of ampmap is: row names are 1 to the total number of amplicons, column names are the same as the data file {ampmap.rda}.}
  \item{celnames}{sample names, a list of character strings}
  \item{all.rle}{the difference between the RM log intensities and the
    median of the RM log intensities across all samples for each probe
    feature. The structure is a list of matrices. Each matrix
    corresponds to an exon. Each row of the matrices corresponds to a base position and each column corresponds to a sample.}
  \item{all.refc}{reference call indicators of 1 or 0. If the RM probe
    intensities is highest, the reference call number is 1. Otherwise
    0. The structure is the same as {all.rle}.}
  \item{all.dm}{the difference between the RM log intensities and the
    maximum of the AM log intensities. The structure is the same as {all.rle}.}
  \item{all.am}{average of the RM intensities and the maximum of the AM
    intensites. The structure is the same as {all.rle}.}
  \item{ampList}{a vector of scalers which is the list of the amplicons you
    want to extract from all amplicons. Default is NA, wich provides values
    for the full list of the amplicons.}
}
\value{
An array with dimensions (1+total number of amplicons)x7x2x(number of
samples).

  Dimension 1:

  The first one corresponds control fragment (AFFX-TagIQ-EX). The others
  correspond to the amplicons whose names have less than 17 letters.

  Dimension 2:

  1) median of rle

  2) IQR of rle

  3) reference call rate

  4) median of dm

  5) MAD of dm

  6) median of am

  7) MAD of am

  Dimension 3:

  1) sense strand

  2) antisense strand
}
\author{Wenyi Wang \email{wwang7@mdanderson.org}
        Yan Xu  
}
\references{
Wang W, Shen P, Thyagarajan S, Lin S, Palm C, Horvath R, Klopstock T,
Cutler D, Pique L, Schrijver I, Davis RW, Mindrinos M, Speed TP, Scharfe
C. Identification of Rare DNA Variants in Mitochondrial Disorders with
Improved Array-based Sequencing. Nucleic Acids Res. 2011 Jan;39(1):44-58.

Shen P, Wang W, Krishnakumar S, Chi AK, Palm C, Enns GM, Davis RW, Speed TP, 
Mindrinos M, Scharfe C. High quality DNA sequence capture of 524 disease candidate genes. 
Proc Natl Acad Sci U S A. 2011 Apr 19;108(16):6549-54.
}
\keyword{amplicon}
\keyword{exon}
\examples{
##load the ampmap,intensity and qc.measures.baselevel data
##  data(ampmap)
##  data(celint)
##  data(qc.measures.baselevel)

##  celnames <- colnames(celint[[1]])
##  y <- qc.measures.baselevel
##  all.rle <- y$all.rle
##  all.refc <- y$all.refc
##  all.dm <- y$all.dm
##  all.am <- y$all.am
##  ampList <- c(0,2,5)
##  allqm <- get.qc.measures.amplevel(ampmap=ampmap,celnames=celnames,
##                                     all.rle=all.rle,
##                                     all.refc=all.refc,
##                                     all.dm=all.dm,
##                                     all.am=all.am,ampList=ampList)
##  str(allqm)
}


