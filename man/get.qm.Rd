\name{get.qm}
\alias{get.qm}
\title{Function for calculating quality measures for an amplicon}
\description{
A function is used to calculate, for each amplicon and each strand,
seven metrics that can be used for quality assessment: the median and
IQR of relative log expression, the reference call rate, the median and
MAD of log intensity ratios and the median and MAD of average intensities.
}
\usage{get.qm(id,rle,refc,dm,am,ampmap,mode.amplevel="average")}
\arguments{
  \item{id}{a scalar as index of an amplicon. Its value is between 0 and
    the total number of amplicons. The id for control fragment
    (AFFX-TagIQ-EX) is 0. The order of querried amplicon matches its
    location in the data file {ampmap.5um.rda}.}
  \item{rle}{the difference between the RM log intensities and the
    median of the RM log intensities across all samples for each probe
    feature. It's a list of matrices, and each matrix corresponds to an
    exon. Each row of the matrices corresponds to a base position and
    each column corresponds to a sample.}
  \item{refc}{reference call indicators of 1 or 0. If the RM probe
    intensities is highest, then the reference call indicator is
    1. Otherwise 0. The structure is the same as the argument {rle}.}
  \item{dm}{the difference between the RM log intensities and the
    maximum of the AM log intensities. The structure is the same as the
    argument {rle}.}
  \item{am}{average of the RM log intensities and the maximum of the
    AM log intensites. The structure is the same as the argument {rle}.}
  \item{ampmap}{Mapping between all amplicons and exons. Structure of ampmap
    is: row names are 1 to the total number of amplicons, column names
    are the same as the data file {ampmap.rda}.}
  \item{mode.amplevel}{choices of mode.amplevel are "average" and "maximum". Default is "average". Mode "average" means the reference call rate is calculated by averaging values for sense strand and antisense strand. Mode "maximum" means the reference call rate is calculated by taking the maximum of values for sense strand and antisense strand.}
}
\value{
  An array with dimensions 7x2x(number of samples).

  Dimension 1:

  1) median of rle

  2) IQR of rle

  3) reference call rate

  4) median of dm

  5) MAD of dm

  6) median of am

  7) MAD of am

  Dimension 2:

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
##  data(qc.measures.baselevel)
##  data(ampmap)
##  y <- qc.measures.baselevel
##  all.rle <- y$all.rle
##  all.refc <- y$all.refc
##  all.dm <- y$all.dm
##  all.am <- y$all.am
##  mode.amplevel <- "average"
##  qm <- get.qm(id=1,rle=all.rle,refc=all.refc,
##               dm=all.dm,am=all.am,ampmap=ampmap,mode.amplevel=mode.amplevel)
##  str(qm)
}


