\name{find.badamps}
\alias{find.badamps}
\title{Function for finding bad amplicons}
\description{
A function is used to find the bad amplicons with the criterion: the
reference call rate R<0.9}
\usage{find.badamps(allqm, mode.amplevel="average")}
\arguments{
  \item{allqm}{An array with dimension (number of amplicons)x7x2x(number
    of samples).

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
  \item{mode.amplevel}{choices of mode are "average" and
    "maximum". Default is "average". Mode "average" means the reference call rate is calculated by averaging values for sense strand and antisense strand, which is the percentage of positions where the intensities from the PM probes are the highest. Mode "maximum" means the reference call rate is calculated by taking the maximum of values for sense strand and antisense strand.}
}
\value{
A matrix. Each column is a vector of 0 or 1 for all amplcons in each
sample. Element 1 refers to a bad amplicon, 0 refers to a good one.
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

\keyword{bad amplicon}

\examples{
  ##load the allqm data
 ##  data(allqm)
 ##  mode.amplevel <- "average"
 ##  badamp <- find.badamps(allqm=allqm,mode.amplevel=mode.amplevel)
 ##  str(badamp)
}


