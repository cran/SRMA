\name{calc.slscore}
\alias{calc.slscore}
\title{Function for calculating Silhouette scores.}
\description{
This function calculates Silhouette width according to genotype information. The median of difference are used to calcuate dissimilarites of a sample between and within clusters.

This function is an internal function called by cross.chip function.
}
\usage{
    calc.slscore(d,genotypes)
}
\arguments{
  \item{d}{matrix of logratios for sense and antisense strands. Samples in rows and 2 columns}
  \item{genotypes}{genotype calls or clusters either 1, or 2, or 3. }
    }
\value{
  a matrix with samples in row and 2 columns representing sense and antisense strand.
      }
\author{Wenyi Wang \email{wwang7@mdanderson.org}
Nianxiang Zhang \email{nzhang@mdanderson.org}
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
\keyword{methods}
\keyword{internal}
\keyword{calc.slscore}
\examples{
d<-matrix(rnorm(100), 50, 2)
g<-sample(c(1:3), 50, replace=TRUE)
summary(calc.slscore(d, g))
}