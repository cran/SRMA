\name{calc.ll}
\alias{calc.ll}
\title{Function for calculating log likelihood for univariatea normal distribution}
\description{
A function is used to calculate the log likelihood function of a normally distributed random variable with known mean and standard deviation. 
This function is an internal function and is called by calc.post, position.post, and calc.pp.bp in this package.
}
\usage{
      calc.ll(x,center,error)
}
\arguments{
  \item{x}{Random variables generated from a univariate normal distribution}
  \item{center}{Mean of the normal distribution}
  \item{error}{Standard deviation of the normal distribution.}
  }
\value{
  log likelihood
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
\keyword{internal}
\examples{
x<-rnorm(1)
calc.ll(x, 0, 1)
}