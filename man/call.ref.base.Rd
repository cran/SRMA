\name{call.ref.base}
\alias{call.ref.base}
\title{Function for a naive base caller used to calculate reference call rate}
\description{
A function is used to call base A, C, G or T according to whose
summation of log probe intensities for both strands is the largest for each base position}
\usage{call.ref.base(st,at)}
\arguments{
  \item{st}{A matrix of raw intensities for sense strand. Each row of the
    matrix corresponds to a base position and four columns correspond to probes that have the center base as A, C, G and T.}
  \item{at}{A matrix of raw intensities for antisense strand. The
    structure is the same as the argument {st}.}
}
\value{
A vector of characters containing base calls (A, C, G or T) for each base position.
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
\keyword{base caller}
\examples{
##load the intensity data
  st<-matrix(rnorm(400), 100, 4, dimnames=list(1:100, c('A','C','G','T')))
  at<-matrix(rnorm(400), 100, 4, dimnames=list(1:100, c('A','C','G','T')))
  basecaller1 <- call.ref.base(st,at)
  str(basecaller1)
}


