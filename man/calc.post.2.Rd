\name{calc.post.2}
\alias{calc.post.2}
\title{Function for calculating posterior probablilities of reference homozygous or heterozygous or alternative homozygous genotypes}
\description{
A function is used to calculate posterior probablility of homozygous or heterozygous with reference or alternative alleles.  
This function is an internal function and is called by position.post, and calc.pp.bp in this package.
This function is called in the second round posteriors are calculated.
}
\usage{
      calc.post.2(d,prior,f1i,f3i,error1,error3)
}
\arguments{
  \item{d}{a data matrix with 2 columns, each row respresents log ratios for a site, a base pair in sense and misense strand of a sample}

  \item{prior}{The prior probablilities for reference homo, hetero, and alternative homo.}
    \item{f1i}{fitted value for site i on sense and antisense strand for reference homozygous, when chip effects were adjusted.}
    \item{f3i}{fitted value for site i on sense and anitsense strand for alternative homozygous, when chip effects were adjusted.}
    \item{error1}{Standard deviation of residuals for reference homozygous, when chip effects were adjusted. }
    \item{error3}{Standard deviation of residuals for alternative homozygous, when chip effects were adjusted. }
  }
\value{
  A matrix with 3 columns, each row represents the posterior probabilities of reference homozygous or heterozygous or alternative homozygous for a sample}

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
prior <- c(0.999,0.0009,0.0001)
n<-100
d<-cbind(rnorm(n, 2, 1), rnorm(n, 3, 2))
f3i<-f1i<-cbind(rnorm(n, 2,0.1), rnorm(n, 3, 0.1))
error1<-error3<- matrix(rep(1:2, 2*n),n,2, byrow=TRUE)
post<-calc.post.2(d, prior, f1i, f3i, error1, error3)
head(post)
}