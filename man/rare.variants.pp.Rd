\name{rare.variants.pp}
\alias{rare.variants.pp}
\title{Function for calculating posterior probablity for SNPs with rare variants.}
\description{
This function calculates the posterior probability for the heterozygous with rare variant/alternative alleles. When the two strands give inconsistent ordering of posterior probabilities, the strand with smaller variance will be used to calculate the loglikelihood and posteriors. If one strand has log intensities lower than 9, then the other strand will be used in the calculation.

This function is an internal function called by cross.chip function.
}
\usage{
    rare.variants.pp(d,ai,bad,center,cov,prior, dbsnpID)
}
\arguments{
  \item{d}{A matrix with observed log ratios.}
  \item{ai}{A matrix with intensities (sigma) of the bases in an exon}
  \item{bad}{A vector indicates samples that do not have good intensity data at this base position}
  \item{center}{vector of 2 numbers indicating the expected log ratios for sense and antisense}
  \item{cov}{2 by 2 variance covariance matrix of log ratios}
  \item{prior}{prior probability of reference homozygous, heterozygous, and alternative homozygous} 
  \item{dbsnpID}{Containing an index numer if the base position is dbSNP, otherwise has a length of 0}
  } 
\value{
  A matrix of 3 columns with each column indicating posterior probability of reference homozygous, heterozygous, and alternative homozygous
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
\keyword{delta}
