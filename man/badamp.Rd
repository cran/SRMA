\name{badamp}
\docType{data}
\alias{badamp}
\title{Matrix for finding out the bad amplicons}
\description{
A matrix. Each column is a vector of 0 or 1 for all amplcons in each
sample. Element 1 refers to a bad amplicon, 0 refers to a good one.
}
\usage{badamp}
\format{
A matrix in which each row corresponds to an amplicon, and each column
corresponds to a sample.
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
\keyword{datasets}
\examples{
data(badamp)
badamp[1:5,]
}
