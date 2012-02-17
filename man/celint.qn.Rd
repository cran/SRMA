\name{celint.qn}
\docType{data}
\alias{celint.qn}
\title{quantile normalized intensities for each base}
\description{
A list of matrices providing quantile normalized intensities for each
base. It's for 5 hapmap samples: NA11840, NA12156, NA12878, NA18507 and NA19625.}
\usage{celint.qn}
\format{
A list. Each element of the list contains a matrix of quantile
normalized intensities for an exon. Each row of the matrices corresponds to a probe feature and each colume corresponds to a sample.
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
data(celint.qn)
is.list(celint.qn)
length(names(celint.qn))
str(celint.qn[[1]])
}


