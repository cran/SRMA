\name{ref.info}
\docType{data}
\alias{ref.info}
\title{Example dataset containing reference allele information that are used in SRMA pipeline of analyzing array-based sequencing array }
\description{
An example reference allele dataset used in SRMA pipeline of analyzing array-based sequencing array. All the columns in the data set are required for the SRMA analysis pipeline.
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
\keyword{ref.info}
\examples{
data('ref.info')
head(ref.info)
}