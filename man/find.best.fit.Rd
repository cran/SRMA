\name{find.best.fit}
\alias{find.best.fit}
\title{Function to select the best model fit}
\description{
   A function to select the best model fit. This function is an internal function called by cross.chip function.
}
\usage{
    find.best.fit(bic,allfit)
}
\arguments{
  \item{bic}{A vector of three BIC values, based on EM clustering result with cluster number k=1, 2, and 3.}
  \item{allfit}{A list of fit objects provided by EM clustering.}
}
\value{
   A scaler, denotes which model is selected as the best fit.
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
\keyword{BIC}
