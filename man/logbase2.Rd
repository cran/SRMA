\name{logbase2}
\alias{logbase2}
\title{Function for taking log-transforming on the base 2 scale}
\description{
A function is used to log-transform x on the base 2 scale.}
\usage{logbase2(x)}
\arguments{
  \item{x}{a vector of numbers}
}
\value{A vector of numbers}
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
\keyword{internal}
\examples{
 x <-c(2,-4,8)
 y <-logbase2(x)
 print(y)
}


