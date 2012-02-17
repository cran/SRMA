\name{call.max}
\alias{call.max}
\title{Function for calling the base with maximum intensity}
\description{
A function is used to select the probe with maximum intensity and
call the base among T, G, C,and A}
\usage{call.max(y)}
\arguments{
  \item{y}{a vector of numbers which contains four intensities
    correspond to probes that differ at the center in T, G, C and A}
}
\value{base name corresponding to the probe with the largest intensity}
\author{Wenyi Wang \email{wwang7@mdanderson.org	}
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
\keyword{internal}
\examples{
 y <-c(3,6,4,8)
 call <- call.max(y) # "A"
}


