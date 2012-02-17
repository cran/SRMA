\name{complement}
\alias{complement}
\title{Base Pairing Function}
\description{
In the canonical Watson-Crick base pairing, A forms a base pair with T,
as and G with C in DNA. This function helps finding the matching base pair(s).
}
\usage{complement(x)}
\arguments{
  \item{x}{input base}
}
\value{
A character string which contains one letter if the input base is in
one-letter format or 2 letters if the input base is in two-letter format.
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
\keyword{base pair}
\keyword{internal}

\examples{
complement("A")
complement("AG")
}


