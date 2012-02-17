\name{ssum}
\alias{ssum}
\title{Average function for two matrices}
\description{
A simple function is used to get an average matrix of two matrices,
element by element. When one element is NA and the other is not, then
the function takes the value from the non-NA element. The two matrices
must at least have same number of elements.
}
\usage{
     ssum(x1,x2)
}
\arguments{
  \item{x1}{a  matrix}
  \item{x2}{a  matrix}
}
\value{
A matrix which elements are 1) the average of the two matrices, 2) the
value of the non-NA element, 3) NA if elemnts from both matrices are NA.
}
\author{Wenyi Wang \email{wwang7@mdanderson.org}
        Yan Xu   
}
\references{
Wang W, Shen P, Thyagarajan S, Lin S, Palm C, Horvath R, Klopstock T,
Cutler D, Pique L, Schrijver I, Davis RW, Mindrinos M, Speed TP, Scharfe
C. Identification of Rare DNA Variants in Mitochondrial Disorders with
Improved Array-based Sequencing. Nucleic Acids Res. 2011 Jan;39(1):44-58.

}
\keyword{matrix}
\keyword{internal}
\examples{
# x1 and x2 both have no missing values
x1<-matrix(c(1:9),3,3)
x2<-matrix(c(rep(2,9)),3,3)
x1
x2
ssum(x1,x2)

#x2 has a missing value
x1<-matrix(c(2,5,NA,1),2,2)
x2<-matrix(c(rep(2,4)),2,2)
x1
x2
ssum(x1,x2)
}



