\name{qn}
\alias{qn}
\title{Quantile normalization function}
\description{
A generic function which normalizes a matrix of probe level intensities using a normalization based upon quantile
}
\usage{qn(cels)}
\arguments{
\item{cels}{probe level intensity, a matrix of intensities where each
  row corresponds to a probe feature and each column corresponds to a
  unit, such as an amplicon, a chip or a sample.}
}
\value{
A matrix of quantile normalized intensities where each row corresponds
to a probe feature and each column corresponds to a unit, such as an amplicon, a chip or a sample.
}
\author{Wenyi Wang \email{wwang7@mdanderson.org}
        Yan Xu  
}
\references{
Bolstad, B.M.,Irizarry R.A., Astrand, M, and Speed, T.P.(2003)
A Comparison of Normalization Method for High Density
Oligonucleotide Array Data Based on Bias and Variance.
Bioinformatics 19(2), pp 185-193.
}
\keyword{quantile normalizaion}
\keyword{internal}
\examples{
a<-matrix(sample(0:100,12),3,4)
a
qn(a)
}


