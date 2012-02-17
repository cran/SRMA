\name{reorg.dat}
\docType{data}
\alias{reorg.dat}
\title{Normalized dataset for downstream analysis}
\description{
An example of normalized array-based sequencing data. The data is the output of normalization module. 
}
\usage{
data(reorg.dat)
}
\value{
  A list with 6 components:
  \item{alldelta.qn}{An array with 4 dimentions containing reorganized log normalized intensity ratios of RM/AM.  The first dimension is  the unit. The second dimension is 3 RM/AM basepairs. The third dimension is sample. The forth dimension is sense or antisense hybridization.}
  \item{allsigma.qn}{An array with 4 dimentions containing reorganized log normalized intensity averages of RM/AM with the same data structure as \code{alldelta.qn}.}
  \item{allpairs.qn}{An array with 3 dimentions containing reorganized base pair information. Its structure is similar to \code{alldelta.qn}, but without the dimension representing samples. Element of the array corresponds to pairs of RM/AM allele, where RM stands for reference match and AM stands for alternative match.}
  \item{ref.info}{A data frame with 3 columns. Fragment is the exon name; FragPos is position of the base; and Ref is reference match base.} 
  \item{gc}{A vector with GC contents of each exon.}
  \item{exonlength}{A vector with lengths of each exon}
}

\author{Wenyi Wang \email{wwang7@mdanderson.org}
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
\keyword{manip}
\examples{
#reorg.dat<- get.delta.sigma.qn(qc.output, ref.info=ref.info)
data(reorg.dat)
str(reorg.dat)
}



