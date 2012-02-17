\name{reorganize.datafiles}
\alias{reorganize.datafiles}
\title{Function for reorganizing normalized data into array formats}
\description{
A function to reorganize the delta, sigma and base pair information into array formats.
}
\usage{
reorganize.datafiles(delta,sigma,pairs,n.bases,n.samples,sample.names,exon.names,ref.info)
}
\arguments{
  \item{delta}{A list of arrays with log quantile normalized intensities RM/AM ratios (delta) for exons.  Each array represents data for an exon with 3 dimention: Units, 3 alternative matches, and samples. }
  \item{sigma}{A list of arrays with log quantile normalized average intensities of RM and AM(sigma) for exons.  Each array represents data for an exon with 3 dimention: Units, 3 alternative matches, and samples.}
  \item{pairs}{a list of matrices with base pair information for exons. Each matrix represents base pair data for an exon with Units in rows and 3 reference match/alternative matche pairs in columns. }
  \item{n.bases}{Total number of bases.}
  \item{n.samples}{Total number of samples}
  \item{sample.names}{Names of the samples}
  \item{exon.names}{Names of exons}
  \item{ref.info}{A data frame with information of reference base (RM), data should have columns 'Fragment','FragPos','Ref'.}
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
#load object from function get.delta.sigma 
#data(all.qn)
#data(ref.info)
#unitnames <- names(all.qn$all.delta.qn) 
#exonnames <- unique(substr(unitnames,1,nchar(unitnames)-5)) 
#n.bases <- 0
#for(i in 1:length(exonnames)){
#  unitindx <- which(regexpr(exonnames[i],unitnames)>0)[1]
#  n.bases <- n.bases+nrow(all.qn$all.delta.qn[[unitindx]])
#}
#reorg.dat<- reorganize.datafiles(delta=all.qn$all.delta.qn,
#                                         sigma=all.qn$all.sigma.qn,
#                                         pairs=all.qn$all.pairs.qn,
#                                         n.bases=n.bases,
#                                         n.samples=length(all.qn$celnames),
#                                         sample.names=all.qn$celnames,
#                                         exon.names=exonnames,
#                                         ref.info=ref.info)

data(reorg.dat)
str(reorg.dat)
}



