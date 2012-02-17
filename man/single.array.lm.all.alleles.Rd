\name{single.array.lm.all.alleles}
\alias{single.array.lm.all.alleles}
\title{A function for adjusting for logratios for all bases on a single array}
\description{
This function is for single-chip adjustement in SRMA algorithms. It processes the data in a sample by sample basis in order to adjust for effects of amplicon length, GC content and average log intensities using piece-wise splines regression. This function takes log ratios for all three alternative alleles.
}
\usage{
    single.array.lm.all.alleles(reorg.Data, verbose=FALSE)
}
\arguments{
	  \item{reorg.Data}{A list with 6 elements:
	  alldelta.qn:array of normalized logratios for samples with 4
	  dimensions: positions,bases(3), samples, strands(2).
	  allsigma.qn: array of normalized log averages for samples with
	  4 dimensions: positions,bases(3), samples, strands(2).
	  allpairs.qn: array of genotypes with 3 dimensions:
	  positions,bases(3), strands(2).
	  ref.info:Information of reference call for all the
	  positions. A dataframe with 3 elements: Fragment(amplicon
	  name), FragPos(Position on the exon), Ref(Reference genotye of
	  the position).
	  gc: The GC content of the amplicons.
	  exonlength: The length of the amplicons.}
	  \item{verbose}{Whether to show verbose message.}
	  }
\value{A list with 5 elements:
  \item{f1}{predicted values for reference homozygous in single chip adjustment}
  \item{f3}{predicted values for alternative homozygous in single chip adjustment}
  \item{mu1}{average predicted values for reference homozygous in single chip adjustment}
  \item{mu3}{average predicted values for alternative homozygous in single chip adjustment}
  \item{errors}{Residuals}
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
\keyword{single.array.lm.all.alleles}
\examples{
data('reorg.dat')
str(reorg.dat)
sd.adj.dat<-single.array.lm.all.alleles(reorg.dat)
str(sd.adj.dat)                                
}