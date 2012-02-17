\name{single.chip.posterior.prob}
\alias{single.chip.posterior.prob}
\title{Wrapper function for adjustment of logratios for all bases using linear model}
\description{
This function is the wrapper function of single chip adjustement in SRMA algorithms. It processes the data in an sample by sample bases to adjust for effect of amplicon length, GC content and average log intensities using piece-wise splines regression. It first calls \code{single.array.lm.all.alleles} to calculate posterior probabilities for all three alternative alleles, then selects one alternative allele and re-calculate posterior probabilities by fitting a piece-wise splines regression model on log ratios of only the selected allele pairs. 
}
\usage{
    single.chip.posterior.prob(reorg.Data, verbose=FALSE)
}
\arguments{
	 \item{reorg.Data}{A list with 6 elements:
	  alldelta.qn: array of normalized logratios for samples with 4
	  dimensions: positions,bases(3), samples, strands(2). 
	  allsigma.qn: array of normalized log averages for samples with
	  4 dimensions: positions,bases(3), samples, strands(2).
	  allpairs.qn: array of genotypes with 3 dimensions:
	  positions,bases(3), strands(2).
	  ref.info: Information of reference allele for all the
	  positions. A dataframe with 3 elements: Fragment(amplicon
	  name), FragPos(Position on the exon), Ref(Reference allele of
	  the position).
	  gc: The GC content of the amplicons.
	  exonlength: The length of the amplicons.}
	  \item{verbose}{Whether to show verbose message.}
	  }
\value{A list with 3 elements:
  \item{post.prob}{A list of posterior probabilities for amplicons. Each component of the list corresponding to an amplicon is an array with 3 dimentions: the first dimentions is bases in the amplicon, the second dimention is samples, and the third dimention is genotypes of homozygous reference, heterozygous, and homozygous alternative allele.}
  \item{allbp.indx}{The index of principal alternative allele chosen based on the posterior probabilities.}
  \item{chipadj.qn}{The model parameters from single base adjustment of logratios. A list of 6 elements:   
  (1) Model fitted values for reference homozygous log ratios.
  (2) Model fitted values for alternative homozygous log ratios (extrapolated).
  (3) Median values of (1)
  (4) Median values of (2)
  (5) Residuals for reference homozygous
  (6) Residuals for alternative homozygous
  }
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
\keyword{single.chip.posterior.prob}
\examples{
data('reorg.dat')
str(reorg.dat)
\dontrun{sd.adj.dat<-single.chip.posterior.prob(reorg.dat)}
\dontshow{data('sd.adj.dat')}
str(sd.adj.dat)                                
}