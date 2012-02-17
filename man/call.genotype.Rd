\name{call.genotype}
\alias{call.genotype}
\title{Wrapper function of performing multi-array analysis}
\description{
This function is the wrapper function for multi-array analysis using SRMA algorithm. In the analysis, the minor allele counts (MAC) are calculated and different approaches are used based on the value of MAC.
}
\usage{
    call.genotype(reorg.Dat, singlechip.pp, dbsnp.info=NULL)
}
\arguments{
	  \item{reorg.Dat}{Normalized data from normalization module.}
	  \item{singlechip.pp}{Posterior probabilities from single array analysis. A list of 3 elements: (1) posterior probabilites (2) indices of principal alternative alleles (3) model paramenters in adjustment of logratios.}
	  \item{dbsnp.info}{A dataframe with dbSNP information, containing column names of 'ExonID', 'BasePos','Ref', 'Alt.Allele', 'rs'.}
	  }
\value{A list with 3 elements:
  \item{output}{A list of results from SRMA for amplicons. Each component of the list corresponds to an amplicon and is list of 5 elements: (1)a matrix of genotype calls for each postion (row) and each sample (column) (2) posterior probabilities of the genotype calls (3) genotype category of the positions, rare variants, ref only, common SNPs w/EM. (4)whether the site is the reference site (5) Silhoutte scores array with 3 dimentions corresponding to position, sample, and strand}
  \item{sample.names}{The names of the samples.}
  \item{refCalls}{A data frame with information of reference base (RM)with columns 'Fragment','FragPos','Ref'.}
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
\keyword{call.genotype}
\examples{
data('reorg.dat')
str(reorg.dat)
\dontrun{sd.adj.dat<-single.chip.posterior.prob(reorg.dat)}
\dontshow{data('sd.adj.dat')}
srma.out<-call.genotype(reorg.dat,sd.adj.dat )
str(srma.out)                                
}