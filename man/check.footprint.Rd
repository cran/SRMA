\name{check.footprint}
\alias{check.footprint}
\title{Function for checking the footprint effect.}
\description{
This function is to check the footprint effects in the resequencing array. Following criteria are used to identify footprint effects:
1. Neighbor location should not be a common SNP
2. Whenever the center location is a homozygous variant call in a sample, the neighbor location is also a homozygous variant call.
3. Too many heterzyous calls in the neighbor of non-homozygous calls

This function is an internal function called by crosschip function.
}
\usage{
    check.footprint(snpsi,exoni,exonname,refCalls,dbsnp.info)
}
\arguments{
  \item{snpsi}{index of bases that are not ref}
  \item{exoni}{A list containing the output information of an exon/amplicon. An element of the output from funciton cross.chip.  }
  \item{exonname}{Name of the exon}
  \item{refCalls}{A dataframe with reference allele information of all bases. }
  \item{dbsnp.info}{A dataframe with information of dbSNPs covered by the array.}
    }
\value{
  Indices of snps without footprint effect.
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
\keyword{internal}
\keyword{footprint}
\examples{
#load('alloutput.rda')  
#load('refCalls.rda')
#load('dbsnp.info.rda')
#exoni <- alloutput[[1]]
#snpsi <- which(!exoni$refsites)
#exonname<-names(alloutput)[1]
#check.footprint(snpsi,exoni,exonname,refCalls,dbsnp.info)
}