\name{calc.pp.bp}
\alias{calc.pp.bp}
\title{Function for calculating posterior probablilities for a list of base positions}
\description{
This function performs the first round of base calls, assigns the genotypes, and picks one alternative allele out of three choices. Single chip posterior probabilities are calculated. Then based on the posteriors for all arrays, one alternative allele will be chosen from the three bases other than reference base. 
}
\usage{
      calc.pp.bp(snps,refCalls,alldelta.qn,allpairs.qn,f1,f3,errors,verbose=F, ...)
}
\arguments{
  \item{snps}{A data frame with column Exon ID and base position information}
  \item{refCalls}{Reference genotype.}
    \item{alldelta.qn}{All logratio data array with 4 dimensions: sites, alternative bases, samples, and directions.}
    \item{allpairs.qn}{All base pair data array with 3 dimensions: sites, alternative bases, and directions.}
    \item{f1}{Fitted values for reference homozygous, when chip effect was adjusted. }
    \item{f3}{Fitted values for alternative homozygous, when chip effect was adjusted. }
    \item{errors}{Standard deviation of residuals when chip effects were adjusted. }
    \item{verbose}{Whether to write verbose messages.}
    \item{...}{Any arguments pass to postscript function, such as file name, location etc.}
  }
\value{
  A list of 3
  \item{all.post}{Posterior probabilites of reference homozygous, heterozygous and alternative homozygousfor a list of base positions}
  \item{all.geno}{Genotype calls for the base positions and samples, either 'Ref', 'Het', or 'Hom'}
  \item{all.bp.indx}{The index of the alternative allele/base chosen, to be used in the second round.}
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
\keyword{natural cubic spline}
\keyword{sigma}
\keyword{delta}\examples{
##load('reorg.dats.rda')
##load("chipadj.r1.qn.rda")
##load('ref.info.rda')
##alldelta.qn<-reorg.dat$alldelta.qn
##allpairs.qn<-reorg.dat$allpairs.qn
##f1 <- chipadj.r1.qn$f1
##f3 <- chipadj.r1.qn$f3
##errors <- chipadj.r1.qn$errors
##snps <- ref.info
##names(snps)[1:2] <- c("ExonID","BasePos")
##output<-calc.pp.bp(snps,ref.info,alldelta.qn,allpairs.qn,f1,f3,errors)
##str(output)                                
}