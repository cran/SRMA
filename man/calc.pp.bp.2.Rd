\name{calc.pp.bp.2}
\alias{calc.pp.bp.2}
\title{Function for second round calculation of posterior probablilities for a list of base positions.}
\description{
This function performs the second round of base calls using the alternative allele chosen in the first round. Single chip posterior probabilities are calculated, then using the posteriors for all arrays, genotype calls are determined.  
}
\usage{
    calc.pp.bp.2(snps,refCalls,alldelta,allpairs,f1,f3,errors1,errors3,allbp.indx,...)
}
\arguments{
  \item{snps}{A data frame with column Exon ID and base position information}
  \item{refCalls}{Reference genotype.}
  \item{alldelta}{All logratio data array with 4 dimensions: sites, alternative bases, samples, and directions.}
  \item{allpairs}{All base pair data array with 3 dimensions: sites, alternative bases, and directions.}
   \item{f1}{Fitted values for reference homozygous, when chip effect was adjusted. }
    \item{f3}{Fitted values for alternative homozygous, when chip effect was adjusted. }
    \item{errors1}{Standard deviation of residuals for reference homozygous,when chip effects were adjusted. }

    \item{errors3}{Standard deviation of residues for alternative homozygous, when chip effects were adjusted using antisense strand data. }
    \item{allbp.indx}{The index of the selected alternative allele.}
    \item{...}{Any arguments pass to postscript function, such as file name, location etc.}
  }
\value{
  A list of 2
  \item{all.post}{Posterior probabilites of reference homozygous, heterozygous and alternative homozygous for the base positions}
  \item{all.geno}{Genotype calls for the base positions and samples, either 'Ref', 'Het', or 'Hom'}
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
#load('reorg.dat.rda')
#load("chipadj.r2.qn.rda")
#load('ref.info.rda')
#alldelta.qn<-reorg.dat$alldelta.qn
#allpairs.qn<-reorg.dat$allpairs.qn
#f1 <- chipadj.r2.qn$f1
#f3 <- chipadj.r2.qn$f3
#errors1 <- chipadj.r2.qn$errors1
#errors3 <- chipadj.r2.qn$errors3
#snps <- ref.info
#names(snps)[1:2] <- c("ExonID","BasePos")
#output<-calc.pp.bp.2(snps,ref.info,alldelta.qn,allpairs.qn,f1,f3,errors1,errors3,allbp.indx)
#str(output)                                
}