\name{position.post}
\alias{position.post}
\title{Function for calculating posterior probablilities and assign genotypes}
\description{
A function is used to calculate posterior probablility and assign reference homozygous, heterozygous and alternative homozygous to each site.  
This function is an internal function and is called by calc.pp.bp in this package.
}
\usage{
      position.post(indx1,ref,het,hom,bp.indx,alldelta,f1,f3,errors1,errors3,round)
}
\arguments{
  \item{indx1}{position or index number of site.}
  \item{ref}{Homozygous reference genotype.}
    \item{het}{Heterozygous genotype.}
    \item{hom}{Alternative homozygous genotype.}
    \item{bp.index}{The index of one of the 3 alternative bases other than the reference.}
    \item{alldelta}{All logratio data array with 4 dimentions: sites, alternative bases, samples, and directions.}
    \item{f1}{Fitted values for reference homozygous when chip effect was adjusted using sense strand data. }
    \item{f3}{Fitted values for alternative homozygous when chip effect was adjusted using antisense strand data. }
    \item{error1}{Standard deviation of residues for reference homozygous when chip effects were adjusted. }
    \item{error3}{Standard deviation of residues for alternative homozygous when chip effects were adjusted. }
    \item{round}{1 or 2, indicating during which iteration the posteriors and genotypes are calculated.}
  }
\value{
  A list of 2
  \item{post}{Posterior probabilites of reference homozygous, heterozygous and alternative homozygous}
  \item{genotype}{Genotype calls for the site and base, either 'Ref', 'Het', or 'Hom'}
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
\keyword{internal}
\examples{
#load('reorg.dat.rda')
#load("chipadj.r1.qn.rda")
#load('ref.info.rda')
#alldelta.qn<-reorg.dat$alldelta.qn
#allpairs.qn<-reorg.dat$allpairs.qn
#f1 <- chipadj.r1.qn$f1
#f3 <- chipadj.r1.qn$f3
#errors <- chipadj.r1.qn$errors
#snps <- ref.info
#names(snps)[1:2] <- c("ExonID","BasePos")
#indx1<-sample(1, 1:nrow(alldelta.qn))
#bp.ind<-sample(1,1:3)
#bad <- is.na(alldelta.qn[indx1,1,,1])
#basepair <- allpairs.qn[indx1,4-bp.ind,2]
#        ref.a <- substr(basepair,1,1)
#        alt.a <- substr(basepair,2,2)
#        ref <- paste(ref.a,ref.a,sep="")
#        het <- basepair
#        hom <- paste(alt.a,alt.a,sep="")
#output <- position.post(indx1=indx1,ref=ref,het=het,hom=hom,bp.indx=bp.indx,
#                                alldelta=alldelta.qn,f1=f1,f3=f3,
#                                errors1=errors,errors3=errors,round=1)
#str(output)                                
}