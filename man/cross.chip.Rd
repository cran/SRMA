\name{cross.chip}
\alias{cross.chip}
\title{Function for cross-chip genotype calling.}
\description{
This is the core function of SRMA algorithms. First kmeans algorithms are used to determine potential category of a position. There are three categories: common SNPs, rare SNPs, and ref only. An EM algorithm is used to generate genotype calls for common SNPs. A Bayesian approach is used to call rare SNPs.
}
\usage{
    cross.chip(indx1,post,ref,het,hom,ref.a,alt.a,frag,fragpos,
                       alldelta,allsigma,f1,f3,mu1,mu3,dbsnp.info,bad,
                       final.bp,plot=FALSE)
}
\arguments{
	\item{indx1}{index of the position}
  \item{post}{A data frame with posterior probability of ref allele homozygous, heterozygous and alternative allele homozygous for chosen base.}
  \item{ref}{Reference allele homozygous genotype.}
    \item{het}{heterozygous genotype.}
  \item{hom}{Alternative allele homozygous genotype.}
  \item{ref.a}{Reference allele genotype.}
  \item{alt.a}{Alternative allele genotype.}
  \item{frag}{Fragment ID that the chosen base is in}
    \item{fragpos}{Position of the chosen base in the fragment }
    \item{alldelta}{All logratio data array with 4 dimensions: sites, alternative bases, samples, and directions.}
  \item{allsigma}{All log average data array with 4 dimensions: sites, alternative bases, samples, and directions.}
    \item{f1}{Fitted values for homozygous reference when chip effects were adjusted.}
    \item{f3}{Fitted values for homozygous alternative when chip effects were adjusted.}
    \item{mu1}{Median of f1 across samples. }
    \item{mu3}{Median of f3 across samples. }
    \item{dbsnp.info}{A data frame with information of dbSNPs}
    \item{bad}{A vector showing the bad samples for the chosen base. }
    \item{final.bp}{The index of the principal alternative allele for chosen base.}
    \item{plot}{Whether plots should be generated}
    }
\value{
  a matrix with samples in row and 2 columns representing sense and antisense strand.
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
\keyword{cross.chip}
\examples{
#load('reorg.dat.rda')
#load("chipadj.r2.qn.rda")
#load('ref.info.rda')
#load('allpost.r2.rda')
#load("dbsnp.info.rda")
#load("bpfrompp1.qn.rda")
#allbp.indx <- bpfrompp1.qn$all.bp.indx
#alldelta.qn<-reorg.dat$alldelta.qn
#allpairs.qn<-reorg.dat$allpairs.qn
#f1 <- chipadj.r2.qn$f1
#f3 <- chipadj.r2.qn$f3
#mu1 <- chipadj.r2.qn$mu1
#mu3 <- chipadj.r2.qn$mu3
#snps <- ref.info
#indx1<-1
#frag <- substr(as.character(unlist(snps$ExonID[i])),1,15)
#fragpos <- snps$BasePos[i]
#indx1 <- which(ref.info$Fragment==frag&ref.info$FragPos==fragpos)
#bad <- is.na(alldelta[indx1,1,,1])
#basepair <- allpairs[indx1,4-allbp.indx[indx1],2]
#ref.a <- substr(basepair,1,1)
#alt.a <- substr(basepair,2,2)
#ref <- paste(ref.a,ref.a,sep="")
#het <- basepair
#hom <- paste(alt.a,alt.a,sep="")
#post <- allpost[i,!bad,]
#names(snps)[1:2] <- c("ExonID","BasePos")
#output<-cross.chip(indx1,post,ref,het,hom,ref.a,alt.a,frag,fragpos,
#                       alldelta,allsigma,f1,f3,mu1,mu3,dbsnp.info,bad,
#                       final.bp,plot=FALSE)
#str(output)                                
}