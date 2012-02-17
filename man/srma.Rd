\name{srma}
\alias{srma}
\title{Function for Sequencing Robust Multi-array Analysis.}
\description{
This function processes the array data in an exon-by-exon bases. Perform cross-chip genotype calling and organize the output of cross-chip calls.
}
\usage{
    srma(snps,allpost,alldelta,allsigma,allpairs,allbp.indx,dbsnp.info,
                  f1,f3,mu1,mu3,refCalls, plot=FALSE, ...)
}
\arguments{
	  \item{snps}{A data frame with column Exon ID and base position information}
  \item{allpost}{An array of posterior probabilities.}
    \item{alldelta}{All log ratio data with 4 dimensions: sites, alternative bases, samples, and directions.}
    \item{allsigma}{All log average data with 4 dimensions: sites, alternative bases, samples, and directions.}
    \item{allpairs}{All base pair data with 3 dimensions: sites, alternative bases, and directions.}
    \item{allbp.indx}{The index of the principal alternative allele.}
    \item{dbsnp.info}{A data frame with information of dbSNPs}
    \item{f1}{Fitted values for homozygous reference when chip effects were adjusted.}
    \item{f3}{Fitted values for homozygous alternative when chip effects were adjusted.}
    \item{mu1}{Median of f1 across samples. }
    \item{mu3}{Median of f3 across samples. }
    \item{plot}{Whether plots should be generated}
    \item{refCalls}{Dataframe with reference allele information}
    \item{...}{Any arguments pass to postscript function, such as file name, location etc.}
    }
\value{A list with 8 elements:
  \item{calls}{Genotype calls}
  \item{posteriors}{Posterior probabilities}
  \item{category}{Categories of genotypes, either common snp, rare snp or ref only}
  \item{refsites}{Number of site are reference hogomozygous}
  \item{slscore}{Silhouette scores}
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
\keyword{srma}
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
#allsigma.qn<-reorg.dat$allsigma.qn
#f1 <- chipadj.r2.qn$f1
#f3 <- chipadj.r2.qn$f3
#mu1 <- chipadj.r2.qn$mu1
#mu3 <- chipadj.r2.qn$mu3
#snps <- ref.info
#names(snps)[1:2] <- c("ExonID","BasePos")
#output<-srma(snps,allpost,alldelta,allsigma,allpairs,allbp.indx,dbsnp.info,
#                  f1,f3,mu1,mu3,plot=FALSE)
#str(output)                                
}