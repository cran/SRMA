\name{correct.lr.singlebp}
\alias{correct.lr.singlebp}
\title{Function for correcting logratios for single chip with determined principal alternative alelles.}
\description{
This function performs correction of logratio data with determined principal alternative alleles. Using all chip information, a principal alternative allele has been chosen. This function will use logratios of determined alternative allele to reference allele, perform linear regression chip-by-chip, in order to adjust for GC content, average log intensity, fragment length, and base pair composition. Only the bases with reference homozygous posterior probability >0.999 are used for robustness.
}
\usage{
      correct.lr.singlebp(i,alldelta,allsigma,fragl=NULL,allpairs,allgc,allbp.indx,allpost, verbose=F, ...)
}
\arguments{
  \item{i}{A sample index indicating which sample will be processed.}

    \item{alldelta}{All logratio data array with 4 dimensions: sites, alternative bases, samples, and directions.}
    \item{allsigma}{All average log2 intensity data. An array with 4 dimensions: snp sites, alternative bases, samples, and directions.}
\item{allpairs}{base pairs for all bases. Its structure is a list
    of matrices. Each row of the matrices corresponds to a base position,
    and the three columns correspond to the three pairs of RM/AM allele,
    where RM stands for reference match and AM stands for alternative match.}
    \item{fragl}{Fragment length. A matrix with identical 2 columns for sense and antisense strands. If NULL, fragment length will not be adjusted.}
    \item{allgc}{GC content.A matrix with identical 2 columns for sense and antisense strands. }
    \item{allbp.indx}{A vector showing the index of the principal alternative allele out of the 3 choices. }
    \item{allpost}{Posterior probabilities of reference homozygous, heterozygous, and alternative homozygous genotypes. }
    \item{verbose}{Whether the intermediate results will be printed}
    \item{...}{Arguments passed to function png}
  }
\value{
  A list of 5 "f1","f3","mu1","mu3","error1","error3"
  \item{f1}{Fitted values for reference homozygous.}
  \item{f3}{Fitted values for alternative homozygous.}
  \item{mu1}{The median of fitted values for reference homozygous}
  \item{mu3}{The median of fitted values for alternative homozygous}
  \item{error1}{The standard deviation of residuals for reference homozygous}
  \item{error3}{The standard deviation of residuals for alternative homozygous}
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
\keyword{delta}
\examples{
#load('reorg.dat.rda')
#load("chipadj.r1.qn.rda")
#load('ref.info.rda')
#load('bpfrompp1.qn.rda')
#fragl<- rep(reorg.dat$exonlength, reorg.dat$exonlength)
#alldelta.qn<-reorg.dat$alldelta.qn
#allsigma.qn<-reorg.dat$allsigma.qn
#allpairs.qn<-reorg.dat$allpairs.qn
#allgc<-rep(reorg.dat$gc, reorg.dat$exonlength)
#allbp.indx <- bpfrompp1.qn$all.bp.indx
#allpost <- bpfrompp1.qn$all.post
##choose sample i
#i <- 1 
#good <- !is.na(alldelta.qn[,1,i,1])
#output <- correct.lr.singlebp(i=1,alldelta=alldelta.qn,allsigma=allsigma.qn,
#                              fragl=fragl,allpairs=allpairs.qn,
#                              allgc=allgc,allbp.indx=allbp.indx,allpost=allpost,
#                              verbose=F) 
#str(output)                                
}