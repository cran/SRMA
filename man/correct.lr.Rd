\name{correct.lr}
\alias{correct.lr}
\title{Function for performing the adjustment of logratios using linear regression}
\description{
A function is used to calculate the model fitted values, predicted values, mean of the fitted and predicted values for the log ratios. The calculations are based on an additive model to delta using natural cubic splines on sigma (df=5), fragment length(df=3), gc content(df=3), and RM/AM basepair as explanatary variables. 
}
\usage{
      correct.lr(i,alldelta.qn,allsigma.qn,fragl=NULL,allpairs.qn,allgc, verbose=FALSE)
}
\arguments{
  \item{i}{sample name or sample index}
  
  Following arguments are components of object from \code{reorganize.datafiles}.
  
  \item{alldelta.qn}{log quantile normalized intensity ratios for each
    base. The structure is a list of arrays of three dimensions. The
    first dimension corresponds to the base positions. The second
    dimension corresponds to the three pairs of RM/AM allele, where RM
    stands for reference match and AM stands for alternative match. And the third dimension corresponds to samples.}
  \item{allsigma.qn}{log quantile normalized average intensities for each
    base. The structure is the same as the argument {alldelta.qn}.}
  \item{fragl}{fragment length, if NULL, fragment length will not be adjusted}
  \item{allpairs.qn}{base pairs for all bases. Its structure is a list
    of matrices. Each row of the matrices corresponds to a base position,
    and the three columns correspond to the three pairs of RM/AM allele,
    where RM stands for reference match and AM stands for alternative match.}
  \item{allgc}{exon level gc content}
  \item{verbose}{defalt value is FALSE.}
}
\value{
  A list of 5 containing:
  (1) Model fitted values for delta.
  (2) Predicted values for delta
  (3) Mean values of the model fitted values for delta
  (4) Mean values of the predicted values for delta
  (5) sigma
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
\keyword{natural cubic spline}
\keyword{sigma}
\keyword{delta}
\examples{
 ## data(reorg.dat)
 ## alldelta.qn<-reorg.dat$alldelta.qn
 ## allsigma.qn<-reorg.dat$allsigma.qn
 ## fragl<-reorg.dat$exonlength
 ## allpairs.qn<-reorg.dat$allpairs.qn
 ## allgc<-rep(reorg.dat$gc, reorg.dat$exonlength)
 ## i <- 5
 ##############
 ## output <-correct.lr(i,all.delta.qn,all.sigma.qn,fragl,allpairs.qn,allgc)
 ## str(output)
}