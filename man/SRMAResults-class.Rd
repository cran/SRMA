\name{SRMAResults-class}
\Rdversion{1.1}
\docType{class}
\alias{SRMAResults-class}
\alias{initialize,SRMAResults-method}
\alias{show,SRMAResults-method}
\alias{summary,SRMAResults-method}

\title{Class "SRMAResults"}
\description{
SRMAResults class is the result of post.process module. It includes the final SNP calls and indel calls from SRMA algorithms.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SRMAResults", footprint, srmaMat, srmaIndelsnps, srmaSnps, sampleNames, slscoreCutoff)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{footprint}:}{Object of class \code{"list"} number of footprint and number of samples with footprint }
    \item{\code{srmaMat}:}{Object of class \code{"list"} matrices (position in row, sample in column) recording SRMA analysis results. (1) base quality (good:T, bad=F) (2) genotype calls (3) quality scores (4) Silhoutte scores. }
    \item{\code{srmaIndelsnps}:}{Object of class \code{"data.frame"} A dataframe with identified insertion/deletion calls. }
    \item{\code{srmaSnps}:}{Object of class \code{"data.frame"} A dataframe with identified SNV calls. }
    \item{\code{sampleNames}:}{Object of class \code{"vector"} A vector of sample names. }
    \item{\code{slscoreCutoff}:}{Object of class \code{"numeric"} A cutoff value for Silhoutte score to filter out unreliable calls. }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "SRMAResults")}: ... }
    \item{show}{\code{signature(object = "SRMAResults")}: ... }
    \item{summary}{\code{signature(object = "SRMAResults")}: ... }
	 }
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
\author{
  Wenyi Wang \email{wwang7@mdanderson.org}
  Nianxiang Zhang \email{nzhang@mdanderson.org}%%  ~~who you are~~
}
\examples{
showClass("SRMAResults")
data('reorg.dat')
\dontrun{sd.adj.dat<-single.chip.posterior.prob(reorg.dat)}
\dontshow{data('sd.adj.dat')}
srma.out<-call.genotype(reorg.dat,sd.adj.dat)
Results<-post.process(srma.out)
Results
summary(Results)
}
\keyword{classes}
