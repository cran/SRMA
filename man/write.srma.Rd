\name{write.srma}
\alias{write.srma}
\title{Function to write final SNV calls and Indel calls to tables}
\description{
This function is to save the SRMA analysis results to tables. Two tables with SNP calls and insertion/deletion calls will be generated.
}
\usage{
    write.srma(x, file='SrmaResults',...)
}
\arguments{
	  \item{x}{An SRMAResults object from post.process function.}
	  \item{file}{Target files to save the results.}
	  \item{...}{Arguments passed to function write.table.}
	  }
\value{Two files tagged with 'SNP' and 'Indel' corresponding to SNP calls and insertion deletion calls.
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
\keyword{write.srma}
\examples{
data('reorg.dat')
\dontrun{sd.adj.dat<-single.chip.posterior.prob(reorg.dat)}
\dontshow{data('sd.adj.dat')}
srma.out<-call.genotype(reorg.dat,sd.adj.dat )
srma.results<-post.process(srma.out)
write.srma(srma.results)
}