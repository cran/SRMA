\name{writeVCF}
\alias{writeVCF}
\title{Output variant calls in VCF format}
\description{
This function is to output variant calls as VCF files.
}
\usage{
    writeVCF(srma.out, reference=NULL, target.file='')
}
\arguments{
	  \item{srma.out}{An SRMAResults object from post.process function.}
	  \item{reference}{Reference to be written into VCF file.}
	  \item{target.file}{A file name of the target VCF file with path.}
	  }
\value{ The output VCF file containing variant calls.
     }
\author{Wenyi Wang \email{wwang7@mdanderson.org}
Nianxiang Zhang \email{nzhang@mdanderson.org}
}
\references{
Wang W, Shen P, Thyagarajan S, Lin S, Palm C, Horvath R, Klopstock T,
Cutler D, Pique L, Schrijver I, Davis RW, Mindrinos M, Speed TP, Scharfe
C. Identification of Rare DNA Variants in Mitochondrial Disorders with
Improved Array-based Sequencing. Nucleic Acids Res. 2011 Jan;39(1):44-58.
}
\keyword{methods}
\keyword{writeVCF}
\examples{
data('reorg.dat')
\dontrun{sd.adj.dat<-single.chip.posterior.prob(reorg.dat)}
data('sd.adj.dat')
srma.out <- call.genotype(reorg.dat,sd.adj.dat )
srma.results <- post.process(srma.out)
writeVCF(srma.results)
}
