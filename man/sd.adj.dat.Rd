\name{sd.adj.dat}
\docType{data}
\alias{sd.adj.dat}
\title{Normalized dataset }
\description{
An example of output data from single chip analysis. The data is the output of \code{\link{single.chip.posterior.prob}}. 
}
\usage{
sd.adj.dat
}
\format{
  A list with 3 components:
  \enumerate{
  \item{post.prob}{Posterior probabilities. A list with length of 2, corresponding to 2 exons. Each element is an array with 3 dimentions: 1) Nucleotide positions in the exon; 2) Samples; 3) 3 genotypes. }
  \item{allbp.indx}{The index of principle alternative alleles.}
  \item{chipadj.qn}{Single chip adjusted normalized data for principal alternative alleles. A list of 6 elements: 1)f1: fitted value for homozygous reference; 2)f3:fitted value for homozygous alternative; 3)mu1: mean of fitted values across all positions for homozygous reference; 4)mu3: mean of fitted values across all positions for homozygous alternative; 5)errors1: residuals for homozygous reference; 2)errors3: residuals for homozygous alternative;}
}
}
\author{Wenyi Wang \email{wenyiw@stanford.edu}
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
\seealso{
       \code{\link{single.chip.posterior.prob}}.
     }
\keyword{manip}


