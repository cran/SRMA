\name{all.qn}
\docType{data}
\alias{all.qn}
\title{Normalized exon level data}
\description{
 a list containing:
 
  1) celnames: sample names of 5 hapmap samples: NA11840, NA12156, NA12878, NA18507 and NA19625.
    
  2) all.logPM.qn: log intensity for perfect match (or reference match) probe only after quantile normalization. The structure is a list of
    matrices. Each matrix corresponds to an exon and one strand, st-sense strand, at-antisense strand. Each row corresponds to a base position, and each column corresponds to a sample.
  
  3) all.delta.qn: log ratios for each base position after normalization. The structure is a list of arrays with three dimensions. The first dimension corresponds to the base positions. The second dimension corresponds to the three pairs of RM/AM allele, where RM stands for reference match and AM stands for alternative match. And the third dimension corresponds to samples.
^M
    4) all.sigma.qn: log average intensities for each base after normalization. The structure is the same as the argument \code{all.delta.qn}.
  
    5) all.pairs.qn: base pairs for all bases. The structure is a list of matrices. Each row of the matrices corresponds to a base position,
    and the three columns correspond to the three pairs of RM/AM allele,
    where RM stands for reference match and AM stands for alternative match.
}
\usage{all.qn}
\format{
A list 
  
  1) A vector with 5 sample names

  2) log normalized intensities: list of 6 matrices corresponding to 3 exons with sense and antisense

  3) log normalized intensity ratios: list of 6 arrays corresponding to 3 exons with sense and antisense

  4) log normalized intensity average: list of 6 arrays corresponding to 3 exons with sense and antisense

  5) base pairs: list of 6 matrices corresponding to 3 exons with sense and antisense
}
\author{Wenyi Wang \email{wwang7@mdanderson.org}
        Nianxiang Zhang  \email{nzhang@mdanderson.org}}
\references{
Wang W, Shen P, Thyagarajan S, Lin S, Palm C, Horvath R, Klopstock T,
Cutler D, Pique L, Schrijver I, Davis RW, Mindrinos M, Speed TP, Scharfe
C. Identification of Rare DNA Variants in Mitochondrial Disorders with
Improved Array-based Sequencing. Nucleic Acids Res. 2011 Jan;39(1):44-58.

Shen P, Wang W, Krishnakumar S, Chi AK, Palm C, Enns GM, Davis RW, Speed TP, 
Mindrinos M, Scharfe C. High quality DNA sequence capture of 524 disease candidate genes. 
Proc Natl Acad Sci U S A. 2011 Apr 19;108(16):6549-54.
}
\keyword{datasets}
\examples{
data(all.qn)
str(all.qn)
}

