\name{find.indels}
\alias{find.indels}
\title{Function to find insertion/deletions.}
\description{
This function is to determine if there is insertion or deletion present in an exon by examining if a string of hetozygotes in the same sample is present.
}
\usage{
    find.indels(subsnps,code)
}
\arguments{
  \item{subsnps}{A dataframe with genotype calls of an exon. The first 3 columns are the reference information, from column 4 are genotype calls. }
  \item{code}{1 or 2, use single letter or two letters for genotype calls and references.  }
      }
\value{
  A list with 2 components.
\item{indelsnps}{List of SNPs with insertion/deletions.}
\item{singlesnps}{List of SNPs that are not indel snps.}
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
\keyword{indels}
\examples{
# load('allsnps.rda'))
# subsnps <- allsnps[allsnps[,1]==allsnps[1,1],,drop=FALSE]
# out <- find.indels(subsnps=subsnps, code=1)
# str(out)
}