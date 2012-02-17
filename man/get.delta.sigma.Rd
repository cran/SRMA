\name{get.delta.sigma}
\alias{get.delta.sigma}
\title{Function for getting delta and sigma for each base}
\description{
A function is used to get the delta and sigma using quantile normalized
intensities for each base
}
\usage{
get.delta.sigma(chipType,dataDir,chipOrder,celint,exonList=NA)
}
\arguments{
  \item{chipType}{chip type name}
  \item{dataDir}{directory for raw data, where CEL files are stored}
  \item{chipOrder}{chip order, a vector of numbers.}
  \item{celint}{raw or quantile normalized intensity matrices for each amplicon obtained by function get.intensity}
  \item{exonList}{a vector of numbers which is the list of the exons you want to extract from all exons. Default is NA, which provides values for the full list of the exons.}
}
\value{
  a list containing:
  \item{celnames}{sample names, a list of character strings.}
  \item{all.logPM.qn}{logarithm of quantile normalized intensity
    matrices for all amplicons. The structure is a list of
    matrices. Each matrix corresponds to one amplicon. Each row refers
    to one base position, and each column refers to a sample.}
  \item{all.delta.qn}{log quantile normalized intensity ratios for each
    base. The structure is a list of arrays of three dimensions. The
    first dimension corresponds to the base positions. The second
    dimension corresponds to the three pairs of RM/AM allele, where RM
    stands for reference match and AM stands for alternative match. And the third dimension corresponds to samples.}
  \item{all.sigma.qn}{log quantile normalized average intensities for each
    base. The structure is the same as the argument {all.delta.qn}.}
  \item{all.pairs.qn}{base pairs for all bases. Its structure is a list
    of matrices. Each row of the matrices corresponds to a base position,
    and the three columns correspond to the three pairs of RM/AM allele,
    where RM stands for reference match and AM stands for alternative match.}
}

\author{Wenyi Wang \email{wwang7@mdanderson.org}
        Yan Xu  
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
\keyword{delta}
\keyword{sigma}
\examples{
# load R packages if not yet loaded
  require(R.utils)
  require(R.filesets)
  require(aroma.core)
  require(aroma.affymetrix)

  data(ampmap)
  data(badamp)

  ##chipType <- "MitoP-2r520651"
  ##dataDir <- "dataSet2"
  ##chipOrder <- 1:5
  ##badamp <- badamp[,chipOrder]
  ##exonList <- 1:2
  ##chip.normalize <- TRUE
  ##amp.normalize <- TRUE

  ##celint <- get.intensity(chipType=chipType,dataDir=dataDir,
  ##                  chipOrder=chipOrder,ampmap=ampmap,
  ##                  badamp=badamp,exonList=exonList,
  ##                  chip.normalize=chip.normalize,
  ##                  amp.normalize=amp.normalize)
  ##data(celint)
  ##y <- get.delta.sigma(chipType=chipType,dataDir=dataDir,
  ##                        chipOrder=chipOrder,celint=celint,
  ##                        exonList=exonList)
  ##str(y)
}




