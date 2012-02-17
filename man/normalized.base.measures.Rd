\name{normalized.base.measures}
\alias{normalized.base.measures}
\title{A high-level function for getting delta and sigma for each base}
\description{
A function is used to get delta and sigma using quantile normalized probe intensities for each base with bad amplicons removed
}
\usage{
normalized.base.measures(chipType,dataDir,chipOrder,ampmap,badamp,exonList,
                         chip.normalize=FALSE,amp.normalize=TRUE)
}
\arguments{
  \item{chipType}{chip type name}
  \item{dataDir}{directory for raw data, where CEL files are stored. e.g., dataDir="DataSet", when all CEL files are under "dataSet/chipType".}
  \item{chipOrder}{chip order, a vector of numbers.}
  \item{ampmap}{a data frame, mapping between amplicon and exons. Structure
    of ampmap is: row names are 1 to the total number of amplicons, column
    names are the same as those of the data file {ampmap}. Default is
    NA, which means that when no quantile normalization is done the
    value of ampmap will not be used.}
   \item{badamp}{a matrix. Each column is a vector of 0 and 1 for all amplicons in each sample. If amplicon is bad, the element is 1, otherwise 0.}
   \item{exonList}{a vector of numbers which is the list of the exons you
    want to extract from all exons. Default is NA, which provides
    values for the full list of the exons.}
   \item{chip.normalize}{a logic, TRUE or FALSE. It's used to specify
    whether quantile normalization will be done for the raw intensities of different chips. Default is FALSE, which means the raw intensities will be used.}
  \item{amp.normalize}{a logic, TRUE or FALSE. It's used to
    specify whether quantile normalization will be done for raw
    intensities of different amplicons. Default is TRUE, which means the
    quantile normalized intensities will be used.}
}
\value{
 a list of 5 containing:
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
\keyword{amplicon}
\keyword{exon}
\examples{
 ## chipType <- "MitoP-2r520651"
 ## dataDir <- "dataSet2"
 ## chipOrder <- 1:5
 ## load the ampmap data for 5 samples
 ## data(ampmap)
 ## exonList <- 1:2
 ## badamp <- quality.control(chipType=chipType,
 ##                           dataDir=dataDir,
 ##                           chipOrder=chipOrder,
 ##                           mode=mode,
 ##                           ampmap=ampmap,
 ##                           exonList=exonList)

 ##delta.sigma.qn <-normalized.base.measures(chipType=chipType,
 ##                                          dataDir=dataDir,
 ##                                          chipOrder=chipOrder,
 ##                                          ampmap=ampmap,
 ##                                          badamp=badamp,
 ##                                          exonList=exonList)
 ##str(delta.sigma.qn[[2]])
}



