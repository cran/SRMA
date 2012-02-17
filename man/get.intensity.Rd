\name{get.intensity}
\alias{get.intensity}
\title{Function for getting raw or quantile normalizated intensities for each base}
\description{
A function is used to get raw or quantile normalized intensities for each base
}
\usage{
get.intensity(chipType,dataDir,chipOrder,ampmap=NA,badamp=NA,exonList=NA,
        chip.normalize=FALSE,amp.normalize=FALSE)
}
\arguments{
  \item{chipType}{chip type name}
  \item{dataDir}{directory for raw data, where CEL files are stored. e.g., dataDir="DataSet", when all CEL files are under "dataSet/chipType".}
  \item{chipOrder}{chip order, a vector of numbers.}
  \item{ampmap}{a data frame, mapping between amplicon and exons. Structure
    of ampmap is: row names are 1 to the total number of amplicons, column
    names are the same as those of the data file {ampmap.rda}. Default is NA, which means that when no quantile normalization is done the value of ampmap will not be used.}
  \item{badamp}{a matrix. Each column is a vector of 0 and 1 for all
    amplicons in each sample. If amplicon is bad, the element is 1,
    otherwise 0. Default is NA, which means that when no
    quantile normalization is done the value of badamp will not be used.}
  \item{exonList}{a vector of numbers which is the list of the exons you
    want to extract from all exons. Default is NA, which provides
    values for the full list of the exons.}
  \item{chip.normalize}{a logic, TRUE or FALSE. It's used to specify
    whether quantile normalization will be done for the raw intensities of different chips. Default is FALSE, which means the raw intensities will be used.}
  \item{amp.normalize}{a logic, TRUE or FALSE. It's used to
    specify whether quantile normalization will be done for raw
    intensities of different amplicons. Default is FALSE, which means the raw intensities will be used.}
}
\value{
A list. Each element of the list contains a matrix of raw intensities
for an exon if the value of "chip.normalize" and "amp.normalize" are
both FALSE, otherwise each element of the list contains a matrix of
quantile normalized intensities for an exon with respect to chips or
amplicons by specifiying "chip.normalize" or "amp.normalize" is TRUE. In all cases, each row of the matrices corresponds to a probe feature and each colume corresponds to a sample.
}
\author{Wenyi Wang \email{wwang7@mdanderson.org}
        Yan Xu  
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
\keyword{quantile normalization}
\examples{
# Input chip type and data directory
  ##chipType <- "MitoP-2r520651"
  ##dataDir <- "dataSet2"

  ##data(badamp)
  ##data(ampmap)
  ##chipOrder <- 1:5
  ##badamp <- badamp[,chipOrder]
  ##exonList <-1:2
  ##chip.normalize <- FALSE
  ##amp.normalize <- FALSE

  ##celint <- get.intensity(chipType=chipType,dataDir=dataDir,
  ##                  chipOrder=chipOrder,ampmap=ampmap,
  ##                  badamp=badamp,exonList=exonList,
  ##                  chip.normalize=chip.normalize,
  ##                  amp.normalize=amp.normalize)

}



