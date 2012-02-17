\name{quality.control}
\alias{quality.control}
\title{A high-level function for identifying bad amplicons using raw intensities from CEL files}
\description{
A function used to identify bad amplicons with the criterion: the reference call rate R<0.9, using raw intensities from CEL files
}
\usage{
quality.control(chipType,dataDir,chipOrder,output.measures=c("R","D","T"),
                mode.baselevel,mode.amplevel,ampmap=NA,exonList=NA)
}
\arguments{
  \item{chipType}{chip type name}
  \item{dataDir}{directory for raw data, where CEL files are stored. e.g., dataDir="DataSet", when all CEL files are under "dataSet/chipType".}
  \item{chipOrder}{chip order, a vector of numbers.}
  \item{output.measures}{the output selection of QC
    measures. "R" represents the reference call rate. "D" represents the median of dm. "T" represents MAD of am. Default is "R","D", and"T". You can choose to input any combination of the three measures.}
  \item{mode.baselevel}{a character, "refcallrate" or
    "average". Default is "refcallrate". Mode "refcallrate" means the
    reference call rate is calculated by a simple base caller for each sample. Mode "average" means the reference call rate is 1 if the RM probe intensities is highest, otherwise 0.}
  \item{mode.amplevel}{a character, "average" or "maximum". Default is "average". Mode "average" means the reference call rate is calculated by averaging values for sense and antisense strand. Mode "maximum" means the reference call rate is calculated by taking the maximum of values for sense and antisense strand.}
  \item{ampmap}{a data frame, mapping between amplicon and exons. Structure
    of ampmap is: row names are 1 to the total number of amplicons, column
    names are the same as those of the data file ampmap.rda. Default is NA, which means that when no quantile normalization is done the value of ampmap will not be used.}
   \item{exonList}{a vector of numbers which is the list of the exons you
    want to extract from all exons. Default is NA, which provides
    values for the full list of the exons.}
}
\value{
A list. The first element of the list contains a matrix. Each column is
a vector of 0 or 1 for all amplcons in each sample. Element 1 refers to
a bad amplicon, 0 refers to a good one. Additional element(s) in the
middle of the
list is chosen by the selection of the argument output.measures. The
default of output.measures is "R","D",and "T". "R","D" and "T" have the
same structure: An array with dimensions (1+total number of
amplicons)x2x(number of samples). The last element is the parameters
from the function arguments including chipType,dataDir,chipOrder, and ampmap.
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
\keyword{reference call rate}
\keyword{quality control}
\examples{
  ##chipType <- "MitoP-2r520651"
  ##dataDir <- "dataSet2"
  ##chipOrder <- 1:5
  ##output.measures <- c("D","T")
  ##mode.baselevel <- "average"
  ##mode.amplevel <- "average"
  ##load the ampmap data for 5 samples
  ##data(ampmap)
  ##exonList <- 1:2

  ##output <- quality.control(chipType=chipType,
  ##                          dataDir=dataDir,
  ##                          chipOrder=chipOrder,
  ##                          output.measures=output.measures,
  ##                          mode.baselevel=mode.baselevel,
  ##                          mode.amplevel=mode.amplevel,
  ##                          ampmap=ampmap,
  ##                          exonList=exonList)
  ##is.list(output)
  ##str(output)
}

