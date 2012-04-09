\name{explore.qc}
\alias{explore.qc}
\title{Function for plotting qc measures for individual sample}
\description{
A function is used to plot quality control measures for individual sample}
\usage{explore.qc(qc,sample.ID,fig.dir,fig.name)}
\arguments{
  \item{qc}{A list which is the output of the function
    quality.control. The first element of the list contains a
    matrix for identifying bad amplicons. Each column is a vector of 0
    or 1 for all amplcons in each sample. Element 1 refers to a bad
    amplicon, 0 refers to a good one. Additional element(s) of the list
    is chosen by the selection of the argument output.measures. This qc
    plot requires the argument output.measures of the function
    quality.control includes "D" and "T". "R" means the reference call
    rate. "D" means median of dm. "T" means median of am. They have the
    same structure: An array with dimensions (1+total number of amplicons)x2x(number of samples).}
  \item{sample.ID}{A scaler, the ID number of the sample to be plotted.}
  \item{fig.dir}{the directory to write the figure}
  \item{fig.name}{the name of the figure, ending with ".eps" or ".ps".}
}
\value{
A plot "D"~"T", where D means median of dm (the difference between the RM
log intensities and the maximum of the AM log intensities) and T means
median of am (average of the RM log intensities and the maximum of the AM log intensites)
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
\keyword{reference call rate}
\keyword{plot}
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

  ##qc <- output
  ##sample.ID <- 2
  ##fig.name <- "explore.qc.eps"

  ##explore.qc(qc=qc,sample.ID=sample.ID,
  ##           fig.dir=fig.dir,fig.name=fig.name)

}


