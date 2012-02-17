\name{amp.spec.qn}
\alias{amp.spec.qn}
\title{Quantile normalization function on amplicons}
\description{
A function which applys quantile normalization across all samples on an
exon if the exon corresponds to one amplicon, with bad exons removed; if
the exon corresponds to multiple amplicons, translate the exon into the
amplicon structure, normalize on the amplicon.
}
\usage{
      amp.spec.qn(exon.i,x,badamp,ampmap)
}
\arguments{
  \item{exon.i}{exon name}
  \item{x}{probe level intensity, a matrix whose structure is: for n
    probe positions, row 1:n for probe1,row n+1:2n for probe2,row
    2n+1:3n for probe3, row 3n+1:4n for probe4, and columns correspond to sample names}
  \item{badamp}{a matrix. Each column is a vector of 0 and 1 for all amplicons in each sample. If amplicon is bad, the element is 1, otherwise 0.}
  \item{ampmap}{a data frame, mapping between amplicons and exons. Structure of map is: row names are 1 to the total number of amplicons, column names are the same as the those of the data file {ampmap.rda}.}
}
\value{
A matrix of values which are quantile normalized intensities on amplicon
level. Each row of the matrix corresponds to a probe feature and each colume corresponds to a sample.
}
\author{Wenyi Wang \email{wwang7@mdanderson.org},
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
\keyword{quantile normalization}
\keyword{amplicon}
\keyword{exon}
\seealso{\code{ampmap}}
\examples{
# load the ampmap and badamp data
#  data(ampmap)
#  data(badamp)

# Input chip type and data directory
 ## chipType <- "MitoP-2r520651"
 ## dataDir <- "dataSet2"

# read in raw cel files
  ##cdf <- AffymetrixCdfFile$byChipType(chipType,tags="byExon");
  ##csR <- AffymetrixCelSet$byName(dataDir, cdf=cdf);
  ## units <- 10:nbrOfUnits(cdf);
  ##if(chipType=="MitoP-2r520651"){
  ##   units <- 10:nbrOfUnits(cdf);
  ##   cellQuartets <- getCellQuartets(cdf,units=units);
  ##}else{
  ##  cellQuartets <- getCellQuartets(cdf);
  ##}

  ##unitnames <- names(cellQuartets)
  ##exons <- unique(c("AFFX-TagIQ-EX",
  ##                 substr(unitnames[regexpr("ENSE",unitnames)>0],1,15)))
  ##if(sum(is.na(exonList))==0){
  ##    subexons <- exons[exonList]
  ##}else{
  ##    subexons <- exons
  ##}
  ##exon.i<-sample(subexons,1)

# read in probe intensities
  ##unitindx <- which(regexpr(i,unitnames)>0)

# Get both sense and antisense data
# Use function extractMatrix to exract the pooled
# intensities from csR
  ##cells.st <- t(cellQuartets[[unitindx[1]]][[1]])
  ##x.st <- extractMatrix(csR,cells=cells.st)
  ##cells.at <- t(cellQuartets[[unitindx[2]]][[1]])
  ##x.at <- extractMatrix(csR,cells=cells.at)

  ##chipOrder<-1:5
  ##x <- rbind(x.st,x.at)[,chipOrder]

  ##badamp<-badamp[,chipOrder]
  ##y <- amp.spec.qn(exon.i,x,badamp,ampmap=ampmap)
  ##str(y)
}


