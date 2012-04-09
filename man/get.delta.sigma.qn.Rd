\name{get.delta.sigma.qn}
\alias{get.delta.sigma.qn}
\title{Wrapper function for getting quantitle normalized delta and sigma for each base}
\description{
A wrapper function to obtain the log ratios (delta) and log average intensities (sigma) using quantile normalized
intensities for each base}
\usage{
get.delta.sigma.qn(qc.out, ref.info)
}
\arguments{
  \item{qc.out}{Output from quality.control module}
  \item{ref.info}{A dataframe containing exon names(Fragment), position (FragPos) and reference allele (Ref)}.
}
\value{
  a list containing:
  \item{alldelta.qn}{An array with 4 dimensions containing log transformed quantile normalized intensity ratios for each
    base. The first dimension corresponds to the base positions. The second
    dimension corresponds to the three pairs of RM/AM allele, where RM
    stands for reference match and AM stands for alternative match. The third dimension corresponds to samples, and the fourth dimentions corresponds to sense and antisense strands.}
  \item{allsigma.qn}{An array with 4 dimentions containing average of log transformed quantile normalized intensities for each
    RM/AM pairs. The structure of the arrays are the same as element delta.}
  \item{allpairs.qn}{An array with 3 dimensions containing the nucleotide base of RM/AM pairs for all bases. The first dimension corresponds to the base positions. The second
    dimension corresponds to the three pairs of RM/AM allele, and the third dimension corresponds to sense and antisense strands.}
  \item{ref.info}{A dataframe with 5 columns containing amplicon
    names(Fragment), position (FragPos), reference allele (Ref),
    chromosome number, and chromosome coordinate for all the bases.}
  \item{gc}{The gc content of amplicon/exon corresponding to the amplicons in ref.info. }
  \item{exonlength}{The length of amplicon/exon corresponding to the amplicons in ref.info.}
}

\author{Wenyi Wang \email{wwang7@mdanderson.org}
        Nianxiang Zhang   \email{nzhang@mdanderson.org}
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

  ###Do not run
  
  ##chipType <- "MitoP-2r520651"
  ##dataDir <- "dataSet2"
  ##chipOrder <- 1:5
  ##output.measures <- c("D","T")
  ##mode.baselevel <- "average"
  ##mode.amplevel <- "average"
  ##data(ampmap)
  ##exonList <- 1:2
  ##load(ref.info)

  ##QC <- quality.control(chipType=chipType,
  ##                          dataDir=dataDir,
  ##                          chipOrder=chipOrder,
  ##                          output.measures=output.measures,
  ##                          mode.baselevel=mode.baselevel,
  ##                          mode.amplevel=mode.amplevel,
  ##                          ampmap=ampmap,
  ##                          exonList=exonList)  
  ##reorg.dat<-get.delta.sigma.qn(QC, ref.info)
  }




