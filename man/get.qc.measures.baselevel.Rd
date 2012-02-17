\name{get.qc.measures.baselevel}
\alias{get.qc.measures.baselevel}
\title{Function for preparing the measures used for quality assessment of DNA targets and resequencing arrays}
\description{
A function is used to calculate, at each base position in each sample,
the reference match (RM) probe intensities; the reference call
indicator; difference between the RM probe intensities and the maximum
of the alternative match (AM) probe intensities; average of RM and AM
intensities; and the difference between the RM intensities and the median of RM intensities for each probe feature. All these measures are used for quality assessment of DNA targets and resequencing arrays.
}
\usage{
get.qc.measures.baselevel(celint,chipType,exonList=NA,celnames,
                          mode.baselevel="refcallrate")}
\arguments{
  \item{celint}{raw or quantile normalized intensities for each probe. A
    list. Each element of the list contains a matrix of intensities
    for an exon, in which each row of the matrices corresponds to a
    probe feature and each column correspons to a sample.}
  \item{chipType}{chip Type name}
  \item{exonList}{a vector of scalers which is the list of the exons you
    want to extract from all exons. Default is NA, wich provides values
    for the full list of the exons.}
  \item{celnames}{sample name, a vector of character strings}
  \item{mode.baselevel}{choices of mode.baselevel are "refcallrate" and
    "average". Default is "refcallrate". Mode "refcallrate" means the
    reference call rate is calculated by a simple base caller for each sample. Mode "average" means the reference call rate is 1 if the RM probe intensities is highest, otherwise 0.}
}
\value{
  a list of 5 containing:
  \item{all.logPM}{RM probe intensities. A list of matrices. Each matrix corresponds to an exon. Each row of the matrices corresponds to a base position and each column corresponds to a sample.}
  \item{all.refc}{reference call indicators of 1 or 0. If the mode is
    "refcallrate", the reference call indicator is obtained by a simple
    base caller for each sample. If the mode is "average" the reference
    indicator is 1 when the RM probe intensities is highest, otherwise 0. The structure is the same as {all.logPM}.}
  \item{all.dm}{the difference between the RM log intensities and the
    maximum of the AM log intensities. The structure is the same as {all.logPM}.}
  \item{all.am}{average of the RM log intensities and the maximum of the
    AM log intensites. The structure is the same as {all.logPM}.}
  \item{all.rle}{the difference between the RM log intensities and the
  median of the RM log intensities across all samples for each probe feature. The structure is the same as {all.logPM}.}
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
 ## load the intensity data
    data(celint)
 ## chipType <- "MitoP-2r520651"
 ## exonList <- 1:2
 ## celnames <- colnames(celint[[1]])
 ## mode.baselevel <- "refcallrate"
 ##y <- get.qc.measures.baselevel(celint=celint,
 ##                               chipType=chipType,
 ##                               exonList=exonList,
 ##                               celnames=celnames,
 ##                               mode.baselevel=mode.baselevel)
 ## str(y[[2]])
}


