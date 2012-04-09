\name{post.process}
\alias{post.process}
\title{Wrapper function for output analysis}
\description{
This function is the wrapper function for analyzing the output of SRMA algorithms. In the analysis, footprint effects will be checked; the probe quality will be accessed; rare SNPs and indels will be identified.
}
\usage{
    post.process(srma.output, dbsnp.info=NULL, sls.cut=0.67, ampmap=NULL)
}
\arguments{
	  \item{srma.output}{Direct output of SRMA, a list of :
	  
	  alldelta.qn:array of normalized logratios for samples with 4 dimensions: positions,bases(3), samples, strands(2). 
	  
	  allsigma.qn:array of normalized log averages for samples with 4 dimensions: positions,bases(3), samples, strands(2).
	  
	  allpairs.qn:array of genotypes with 3 dimensions: positions,bases(3), strands(2).
	  
	  ref.info:Information of reference call for all the positions. A dataframe with 3 elements: Fragment(amplicon name), FragPos(Position on the exon), Ref(Reference genotye of the position).
	  
	  gc:The GC content of the amplicons.
	  
	  exonlength:The length of the amplicons.
	  }
	  \item{dbsnp.info}{A dataframe with information of dbSNPs containing column names
	  of 'ExonID', 'BasePos','Ref', 'Alt.Allele', 'rs'.}
	  
	  \item{sls.cut}{Cutoff value for Silhoutte scores to filter out
	  unreliable calls}
	  
	  \item{ampmap}{A dataframe mapping exons to amplicons with column names
	  of 'exon','exon_len','exon_start', 'exon_end', 'amplicon','amp_len',  'amp_start' , 'amp_end'.}
	  }
\value{An SRMAResults object containing a list of 6 elements:
 \item{footprint}{A list of number of footprint and number of samples with footprint.}
	  \item{srmaMat}{A list of matrices (position in row, sample in column) recording SRMA analysis results. (1) base quality (good:T, bad=F) (2) genotype calls (3) quality scores (4) Silhoutte scores.}
	  \item{srmaIndelsnps}{A dataframe with identified insertion/deletion calls.}
	  	  \item{srmaSnps}{A dataframe with identified SNV calls.}
	  	  \item{sampleNames}{A vector of sample names.}
	  	  \item{slscoreCutoff}{A cutoff value for Silhoutte
		    score to filter out unreliable calls.}
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
\keyword{post.process}
\examples{
data('reorg.dat')
\dontrun{sd.adj.dat<-single.chip.posterior.prob(reorg.dat)}
data('sd.adj.dat')
srma.out<-call.genotype(reorg.dat,sd.adj.dat )
Results<-post.process(srma.out)
Results
writeVCF(Results)                               
}