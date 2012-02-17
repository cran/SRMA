#
#Author: Nianxiang Zhang UTMDACC Oct 17 2011
#
library(methods)
setClass("SRMAResults", representation(
				footprint="list", 
				srmaMat="list", 
				srmaIndelsnps="data.frame",
				srmaSnps="data.frame",
				sampleNames='vector',
				slscoreCutoff='numeric'
				))

setMethod("initialize", "SRMAResults", 
		function(.Object,footprint,srmaMat,srmaIndelsnps,srmaSnps,sampleNames, slscoreCutoff ) 
		{
			.Object@footprint <- footprint
			.Object@srmaMat <- srmaMat 
			.Object@srmaIndelsnps <- srmaIndelsnps
			.Object@srmaSnps <- srmaSnps
			.Object@sampleNames <- sampleNames
			.Object@slscoreCutoff<-slscoreCutoff
			.Object
		})

setMethod("show",    signature(object = "SRMAResults"),
		function(object) 
		{.local <- function (x){
		cat('Analysis Results From SRMA Analysis Pipeline: \n-----------------------------------------------------\n-----------------------------------------------------\n QC: \n')
			qc.m<-x@srmaMat$srma.qc
			qc.pct<-sum(!qc.m)/(nrow(qc.m)*ncol(qc.m)) *100
			s.qc<-cbind(nrow(qc.m)*ncol(qc.m),sum(!qc.m),qc.pct)
			colnames(s.qc)<-c('Total', 'QC.Failed', 'Percent.Failed')
			rownames(s.qc)<-'Base Number'
			print(s.qc)
			cat('-----------------------------------------------------\n Footprint Effects: \n     Bases:', x@footprint$base.n, '\n     Samples:', x@footprint$sample.n, '\n     ')
			ft.pct<-x@footprint$base.n/(nrow(qc.m)*ncol(qc.m)-sum(is.na(qc.m)))*100
			cat(round(ft.pct,4), 'percent of bases are affected by footprint effects. \n-----------------------------------------------------\n Silhouette Scores: \n')
			slscore<-x@srmaMat$srma.slscore
			n.low<-sum(slscore<x@slscoreCutoff,na.rm=T)
			n.all<-(nrow(qc.m)*ncol(qc.m)-sum(!qc.m))
			perc.low<-sum(slscore<x@slscoreCutoff,na.rm=T)/(nrow(qc.m)*ncol(qc.m)-sum(!qc.m))
			sls.t<-cbind(n.low, n.all, perc.low)
			colnames(sls.t)<-c('Low', 'Total', 'Percent')	
						rownames(sls.t)<-'Base Number'

			cat('     Summary of Silhouette scores (cutoff=', x@slscoreCutoff, ').\n') 
			print(sls.t)
			allsnps3<-x@srmaSnps
			cat('-----------------------------------------------------\n Single nucleotide variations: \n     ', nrow(allsnps3), 'SNPs are identified. \n')
			if(nrow(allsnps3)!=0) {
			cat('     ',sum(allsnps3[,4:ncol(allsnps3)]!='N', na.rm=T), 'SNPs are called in', ncol(qc.m), 'samples. \n     ')
			snp.pct<-sum(allsnps3[,4:ncol(allsnps3)]!='N', na.rm=T)/sum(qc.m)
			cat(round(snp.pct*100, 4), 'percent of calls are SNPs. \n')
			}
			indelsnps<-x@srmaIndelsnps
			cat('-----------------------------------------------------\n Insertion/Deletions: \n     ', dim(indelsnps)[1], 'insertion/deletions are identified. \n     ')
			cat('     Insertion/deletions are identified in', length(unique(indelsnps[,1])), 'exons. \n') 
			}
		.local(object)
		})

	
setMethod("summary",    signature(object = "SRMAResults"),
		function(object) 
		{.local <- function (x){
		cat('Analysis Results From SRMA Analysis Pipline: \n-----------------------------------------------------\n-----------------------------------------------------\n QC: \n')
			qc.m<-x@srmaMat$srma.qc
			qc.pct<-sum(!qc.m)/(nrow(qc.m)*ncol(qc.m)) *100
			s.qc<-cbind(nrow(qc.m)*ncol(qc.m),sum(!qc.m),qc.pct)
			colnames(s.qc)<-c('Total', 'QC.Failed', 'Percent.Failed')
			rownames(s.qc)<-'Base Number'
			print(s.qc)
			cat('-----------------------------------------------------\n Footprint Effects: \n     Bases:', x@footprint$base.n, '\n     Samples:', x@footprint$sample.n, '\n     ')
			ft.pct<-x@footprint$base.n/(nrow(qc.m)*ncol(qc.m)-sum(is.na(qc.m)))*100
			cat(round(ft.pct,4), 'percent of bases are affected by footprint effects. \n-----------------------------------------------------\n Silhouette Scores: \n')
			slscore<-x@srmaMat$srma.slscore
			n.low<-sum(slscore<x@slscoreCutoff,na.rm=T)
			n.all<-(nrow(qc.m)*ncol(qc.m)-sum(!qc.m))
			perc.low<-sum(slscore<x@slscoreCutoff,na.rm=T)/(nrow(qc.m)*ncol(qc.m)-sum(!qc.m))
			sls.t<-cbind(n.low, n.all, perc.low)
			colnames(sls.t)<-c('Low', 'Total', 'Percent')	
						rownames(sls.t)<-'Base Number'

			cat('     Summary of Silhouette scores (cutoff=', x@slscoreCutoff, ').\n') 
			print(sls.t)
			allsnps3<-x@srmaSnps
			cat('-----------------------------------------------------\n Single nucleotide variations: \n     ', nrow(allsnps3), 'SNPs are identified. \n')
			if(nrow(allsnps3)!=0) {
			cat('     ',sum(allsnps3[,4:ncol(allsnps3)]!='N', na.rm=T), 'SNPs are called in', ncol(qc.m), 'samples. \n     ')
			snp.pct<-sum(allsnps3[,4:ncol(allsnps3)]!='N', na.rm=T)/sum(qc.m)
			cat(round(snp.pct*100, 4), 'percent of calls are SNPs. \n')
			}
			indelsnps<-x@srmaIndelsnps
			cat('-----------------------------------------------------\n Insertion/Deletions: \n     ', dim(indelsnps)[1], 'insertion/deletions are identified. \n     ')
			cat('     Insertion/deletions are identified in', length(unique(indelsnps[,1])), 'exons. \n') 
			}
		.local(object)
		})

		