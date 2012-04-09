writeVCF<-function(srma.out, reference=NULL, target.file=''){
    if(class(srma.out)!='SRMAResults') stop('An object of SRMAResults class is required.')
    ##header
    cat('##fileformat=VCFv4.0\n', target.file, append=F, sep='')
    cat(paste('##fileDate=', as.character(as.Date(Sys.time()), format='%Y%m%d') , '\n',sep=''), target.file, append=T, sep='')
    cat('##source=SRMA\n', target.file, append=T, sep='')
    if(!is.null(reference)) cat(    paste('##reference=', reference, sep=''), target.file, append=T, sep='')
    cat( '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', target.file, append=T, sep='')
    ##body
    snps<-srma.out@srmaSnps
    stopifnot(all((c('Chromosome','Position') %in% colnames(snps)))) 
    if(nrow(snps)==0) break
    snps<-snps[!is.na(snps$Chromosome),]
    b1<-snps$Chromosome
    b2<-snps$Position
    b3<-paste(snps$Fragment, snps$FragPos, sep='_')
    b4<-toupper(as.character(snps$Ref))
    sn<-srma.out@sampleNames
    gt.temp<-as.matrix(snps[,sn])
    
    Ref<-b4
    alt<-gt<-NULL
    for(i in 1:nrow(gt.temp)){
      pos.i<-gt.temp[i,]
      ref<-Ref[i]
      allele1<-ifelse(substr(pos.i,1,1)==ref, 0, 1)
      allele1[is.na(allele1)]<-'.'
      allele2<-ifelse(substr(pos.i,2,2)==ref, 0, 1)
      allele2[is.na(allele2)]<-'.'
      gt.i<-paste(allele1, allele2, sep='/')
      allele.all<-na.exclude(unique(c(substr(pos.i,1,1),substr(pos.i,2,2) )))
      alt.i<-allele.all[allele.all!=ref]
      gt<-rbind(gt, gt.i)
      alt<-c(alt, alt.i)
    }
    colnames(gt)<-colnames(gt.temp)
    b5<-alt
    b6<-b7<-b8<-rep('.', nrow(snps))
    b9<-rep('GT', nrow(snps))
    mand<-cbind(b1,b2,b3,b4,b5,b6,b7,b8,b9)
    colnames(mand)<-c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
    vcfbody<-cbind(mand, gt)
    cat(colnames(vcfbody), file=target.file, append=T, sep='\t')
    cat('\n', file=target.file, append=T, sep='')
    write.table(vcfbody, file=target.file, append=T, quote=FALSE, sep='\t', na='.', row.names=FALSE, col.names=FALSE)
  }
