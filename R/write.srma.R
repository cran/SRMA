write.srma<-function(x, file='SrmaResults',...) {
  stopifnot(class(x)=="SRMAResults")
  snp<-x@srmaSnps    
  Indel<-x@srmaIndelsnps
  snp.f<-ind.f<-temp.f<-unlist(strsplit(file, split='/'))
  fn<-temp.f[length(temp.f)]
  snp.fn<-paste('SNP', fn, sep='-')
  ind.fn<-paste('Indel', fn, sep='-')
  snp.f[length(temp.f)]<-snp.fn
  ind.f[length(temp.f)]<-ind.fn
  write.table(snp, file=paste(snp.f, collapse='/'), ...)
  write.table(Indel, file=paste(ind.f, collapse='/'), ...)
}
