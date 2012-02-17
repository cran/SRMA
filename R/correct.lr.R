correct.lr <- function(i,alldelta.qn,allsigma.qn,fragl=NULL,allpairs.qn,allgc,verbose=FALSE){
  
  if(verbose) cat('Processing sample', dimnames(alldelta.qn)[[3]][i], '... \n')
  
  x <- as.vector(allsigma.qn[!is.na(alldelta.qn[,1,i,1]),,i,]) #average log intensity
  y <- as.vector(alldelta.qn[!is.na(alldelta.qn[,1,i,1]),,i,])  #logratios
  if(is.null(fragl)) z1<-0 else z1 <- fragl[as.vector(!is.na(alldelta.qn[,1,i,1]))]
  z2 <- allgc[as.vector(!is.na(alldelta.qn[,1,i,1]))]
  bp <- as.vector(allpairs.qn[!is.na(alldelta.qn[,1,i,1]),,])
  dat <- data.frame(y=y,x=x,z1=z1,z2=z2,bp=as.factor(bp))
  new.bp<-paste(substr(bp,2,2), substr(bp,1,1), sep="")
  new.dat<-cbind(dat[,1:4], bp=new.bp)
  require(splines)
  
  if(nrow(dat)>100000) {
    if(verbose) cat(paste('The data have',nrow(dat) , 'observations, 100,000 random samples will be used to build the model and correct the data. \n'))
    samp.ind<-sample(1:nrow(dat), 100000)
  }else samp.ind<-1:nrow(dat)
  mod.dat<-dat[samp.ind,]
  if(verbose) cat('Fitting model. \n')
  if(length(unique(z2))<5) {
    fit1<-lm(y~ns(x,df=5)+I(bp),data=mod.dat)  #this is added to make sure the test data pass 7/26/11
  } else if(is.null(fragl)) {
    fit1 <- lm(y~ns(x,df=5)+ns(z2,df=3)+I(bp),data=mod.dat)
  } else {
    fit1 <- lm(y~ns(x,df=5)+ns(z1,df=3)+ns(z2,df=3)+I(bp),data=mod.dat) 
  }
  if(verbose) cat('Predicting. \n')
  f1<-predict(fit1, newdata=dat)
  f3<- -predict(fit1, newdata=new.dat)
  mu1 <- mean(f1)
  mu3 <- mean(f3)
  error <- rep(summary(fit1)$sigma,length(y))
  output <- list(f1,f3,mu1,mu3,error)
  names(output) <- c("f1","f3","mu1","mu3","error")
  if(verbose) cat('Done! \n')
  return(output)
}

