correct.lr.singlebp <- function(i,alldelta,allsigma,fragl=NULL,allpairs,allgc,allbp.indx,
                                allpost, verbose=F, ...){
##add verbose and plot option
  ##if(Plot && !file.exists('modifiedfigs')) dir.create('modifiedfigs')
  require(nlme) #added
##  good <- allbad[,i]==0
  good <- !is.na(alldelta[,1,i,1])
  indx.s <- cbind(1:sum(good),allbp.indx[good])
  indx.a <- cbind(1:sum(good),4-allbp.indx[good])
  post.ref <- allpost[good,i,,1][indx.s]
  ##sum(post.ref>0.999)/length(post.ref)
  ##[1] 0.954992
  post.het <- allpost[good,i,,2][indx.s]
  post.hom <- allpost[good,i,,3][indx.s]

  ####
  ## now only take the reference base pairs with post>0.999
  ####
  ref.indx <- which(post.ref>0.999)
  
  x <- as.vector(cbind(allsigma[good,,i,1][indx.s][ref.indx],
                       allsigma[good,,i,2][indx.a][ref.indx]))
  y <- as.vector(cbind(alldelta[good,,i,1][indx.s][ref.indx],
                       alldelta[good,,i,2][indx.a][ref.indx]))

  z1 <- as.vector(fragl[good,][ref.indx,])
  z2 <- as.vector(allgc[good,][ref.indx,])
  bp <- as.vector(cbind(allpairs[good,,1][indx.s][ref.indx],
                        allpairs[good,,2][indx.a][ref.indx]))

  dat <- data.frame(y=y,x=x,z1=z1,z2=z2,bp=as.factor(bp))  
  
  if(length(unique(z2))<5) {
    fit0 <- lm(y~I(bp)+ns(x,df=5),data=dat)
    fit1 <- gls(y~I(bp)+ns(x,df=5),data=dat[x>8,])
  } else if(is.null(fragl)) {
    fit0 <- lm(y~I(bp)+ns(x,df=5)+ns(z2,df=3),data=dat)
    fit1 <- gls(y~I(bp)+ns(x,df=5)+ns(z2,df=3),data=dat[x>8,])
  } else {
    fit0 <- lm(y~I(bp)+ns(x,df=5)+ns(z1,df=3)+ns(z2,df=3),data=dat)
    fit1 <- gls(y~I(bp)+ns(x,df=5)+ns(z1,df=3)+ns(z2,df=3),data=dat[x>8,])
  }  
  fit2 <- try(update(fit1,weights=varConstPower(form=~fitted(.))))
  if(class(fit2)=='try-error') fit2 <- fit1 #added by NZ
  
#######moved out from plot chunk  
  a <- fitted(fit0)
  b <- studres(fit0)^(2/3)
  tmp1 <- a>0&!is.na(b)
  tmp2 <- lowess(x=log(a[tmp1]),y=b[tmp1],f=1/3)
  
  if(verbose){ 
    sd1 <- (min(tmp2$y))^(3/2)
    sd2 <- (tmp2$y[order(tmp2$x)][length(tmp2$y)])^(3/2)
    cat('Identify the range of heteroscadasticity: \n sd1=', sd1, '\n sd2=',sd2, '\n sd2/sd1=', sd2/sd1, '\n Call for weighting when sd2/sd1>1.5 \n')
    ##[1] 2.43495
    ##this ratio is greater than 1.5, therefore calls for weighting
  }

  ##Now make prediction for all good bases, some were not used in the fitting because of their low posterior probabilities for being the reference
  x.new <- as.vector(cbind(allsigma[good,,i,1][indx.s],
                           allsigma[good,,i,2][indx.a]))
  y.new <- as.vector(cbind(alldelta[good,,i,1][indx.s],
                           alldelta[good,,i,2][indx.a]))
  
  z1.new <- as.vector(fragl[good,])
  z2.new <- as.vector(allgc[good,])
  bp.new <- as.vector(cbind(allpairs[good,,1][indx.s],
                            allpairs[good,,2][indx.a]))
  
  newdat1 <- data.frame(y=y.new,x=x.new,z1=z1.new,z2=z2.new,bp=as.factor(bp.new))
  f1 <- predict(fit2,newdata=newdat1)
  mu1 <- median(f1)
  ## power of the mean, theta is the power
  if(verbose) print(fit2)
  const <- unlist(fit2$modelStruct$varStruct)[1]
  if(verbose) print(const)
  theta <- unlist(fit2$modelStruct$varStruct)[2]
  error1 <- summary(fit2)$sigma*(exp(const)+abs(f1)^theta)
 # #################
 # bp.new2 <- sapply(bp.new,function(a){
 #   a1 <- unlist(strsplit(a,split=""))
 #   a2 <- paste(a1[2],a1[1],sep="")
 #   return(a2)
 # })
 ########### Being replaced by
  bp.new2<-paste(substr(bp.new,2,2), substr(bp.new,1,1), sep="")
  
  newdat2 <- data.frame(y=y.new,x=x.new,z1=z1.new,z2=z2.new,bp=as.factor(bp.new2))
  f3 <- -predict(fit2,newdata=newdat2)
  mu3 <- median(f3)
  error3 <- summary(fit2)$sigma*(exp(const)+abs(f3)^theta)
  
  output <- list(f1,f3,mu1,mu3,error1,error3)
  names(output) <- c("f1","f3","mu1","mu3","error1","error3")
  return(output)
}
