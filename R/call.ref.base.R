########
## Base caller 1:
##-----------
##a) form the following 4 numbers from the 8 intensities
##
##fA+rT, fC+rG, fG+rC, fT+rA.
##
##b) call base A, C, G or T according to whether the 1st, 2nd, 3rd or 4th of
##these numbers is largest.
##
##Note: this might be orginal scale or log scale. To be decided.
#########

call.ref.base <- function(st,at){
  ##n1 <- logbase2(pri[sample.id,,4])+logbase2(pri[sample.id,,5])
  ##n2 <- logbase2(pri[sample.id,,3])+logbase2(pri[sample.id,,6])
  ##n3 <- logbase2(pri[sample.id,,2])+logbase2(pri[sample.id,,7])
  ##n4 <- logbase2(pri[sample.id,,1])+logbase2(pri[sample.id,,8])

  n1 <- logbase2(st[,1])+logbase2(at[,4])
  n2 <- logbase2(st[,2])+logbase2(at[,3])
  n3 <- logbase2(st[,3])+logbase2(at[,2])
  n4 <- logbase2(st[,4])+logbase2(at[,1])

  bc1 <- apply(cbind(n1,n2,n3,n4),1,call.max)
  return(bc1)
}
