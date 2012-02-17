call.max <- function(y){
 call <- c("T","G","C","A")
 x <- which(y==max(y))
 z <- length(x)
 if(z>1){
   x <- sample(x,1)
 }
 return(call[x])
}
