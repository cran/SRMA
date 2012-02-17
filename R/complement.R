complement <- function(x){
    n <- nchar(x)
  if(n==1){
    output <- ifelse(x=="A","T",ifelse(x=="C","G",ifelse(x=="G","C","A")))
  }
  ##if x contains two characters, output the base pairs with no interval
  if(n==2){
    x1 <- substr(x,1,1)
    x2 <- substr(x,2,2)
    output <- paste(ifelse(x1=="A","T",ifelse(x1=="C","G",ifelse(x1=="G","C","A"))),ifelse(x2=="A","T",ifelse(x2=="C","G",ifelse(x2=="G","C","A"))),  sep="")
  }
  return(output)
}
