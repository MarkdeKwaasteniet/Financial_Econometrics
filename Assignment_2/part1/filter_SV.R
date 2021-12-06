#This R code contains the function that needs to be minimized,
#with respect to ft at each time period to obtain an approximation for
#the filtered volatility of the SV model. 
#For a more detailed explanation see the Lecture Notes.

filter_SV <- function(yt,ft,ft1,theta){
  
  omega <- theta[1]   
  beta <- theta[2]   
  sig2f <- theta[3]  
  
  output <- yt^2*exp(-ft)+3*ft+(ft-omega-beta*ft1)^2/sig2f
  return(output)
}