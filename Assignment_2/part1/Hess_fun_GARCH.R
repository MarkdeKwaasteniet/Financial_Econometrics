# This file contains the average log-likelihood function without the link fuctions
# for the paramters, which can be used to obtain the Hessian matrix. 


Hess_fun_GARCH <- function(par,x){
  n <- length(x)
  
  omega <- par[1]                
  alpha <- par[2]
  beta <- par[3]
  
  sig2 <- rep(0,n)
  sig2[1] <- var(x) 
  
  for(t in 2:n){
    sig2[t] <- omega + alpha*x[t-1]^2 + beta*sig2[t-1]
  }
  
  l <- -(1/2)*log(2*pi) - (1/2)*log(sig2) - (1/2)*x^2/sig2
  llik <- mean(l) 
  
  return(llik) 
}