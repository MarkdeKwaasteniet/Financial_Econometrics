### This file contains the code to obtain simulated moments from an SV model
# The input are parameter values and the error vectors et and vt
# For more details see the Lecture Notes.

sim_m_SV <- function(e,par){
  
  omega <- par[1]      
  beta <- exp(par[2])/(1+exp(par[2]))
  sig2f <- exp(par[3]) 
  
  H <- length(e[,1])
  
  epsilon <- e[,1] 
  eta <- sqrt(sig2f)*e[,2]
  
  x <- rep(0,H) 
  f <- rep(0,H) 
  
  f[1] <- omega/(1-beta)
  x[1] = exp(f[1]/2) * epsilon[1]
  
  for(t in 2:H){
    f[t] <- omega + beta * f[t-1] + eta[t]  # state equation
    x[t] <- exp(f[t]/2) * epsilon[t]       # observation equation
  }
  
  xa <- abs(x) 
  
  output <- c(var(x),kurtosis(x),cor(xa[2:H],xa[1:(H-1)]))
  return(output)
}
