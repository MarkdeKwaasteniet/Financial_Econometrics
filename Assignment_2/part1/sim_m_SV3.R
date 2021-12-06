### This file contains the code to obtain simulated moments from an SV model
# The input are parameter values and the error vectors et and vt
# For more details see the Lecture Notes.

sim_m_SV3 <- function(e,par){
  
  omega <- par[1]      
  beta <- exp(par[2])/(1+exp(par[2]))
  beta1 <- exp(par[3])/(1+exp(par[3]))
  sig2f <- exp(par[4]) 
  
  H <- length(e[,1])
  
  epsilon <- e[,1] 
  eta <- sqrt(sig2f)*e[,2]
  
  x <- rep(0,H) 
  f <- rep(0,H) 
  
  f[1] <- omega/(1-beta-beta1)
  f[2] <- omega/(1-beta-beta1)
  x[1] = exp(f[1]/2) * epsilon[1]
  
  for(t in 3:H){
    f[t] <- omega + beta * f[t-1] + beta1 * f[t-2] + eta[t]  # state equation
    x[t] <- exp(f[t]/2) * epsilon[t]       # observation equation
  }
  
  xa <- abs(x)
  mean_x <- mean(xa)
  # obtain the autocovariance function
  acv15 <-acf(xa, lag.max = 15, type = "covariance", plot = F)$acf
  
  output <- c(mean_x, acv15)
  return(output)
}
