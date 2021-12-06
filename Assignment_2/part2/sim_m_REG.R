sim_m_REG <- function(e,x,par){
  
  a0 <- par[1]
  a1 <- exp(par[2])/(1+exp(par[2]))
  s_eta <- exp(par[3]) 
  s_eps <- exp(par[4])
  
  output <- 0
  
  n <- length(x)
  H <- length(e[,1])
  M <- H/n
  
  eta <- sqrt(s_eta)*e[,1]
  eps <- sqrt(s_eps)*e[,2]
  
  for(m in 1:M){
    
    b <- rep(0,n)
    f <- rep(0, n)
    f[1] <- a0/(1-a1)
    
    for(t in 2:n){
      f[t] <- a0+a1*f[t-1]+eta[(m-1)*n+t]
      b[t] <- exp(f[t])
      
    }
    
    y <- b*x+eps[((m-1)*n+1):(m*n)]
    
    hb <- cov(y,x)/var(x)
    yr <- y-hb*x
    
    xy <- yr*x
    acvfxy <- acf(xy, lag.max=15, type ="covariance", plot=F)$acf[-1]
    
    output <- output+c(var(yr),hb,acvfxy)/M
  }
   return(output)
  
}

