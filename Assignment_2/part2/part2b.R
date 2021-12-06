### Financial Engineering ####
# Assignment part 2b: dynamic regression and CAPM
# Daniel van Hanswijk(2726843)
# Mark de Kwaasteniet(2649271)
# Martijn van Welsenes(2713315)
# Sem Wierdsma(2670919)
##############################

# Clean workspace and set workdirectory ----
rm(list=ls())    

# Change to your directory
#setwd("~/Documents/Econometrics/Financial Econometrics/Assignment part 2")


# Load packages and modules ----
library(quantmod) 
library(quadprog)

source("llik_OD_regression.R")
source("llik_fun_GARCH.R")
source("sim_m_REG.R") 



### QUESTION 1 ####
# Plot the log-returns of each asset as well as the market log-returns. 
# Estimate the betas of the CAPM model for each of the three assets using OLS.

# load the stock market prices from data files 
p_market <- scan("./data/market.txt") # load SP500
p_msft <- scan("./data/MSFT.txt") # load Microsoft
p_bac <- scan("./data/BAC.txt") # load Bank of America
p_xom <- scan("./data/XOM.txt")  # load Exxon Mobil


# obtain log-returns of stock market prices
r_market <- diff(log(p_market)) * 100
r_msft <- diff(log(p_msft)) * 100
r_bac <- diff(log(p_bac)) * 100
r_xom <- diff(log(p_xom)) * 100

# combine all log-returns in one matrix
x <- cbind(r_msft, r_bac, r_xom, r_market)


# Plot the log-returns of the stock market prices
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(r_msft,type="l",main = "Log-return MSFT ",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(r_bac,type="l",main = "Log-return BAC",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(r_xom,type="l",main = "Log-return XOM",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(r_market,type="l",main = "Log-return market",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")


n <- length(x[1,])-1
beta_ols <- rep(0, n)

# calculate the beta estimates using OLS
# loop through all individual stock, compare with market proxy
for (i in 1:n){
  xt = r_market
  yt = x[, i]
  # beta_ols[i] <- 1/(r_market %*% r_market) * (r_market %*% x[, i])
  beta_ols[i] <- cov(xt,yt)/var(xt)
}

# print beta estimates
beta_ols

### QUESTION 3 ####
# Estimate by ML a CAPM model with an observation-driven 
# dynamic coefficient beta_t for each of the assets. 

# initialize parameter matrices
omega_hat <- rep(0, n)
phi_hat <- rep(0, n)
alpha_hat <- rep(0, n)
sigma2_hat <- rep(0, n)

# Estimate CAPM model parameters for the three assets 
for (i in 1:n){
  xt = r_market
  yt = x[, i]
  
  a <- 0.2/sd(xt*yt) # initial value for alpha
  phi <- 0.9  # initial value for beta
  omega <- (cov(xt,yt)/var(xt))*(1-phi) # initial value for omega
  sig2 <- var(yt)
  
  par_ini <- c(omega,log(phi/(1-phi)),log(a),log(sig2))
  
  # estimate parameters using ML
  est <- optim(par=par_ini,fn=function(par)-llik_OD_regression(yt,xt,par), method = "BFGS")
  
  omega_hat[i] <- est$par[1]
  phi_hat[i] <- exp(est$par[2])/(1+exp(est$par[2]))
  alpha_hat[i] <- exp(est$par[3])
  sigma2_hat[i] <- exp(est$par[4])
}

# print parameters
omega_hat
phi_hat
alpha_hat
sigma2_hat


# Now estimate the time varying betas
t_len <- length(r_market)
beta_t <- matrix(0, nrow=t_len, ncol=n)

# initialize beta at time 1
beta_t[1, 1] <- omega_hat[1]/(1-phi_hat[1])
beta_t[1, 2] <- omega_hat[2]/(1-phi_hat[2])
beta_t[1, 3] <- omega_hat[3]/(1-phi_hat[3])

# simulate the dynamic coefficient beta
for (i in 1:n){
  xt = r_market
  yt = x[, i]
  
  # for asset i, simulate the dynamic beta
  for (t in 2:t_len){
    beta_t[t, i] <- omega_hat[i]+phi_hat[i]*beta_t[t-1, i]+alpha_hat[i]*(yt[t-1]-beta_t[t-1, i]*xt[t-1])*xt[t-1];
  }
}

# plot the betas of the three assets
par(mfrow=c(3,1),mar=c(4.1,4.1,1.1,2.1))
plot(beta_t[,1],type="l",main = "dynamic coefficient beta MSFT",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(beta_t[,2],type="l",main = "dynamic coefficient beta BAC",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(beta_t[,3],type="l",main = "dynamic coefficient beta XOM",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")



### QUESTION 4 ####
# Obtain the exposition to the market of each asset at time T + 1


beta_pred <- rep(0, n)
# predict beta at T+1 using data at time T (this is t, since the loop ended at t=T at the prev question)
for (i in 1:n){
  xt = r_market[t]
  yt = x[t, i]

  beta_pred[i] <- omega_hat[i]+phi_hat[i]*beta_t[t, i]+alpha_hat[i]*(yt-beta_t[t, i]*xt)*xt;
}

beta_pred


### QUESTION 5 ####
# For each asset, obtain the time-varying beta based on the bivariate CCC model.
#  
# 

beta_t_ccc <- matrix(0, nrow=t_len, ncol=n)
for (i in 1:n){
  # Estimate parameters of asset i
  
  # initialize parameters
  alpha_ini <- 0.2
  beta_ini <- 0.6
  omega_ini <- var(x[,i]) * (1-alpha_ini - beta_ini)
  par_ini <- c(log(omega_ini),
               log(alpha_ini/(1-alpha_ini)),
               log(beta_ini/(1-beta_ini)))
  
  # maximize log likelihood to obtain the parameter estimates of the CCC model
  est1 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,i]), method = "BFGS")
  omega_hat1 <- exp(est1$par[1])  # took the log in the loglikelihood function
  alpha_hat1 <- exp(est1$par[2])/(1+exp(est1$par[2])) # make a logistic in the llik function
  beta_hat1 <- exp(est1$par[3])/(1+exp(est1$par[3]))
  
  # Estimate parameters of the market 
  # KEEP THIS CONSTANT: THIS IS THE MARKET 
  alpha_ini <- 0.2
  beta_ini <- 0.6
  omega_ini <- var(x[,4]) * (1-alpha_ini - beta_ini)
  par_ini <- c(log(omega_ini),
               log(alpha_ini/(1-alpha_ini)),
               log(beta_ini/(1-beta_ini)))
  
  est2 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,4]), method = "BFGS")
  omega_hat2 <- exp(est2$par[1])  # took the log in the loglikelihood function
  alpha_hat2 <- exp(est2$par[2])/(1+exp(est2$par[2])) # make a logistic in the llik function
  beta_hat2 <- exp(est2$par[3])/(1+exp(est2$par[3]))
  
  # Find the estimated conditional variances and covariances 
  # initialize 
  t_len <- length(x[,1])
  s1 <- rep(0, n)
  s2 <- rep(0, n)
  
  #as an initial condition for the variance, use the sample variance
  s1[1] <- var(x[,i])
  s2[1] <- var(x[,4])
  
  # obtain conditional variance of asset i usign CCC model
  for (t in 2:t_len){
    s1[t] <- omega_hat1 + alpha_hat1*x[t-1, i]^2 + beta_hat1*s1[t-1]
    s2[t] <- omega_hat2 + alpha_hat2*x[t-1, 4]^2 + beta_hat2*s2[t-1]
  }
  
  
  # obtain the standardized errors from the log-return series
  # obtain residuals
  e1 <- x[,i]/sqrt(s1)
  e2 <- x[,4]/sqrt(s2)
  
  # calculate the correlation between the residuals of the first and second series
  r <- cor(e1,e2)
  s12 <- r*sqrt(s1)*sqrt(s2)
  
  
  # calculate the time-varying beta
  for (t in 1:t_len){
    beta_t_ccc[t, i] <- s12[t] / s2[t]
  }
  
}

# print results
par(mfrow=c(3,1),mar=c(4.1,4.1,1.1,2.1))
plot(beta_t_ccc[,1],type="l",main = "dynamic coefficient beta MSFT (CCC)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(beta_t_ccc[,2],type="l",main = "dynamic coefficient beta BAC (CCC)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(beta_t_ccc[,3],type="l",main = "dynamic coefficient beta XOM (CCC)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

### QUESTION 6 ####
# Estimate by indirect inference a dynamic CAPM model for each of the three assets.

a0_hat <- rep(0, n)
a1_hat <- rep(0, n)
s_eta_hat <- rep(0, n)
s_eps_hat <- rep(0, n)

for (i in 1:n){
  xt = r_market
  yt = x[, i]
  
  #Obtain OLS estimate for beta and calculate residuals yr
  hb <- cov(yt,xt)/var(xt)
  yr <- yt-hb*xt
  
  xy <- yr*xt
  
  # calculate the auto-correlation function
  acvfxy <- acf(xy, lag.max=15, type ="covariance", plot=F)$acf[-1]
  
  sample_m <- c(var(yr),hb,acvfxy)
  
  M <- 20
  H <- M*length(xt)
  
  #Generate errors for the simulation
  set.seed(123)
  eta <- rnorm(H)
  eps <- rnorm(H)
  e <- cbind(eta,eps)
  
  #Choose initial parameter values for the optimization
  a1 <- 0.95 
  a0 <- cov(xt,yt)/var(xt)*(1-a1)
  s_eta <- 0.2
  s_eps <- var(yr)
  
  par_ini <- c(a0, log(a1/(1-a1)), log(s_eta), log(s_eps))
  
  #Minimize criterion function
  est <- optim(par=par_ini,fn=function(par) mean((sim_m_REG(e,xt,par)-sample_m)^2), method = "BFGS")

  #Obtain parameter estimate using the link functions
  a0_hat[i] <- est$par[1]
  a1_hat[i] <- exp(est$par[2])/(1+exp(est$par[2]))
  s_eta_hat[i] <- exp(est$par[3])
  s_eps_hat[i] <- exp(est$par[4])
}
#show results
a0_hat
a1_hat
s_eta_hat
s_eps_hat
