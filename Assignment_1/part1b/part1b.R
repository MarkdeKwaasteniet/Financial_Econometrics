### Financial Engineering ####
# Assignment part 1b: Portfolio management with multivariate GARCH models
# Daniel van Hanswijk(2726843)
# Mark de Kwaasteniet(2649271)
# Martijn van Welsenes(2713315)
# Sem Wierdsma(2670919)
##############################

# Clean workspace and set workdirectory ----
rm(list=ls())    
setwd("C:/Users/marti/OneDrive/Bureaublad/Universiteit van Amsterdam/Econometrie/Fin_eco/Assignment 1/R_code_assignment_1/part1b")



# Load packages and modules ----
library(quantmod) 
library(quadprog)
source("llik_fun_GARCH.R")
source("max_SR_portfolio.R")
source("llik_CT_sDVECH.R")


### QUESTION 1 ####
# Estimate a bivariate CCC model using the equation by equation approach
# >> 1. estimate a univariate GARCH model for each series of log returns
#    2. obtain the standardized errors from each of these series
#    3. estimate the correlation matrix from the residuals


# load the US and Honk Kong stock market prices from data files 
p_sp500 <- scan("./data/SeP500.txt")
p_hsi <- scan("./data/HSI.txt")

# obtain log-returns of US and Honk Kong stock market
r_sp500 <- diff(log(p_sp500)) * 100
r_hsi <- diff(log(p_hsi)) * 100

# combine the two log return series in a single matrix 
x <- cbind(r_sp500, r_hsi)



# Estimate a univariate GARCH model for the SP500 log-returns
# initialize parameters
alpha_ini <- 0.2
beta_ini <- 0.6
omega_ini <- var(x[,1]) * (1-alpha_ini - beta_ini)
par_ini <- c(log(omega_ini),
             log(alpha_ini/(1-alpha_ini)),
             log(beta_ini/(1-beta_ini)))

# maximize log likelihood to obtain the parameter estimates
est1 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,1]), method = "BFGS")
omega_hat1 <- exp(est1$par[1])  # took the log in the loglikelihood function
alpha_hat1 <- exp(est1$par[2])/(1+exp(est1$par[2])) # make a logistic in the llik function
beta_hat1 <- exp(est1$par[3])/(1+exp(est1$par[3]))

# Repeat for the HSI log-returns
alpha_ini <- 0.2
beta_ini <- 0.6
omega_ini <- var(x[,2]) * (1-alpha_ini - beta_ini)
par_ini <- c(log(omega_ini),
             log(alpha_ini/(1-alpha_ini)),
             log(beta_ini/(1-beta_ini)))

est2 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,2]), method = "BFGS")
omega_hat2 <- exp(est2$par[1])  # took the log in the loglikelihood function
alpha_hat2 <- exp(est2$par[2])/(1+exp(est2$par[2])) # make a logistic in the llik function
beta_hat2 <- exp(est2$par[3])/(1+exp(est2$par[3]))

# now we have the parameter estimates, which are:
# sp500: omega_hat1, alpha_hat1 & omega_hat1
# hsi: omega_hat2, alpha_hat2 & omega_hat2

# Find the estimated conditional variances and covariances 
# initialize 
n <- length(x[,1])
s1 <- rep(0, n)
s2 <- rep(0, n)

#as an initial condition for the variance, use the sample variance
s1[1] <- var(x[,1])
s2[1] <- var(x[,2])


for (t in 2:n){
  s1[t] <- omega_hat1 + alpha_hat1*x[t-1, 1]^2 + beta_hat1*s1[t-1]
  s2[t] <- omega_hat2 + alpha_hat2*x[t-1, 2]^2 + beta_hat2*s2[t-1]
}


# obtain the standardized errors from the log-return series
# obtain residuals
e1 <- x[,1]/sqrt(s1)
e2 <- x[,2]/sqrt(s2)

# calculate the correlation between the residuals of the first and second series
r <- cor(e1,e2)
s12 <- r*sqrt(s1)*sqrt(s2)


# Plot the estimated conditional variances, covariances and correlation
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(s1,type="l",main = "Conditional variance S&P 500",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s2,type="l",main = "Conditional variance HSI",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s12,type="l",main = "Conditional covariance",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(rep(r,n),type="l",main = "Conditional correlation",ylab="",xlab="")
abline(h=r)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")





### QUESTION 2 ####
#  Obtain the conditional variance and the a-Var for the portfolio of the bank
#
#

# Obtain the conditional variance of the portfolio consisting op 70% sp500 and 30% hsi
k1 = 0.7  # weight sp500
k2 = 0.3  # weight hsi

# initialize variance, return and a-VaRportfolio p
sp <- rep(0, n)
xp <- rep(0, n)
VaR <- rep(0, n)

for (t in 1:n){
  xp[t] <- k1*r_sp500[t] + k2*r_hsi[t]
  
  # Calculate variance
  sp[t] <- k1^2 * s1[t] + k2^2 * s2[t] + 2*k1*k2*s12[t]  # slide 25 
  
  # calculate the a-VaR
  VaR[t] <- qnorm(0.01, 0, sqrt(sp[t]))
}

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(sp,type="l",main = "Conditional Variance",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(VaR,type="l",main = "a-VaR",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")




### QUESTION 3 ####
#  Obtain the optimal portfolio weights in terms of the Sharpe ratio
# 
#
#


# initialize weights matrix
kt <- matrix(0,nrow=n,ncol=2)

# set the conditional mean as the sample mean of the log returns
mu1 <- mean(x[,1])
mu2 <- mean(x[,2])
mut <- cbind(mu1,mu2) # combine means together in a matrix

# Obtain the optimal portfolio weights at each point in time

for (t in 1:n){
  # Obtain the conditional variance matrix
  SIGMAt <- cbind(c(s1[t],s12[t]),c(s12[t],s2[t]))
  
  # obtain the optimal weights
  kt[t,] <- max_SR_portfolio(mut,SIGMAt)
}

# Predict the optimal weight at time T+1
# predict the conditional variance, covariance and correlation at T+1 
# Using data at time T

s1_pred <- omega_hat1 + alpha_hat1*x[t, 1]^2 + beta_hat1*s1[t]
s2_pred <- omega_hat2 + alpha_hat2*x[t, 2]^2 + beta_hat2*s2[t]
s12_pred <- r*sqrt(s1_pred)*sqrt(s2_pred)

# Obtain the conditional variance matrix
SIGMAt_pred <- cbind(c(s1_pred,s12_pred),c(s12_pred,s2_pred))

# Obtain the optimal portfolio weighs
kt_pred <- max_SR_portfolio(mut,SIGMAt_pred)

# Plot the optimal weights
par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(kt[,1], type="l",main = "Portfolio weight S&P 500",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(kt[,2],type="l",main = "Portfolio weight HSI",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")



### QUESTION 4 ####
#  Repeat Question 1 & 2, now estimating a bivariate sDVECH model
#
#

# initialize parameters
alpha_ini <- 0.2
beta_ini <- 0.2
par_ini <- c(log(alpha_ini/(1-alpha_ini)),log(beta_ini/(1-beta_ini)))

est <- optim(par=par_ini, fn=function(par) -llik_CT_sDVECH(par, x), method="BFGS")

# map the parameters back to the true value
a_hat <- exp(est$par[1])/(1+exp(est$par[1]))
b_hat <- exp(est$par[2])/(1+exp(est$par[2]))

# make VECH matrix

VECHt <- matrix(0, nrow=n, ncol=3)

C <- cov(x)
VECHt[1,] <- c(C[1,1], C[1,2], C[2,2])

# Simulate for each timestep

for(t in 2:n){
  VECHt[t,1] <- C[1,1]*(1-a_hat-b_hat)+b_hat*VECHt[t-1,1]+a_hat*x[t-1,1]^2
  VECHt[t,3] <- C[2,2]*(1-a_hat-b_hat)+b_hat*VECHt[t-1,3]+a_hat*x[t-1,2]^2
  VECHt[t,2] <- C[1,2]*(1-a_hat-b_hat)+b_hat*VECHt[t-1,2]+a_hat*x[t-1,1]*x[t-1,2]
}

# Calculate the conditional correlation
sd1 <- sqrt(VECHt[,1])
sd2 <- sqrt(VECHt[,3])
corrt <- VECHt[,2]/(sd1*sd2)

# portfolio weights SP 500 and HSI are 0.7 and 0.3 respectively
k1 = 0.7
k2 = 0.3

# initialize variance, return and a-VaRportfolio p
sp <- rep(0, n)
xp <- rep(0, n)
VaR <- rep(0, n)



for (t in 1:n){
  # Calculate portfolio log return
  xp[t] <- k1*r_sp500[t] + k2*r_hsi[t]
  
  # Calculate variance
  sp[t] <- k1^2 * VECHt[t,1] + k2^2 * VECHt[t,3] + 2*k1*k2*VECHt[t,2]  # slide 25 
  
  # calculate the a-VaR
  # z-score at 1% precentile using qnorm
  VaR[t] <- qnorm(0.01, 0, sqrt(sp[t]))
}


# Plot the estimated conditional variances, covariances and correlation
par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(VECHt[,1],type="l",main = "Conditional variance S&P 500",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(VECHt[,3],type="l",main = "Conditional variance HSI",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(VECHt[,2],type="l",main = "Conditional covariance",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(corrt,type="l",main = "Conditional correlation",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

# Plot the conditional Var and a-VaR
par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(sp,type="l",main = "Conditional Variance",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(VaR,type="l",main = "a-VaR",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")



### QUESTION 5 ####
#  Obtain forecasts of the volatility of the banks portfolio for the next 52 weeks
#
#

# make VECH prediction matrix
weeks = 52
VECHt_pred <- matrix(0, nrow=weeks, ncol=3)
C <- cov(x)

# Predict the next time step
VECHt_pred[1,1] <- C[1,1]*(1-a_hat-b_hat)+b_hat*VECHt[t,1]+a_hat*x[t,1]^2
VECHt_pred[1,3] <- C[2,2]*(1-a_hat-b_hat)+b_hat*VECHt[t,3]+a_hat*x[t,2]^2
VECHt_pred[1,2] <- C[1,2]*(1-a_hat-b_hat)+b_hat*VECHt[t,2]+a_hat*x[t,1]*x[t,2]

# Predict for the  for each timestep, assuming that the expected returns squared
# is the variance 


for(t in 2:weeks){
  VECHt_pred[t,1] <- C[1,1]*(1-a_hat-b_hat) + ( b_hat+ a_hat)*VECHt_pred[t-1,1] 
  VECHt_pred[t,3] <- C[2,2]*(1-a_hat-b_hat) + (a_hat+ b_hat)*VECHt_pred[t-1,3]
  VECHt_pred[t,2] <- C[1,2]*(1-a_hat-b_hat) + ( a_hat+ b_hat)*VECHt_pred[t-1,2]
}

# portfolio weights SP 500 and HSI are 0.7 and 0.3 respectively
k1 = 0.7
k2 = 0.3

# initialize variance, return and a-VaRportfolio p
sp <- rep(0, weeks)

for (t in 1:weeks){

  sp[t] <- k1^2 * VECHt_pred[t,1] + k2^2 * VECHt_pred[t,3] + 2*k1*k2*VECHt_pred[t,2]  # slide 25 
  
}


par(mfrow=c(1,1),mar=c(4.1,4.1,1.1,2.1))
plot(sp,type="l",main = "Conditional variance portfolio",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")



### QUESTION 6 ####
#  Estimate a 3 dimensional CCC model for the three market returns
#
#

# load the DAX data
p_dax <- scan("./data/DAX.txt")

# obtain log-returns of US and Honk Kong stock market
r_dax <- diff(log(p_dax)) * 100

# combine the now three log return series in a single matrix 
x <- cbind(r_sp500, r_hsi, r_dax)


# Repeat Question 1, now also add estimated parameters for the DAX
# Estimate a univariate GARCH model for the SP500 log-returns
# initialize parameters
alpha_ini <- 0.2
beta_ini <- 0.6
omega_ini <- var(x[,1]) * (1-alpha_ini - beta_ini)
par_ini <- c(log(omega_ini),
             log(alpha_ini/(1-alpha_ini)),
             log(beta_ini/(1-beta_ini)))  # use exp to make sure that positive bla bla bla

# maximize log likelihood to obtain the parameter estimates
est1 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,1]), method = "BFGS")
omega_hat1 <- exp(est1$par[1])  # took the log in the loglikelihood function
alpha_hat1 <- exp(est1$par[2])/(1+exp(est1$par[2])) # make a logistic in the llik function
beta_hat1 <- exp(est1$par[3])/(1+exp(est1$par[3]))

# Repeat for the HSI log-returns
alpha_ini <- 0.2
beta_ini <- 0.6
omega_ini <- var(x[,2]) * (1-alpha_ini - beta_ini)
par_ini <- c(log(omega_ini),
             log(alpha_ini/(1-alpha_ini)),
             log(beta_ini/(1-beta_ini)))

est2 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,2]), method = "BFGS")
omega_hat2 <- exp(est2$par[1])  # took the log in the loglikelihood function
alpha_hat2 <- exp(est2$par[2])/(1+exp(est2$par[2])) # make a logistic in the llik function
beta_hat2 <- exp(est2$par[3])/(1+exp(est2$par[3]))

# Repeat for the DAX log-returns
alpha_ini <- 0.2
beta_ini <- 0.6
omega_ini <- var(x[,3]) * (1-alpha_ini - beta_ini)
par_ini <- c(log(omega_ini), 
             log(alpha_ini/(1-alpha_ini)),
             log(beta_ini/(1-beta_ini)))

est3 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,x[,3]), method = "BFGS")
(omega_hat3 <- exp(est3$par[1]))  # took the log in the loglikelihood function
(alpha_hat3 <- exp(est3$par[2])/(1+exp(est3$par[2]))) # make a logistic in the llik function
(beta_hat3 <- exp(est3$par[3])/(1+exp(est3$par[3])))

# Find the estimated conditional variances and covariances 
# initialize 
n <- length(x[,1])
s1 <- rep(0, n)
s2 <- rep(0, n)
s3 <- rep(0, n)

#as an initial condition for the variance, use the sample variance
s1[1] <- var(x[,1])
s2[1] <- var(x[,2])
s3[1] <- var(x[,3])


for (t in 2:n){
  s1[t] <- omega_hat1 + alpha_hat1*x[t-1, 1]^2 + beta_hat1*s1[t-1]
  s2[t] <- omega_hat2 + alpha_hat2*x[t-1, 2]^2 + beta_hat2*s2[t-1]
  s3[t] <- omega_hat3 + alpha_hat3*x[t-1, 3]^2 + beta_hat3*s3[t-1]
}

# obtain the standardized errors from the log-return series
# obtain residuals
e1 <- x[,1]/sqrt(s1)
e2 <- x[,2]/sqrt(s2)
e3 <- x[,3]/sqrt(s3)

# calculate the correlation between the residuals of the first and second series
r12 <- cor(e1, e2)
r13 <- cor(e1, e3)
r23 <- cor(e2, e3)

s12 <- r12*sqrt(s1)*sqrt(s2)
s13 <- r13*sqrt(s1)*sqrt(s3)
s23 <- r23*sqrt(s2)*sqrt(s3)


mu1 <- mean(x[,1])
mu2 <- mean(x[,2])
mu3 <- mean(x[,3])
mut <- cbind(mu1,mu2,mu3) # combine means together in a matrix

s1_pred <- omega_hat1 + alpha_hat1*x[t, 1]^2 + beta_hat1*s1[t]
s2_pred <- omega_hat2 + alpha_hat2*x[t, 2]^2 + beta_hat2*s2[t]
s3_pred <- omega_hat3 + alpha_hat3*x[t, 3]^2 + beta_hat3*s3[t]


s12_pred <- r12*sqrt(s1_pred)*sqrt(s2_pred)
s13_pred <- r13*sqrt(s1_pred)*sqrt(s3_pred)
s23_pred <- r23*sqrt(s2_pred)*sqrt(s3_pred)

# Obtain the conditional variance matrix
SIGMAt_pred <- cbind(c(s1_pred,s12_pred,s13_pred),
                     c(s12_pred,s2_pred,s23_pred),
                     c(s13_pred, s23_pred, s3_pred))


# Obtain the conditional variance matrix
# Obtain the optimal portfolio weights
kt_pred <- max_SR_portfolio(mut,SIGMAt_pred)


par(mfrow=c(3,3),mar=c(4.1,4.1,1.1,2.1))
plot(s1,type="l",main = "Conditional variance S&P 500",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(NULL,type="l",main = "",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(NULL,type="l",main = "",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s12,type="l",main = "Conditional covariance Sp500 and HSI",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s2,type="l",main = "Conditional variance HSI",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(NULL,type="l",main = "",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s13,type="l",main = "Conditional covariance SP500 and DAX",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s23,type="l",main = "Conditional covariance HSI & DAX",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(s3,type="l",main = "Conditional variance DAX",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

