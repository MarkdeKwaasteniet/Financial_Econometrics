### Financial Engineering ####
# Assignment part 1b: Portfolio management with multivariate GARCH models
# Daniel van Hanswijk(2726843)
# Mark de Kwaasteniet(2649271)
# Martijn van Welsenes(2713315)
# Sem Wierdsma(2670919)
##############################

# Clean workspace and set workdirectory ----
rm(list=ls()) 
#Change this to the appropriate working directory
setwd("C:/Users/marti/OneDrive/Bureaublad/Universiteit van Amsterdam/Econometrie/Fin_eco/Assignment 1/R_code_assignment_1/part1a")

########################
############ Problem 1
########################

SeP500 <- scan("SeP500.txt")
log_return <- 100 *diff(log(SeP500))
x <- log_return
x2 = x^2 #Obtain the squared log-returns
n <- length(x)

#Plot the index, log-return and squared log-return
par(mfrow=c(3,1),mar=c(4.1,4.1,1.1,2.1)) 
plot(SeP500, type = "l", main="S&P500", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(x, type = "l", main="log return S&P500", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(x2, type = "l", main=" squared log return S&P500", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

#Plot the ACF of log returns and squared log returns
par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1)) 
acf(x,main="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
title("ACF  log-returns SeP500", line = 0.3)
acf(x2, main="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
title("ACF squared log-returns SeP500", line = 0.3)


########################
############ Problem 2
########################


#At time T the price is increasing. Looking at the AFC we see that there is a
#small negative correlation for 1 time lag, which is just about more significant
#than 5%. This means that at the next time T+1 the index will likely decrease
#and thus expect negative returns. ALtough log returns are weakly correlated


########################
############ Problem 3
########################

source("llik_functions/llik_fun_GARCH_1_3.R")
a <- 0.01 # initial value for alpha
b1 <- 0.9  # initial value for beta1
b2 <- 0.01 # initial value for beta2
b3 <- 0.01 # initial value for beta3
omega <- var(x)*(1-a-b1-b2-b3) # initial value for omega
#Use logarithmic-link function for the parameters
par_ini <- c(log(omega),log(a/(1-a)),log(b1/(1-b1)),log(b2/(1-b2)),log(b3/(1-b3)))

est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH_1_3(par,x), method = "BFGS")

#Use logarithmic link function backwards again
omega_hat <- exp(est$par[1])
alpha_hat <- exp(est$par[2])/(1+exp(est$par[2]))
beta1_hat <- exp(est$par[3])/(1+exp(est$par[3]))
beta2_hat <- exp(est$par[4])/(1+exp(est$par[4]))
beta3_hat <- exp(est$par[5])/(1+exp(est$par[5]))

theta_1_3 <- c(omega_hat,alpha_hat,beta1_hat, beta2_hat, beta3_hat)
cat("The parameter estimates are:", round(theta_1_3,4) )

llik_1_3 <- -est$value*length(x) #Calculate log likelihood
aic_1_3 <- 2*length(par_ini)-2*llik_1_3 #Obtain AIC and BIC already for problem 5
bic_1_3 <- log(n)*length(par_ini)-2*llik_1_3
cat("The log-likelihood value is:", llik_1_3)

cat("Exit flag:", est$convergence) # zero indicates succesfull optimization

#Calculate the estimated conditional variance with the updating equation
n <- length(x)
sig2 <- rep(0,n)
sig2[3] <- var(x) #initialize volatility at unconditional variance
sig2[2] <- var(x) # we need 3 timesteps back because GARCH(1,3)
sig2[1] <- var(x)
for(t in 4:n){
  sig2[t] <- (theta_1_3[1] + theta_1_3[2]*x[t-1]^2 + theta_1_3[3]*sig2[t-1]
              + theta_1_3[4]*sig2[t-2] + theta_1_3[5]*sig2[t-3])
}

#Plot estimated conditional variance
plot(sig2,type="l", main="Estimated conditional variance",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################
############ Problem 4
########################

u = x/sqrt(sig2) # obtain residuals
par(mfrow=c(1,2),mar=c(4.1,4.1,1.1,2.1))
plot(acf(u^2,plot=F)[1:20]) # Plot the acf of residuals squared
abline(h=0)
axis(1, at=c(1,25)) #start label at 1

qqnorm(u)#Plot the qq plot of the residuals
qqline(u) #add 45degree line

library(tseries)
jarque.bera.test(u)


########################
############ Problem 5
########################

#Set initial values for all models to come
a1 <- 0.01 # initial value for alpha1
a2 <- 0.01 # initial value for alpha2
b1 <- 0.79  # initial value for beta1
b2 <- 0.01 # initial value for beta2
b3 <- 0.01 # initial value for beta3

################################################################################
# GARCH(1,1)
source("llik_functions/llik_fun_GARCH_1_1.R")
omega <- var(x)*(1-a1-b1) # initial value for omega
par_ini <- c(log(omega),log(a1/(1-a1)),log(b1/(1-b1)))

est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH_1_1(par,x), method = "BFGS")
llik_1_1 = -est$value*length(x)

omega_hat <- exp(est$par[1])
alpha_hat <- exp(est$par[2])/(1+exp(est$par[2]))
beta_hat <- exp(est$par[3])/(1+exp(est$par[3]))

theta_1_1 <- c(omega_hat,alpha_hat,beta_hat)
cat("The parameter estimates are:", theta_1_1)

aic_1_1 <- 2*length(par_ini)-2*llik_1_1
bic_1_1 <- log(n)*length(par_ini)-2*llik_1_1


################################################################################
#GARCH(1,2)
source("llik_functions/llik_fun_GARCH_1_2.R")
omega <- var(x)*(1-a1-b1-b2) # initial value for omega
par_ini <- c(log(omega),log(a1/(1-a1)),log(b1/(1-b1)), log(b2/(1-b2)))

est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH_1_2(par,x), method = "BFGS")
llik_1_2 = -est$value*length(x)

omega_hat <- exp(est$par[1])
alpha_hat <- exp(est$par[2])/(1+exp(est$par[2]))
beta1_hat <- exp(est$par[3])/(1+exp(est$par[3]))
beta2_hat <- exp(est$par[4])/(1+exp(est$par[4]))

theta_1_2 <- c(omega_hat,alpha_hat,beta1_hat, beta2_hat)
cat("The parameter estimates are:", theta_1_2)

aic_1_2 <- 2*length(par_ini)-2*llik_1_2
bic_1_2 <- log(n)*length(par_ini)-2*llik_1_2

################################################################################
##GARCH(2,1)
source("llik_functions/llik_fun_GARCH_2_1.R")
omega <- var(x)*(1-a1-a2-b1) # initial value for omega
par_ini <- c(log(omega),log(a1/(1-a1)),log(a2/(1-a2)), log(b1/(1-b1)))

est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH_2_1(par,x), method = "BFGS")
llik_2_1 = -est$value*length(x)

omega_hat <- exp(est$par[1])
alpha1_hat <- exp(est$par[2])/(1+exp(est$par[2]))
alpha2_hat <- exp(est$par[3])/(1+exp(est$par[3]))
beta_hat <- exp(est$par[4])/(1+exp(est$par[4]))

theta_2_1 <- c(omega_hat,alpha1_hat,alpha2_hat, beta_hat)
cat("The parameter estimates are:", theta_2_1)


aic_2_1 <- 2*length(par_ini)-2*llik_2_1
bic_2_1 <- log(n)*length(par_ini)-2*llik_2_1

################################################################################
##GARCH(2,2)
source("llik_functions/llik_fun_GARCH_2_2.R")
omega <- var(x)*(1-a1-a2-b1-b2) # initial value for omega
par_ini <- c(log(omega),log(a1/(1-a1)),log(a2/(1-a2)), log(b1/(1-b1)), log(b2/(1-b2)))

est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH_2_2(par,x), method = "BFGS")
llik_2_2 = -est$value*length(x)

omega_hat <- exp(est$par[1])
alpha1_hat <- exp(est$par[2])/(1+exp(est$par[2]))
alpha2_hat <- exp(est$par[3])/(1+exp(est$par[3]))
beta1_hat <- exp(est$par[4])/(1+exp(est$par[4]))
beta2_hat <- exp(est$par[5])/(1+exp(est$par[5]))

theta_2_2 <- c(omega_hat,alpha1_hat,alpha2_hat, beta1_hat, beta2_hat)
cat("The parameter estimates are:", theta_2_2)

aic_2_2 <- 2*length(par_ini)-2*llik_2_2
bic_2_2 <- log(n)*length(par_ini)-2*llik_2_2

cat('The AIC values are ')
cat('GARCH(1,1)', aic_1_1)
cat('GARCH(1,2)', aic_1_2)
cat('GARCH(1,3)', aic_1_3)
cat('GARCH(2,1)', aic_2_1)
cat('GARCH(2,2)', aic_2_2)

cat('The BIC values are ')
cat('GARCH(1,1)', bic_1_1)
cat('GARCH(1,2)', bic_1_2)
cat('GARCH(1,3)', bic_1_3)
cat('GARCH(2,1)', bic_2_1)
cat('GARCH(2,2)', bic_2_2)

#The model with lowest AIC and also BIC is the GARCH(1,1) model. 
#Increasing p raises the infromation criterion and also increasing q 
# raises the infromation criterion. So the model with lowest p,q in this case
# is the 'best' model

par(mfrow=c(3,2),mar=c(4.1,4.1,1.1,2.1))

# NEXT we calculate the estimated conditional variance of all plots
################################################################################
##GARCH(1,1)
sig2 <- rep(0,n)
sig2[1] <- var(x) #initialize volatility at unconditional variance

for(t in 2:n){
  sig2[t] <- (theta_1_1[1] + theta_1_1[2]*x[t-1]^2 + theta_1_1[3]*sig2[t-1])
}
plot(sig2,type="l", main="GARCH(1,1)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

################################################################################
##GARCH(1,2)
sig2 <- rep(0,n)
sig2[1] <- var(x) #initialize volatility at unconditional variance
sig2[2] <- var(x)

for(t in 3:n){
  sig2[t] <- (theta_1_2[1] + theta_1_2[2]*x[t-1]^2 + 
                theta_1_2[3]*sig2[t-1] + theta_1_2[4]*sig2[t-2])
}
plot(sig2,type="l", main="GARCH(1,2)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

################################################################################
##GARCH(1,3)
sig2 <- rep(0,n)
sig2[1] <- var(x) #initialize volatility at unconditional variance
sig2[2] <- var(x)
sig2[3] <- var(x)

for(t in 4:n){
  sig2[t] <- (theta_1_3[1] + theta_1_3[2]*x[t-1]^2 + theta_1_3[3]*sig2[t-1]
              + theta_1_3[4]*sig2[t-2] + theta_1_3[4]*sig2[t-3])
}

plot(sig2,type="l", main="GARCH(1,3)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

################################################################################
##GARCH(2,1)
sig2 <- rep(0,n)
sig2[1] <- var(x) #initialize volatility at unconditional variance
sig2[2] <- var(x)

for(t in 3:n){
  sig2[t] <- (theta_2_1[1] + theta_2_1[2]*x[t-1]^2 + 
                theta_2_1[3]*x[t-2]^2 + theta_2_1[4]*sig2[t-1])
}

plot(sig2,type="l", main="GARCH(2,1)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

################################################################################
##GARCH(2,2)
sig2 <- rep(0,n)
sig2[1] <- var(x) #initialize volatility at unconditional variance
sig2[2] <- var(x)

for(t in 3:n){
  sig2[t] <- (theta_2_2[1] + theta_2_2[2]*x[t-1]^2 + theta_2_2[3]*x[t-2]^2 + 
                theta_2_2[4]*sig2[t-1] + theta_2_2[5]*sig2[t-2])
}

plot(sig2,type="l", main="GARCH(2,2)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################
############ Problem 6
########################

sig2 <- rep(0,n)
sig2[1] <- var(x) #initialize volatility at unconditional variance

for(t in 2:n){
  sig2[t] <- (theta_1_1[1] + theta_1_1[2]*x[t-1]^2 + theta_1_1[3]*sig2[t-1])
}

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(x,type="l", main="log-return",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(sig2,type="l",col=2, main="Estimated conditional variance",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

sig2_next = theta_1_1[1] + theta_1_1[2]*x[n]^2 + theta_1_1[3]*sig2[n]
drop_prob <- pnorm(-5,0,sqrt(sig2_next))
cat('The probability of dropping 5% is', drop_prob*100, '%')
#probably the line will not drop with more than 5% 
#next week because, at the moment of t+1 the value is 
#0.92%

#7

library("tseries")
u<-x/sqrt(sig2) #obtain residuals
K <- 1/n * sum(u**4)/var(u)^2
hist(u)
