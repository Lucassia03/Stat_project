# Bayesian logistic regression and Poisson regression
# Random Walk Metropolis-Hastings

rm(list=ls())

library(mvtnorm)
library(AER)
library(COUNT)

# Definition of the function performing Markov-Chain Monte-Carlo
# for Bayesian Logistic and Poisson regression. It produces 
# a Markov Chain whose limit distribution is the posterior distribution
# of the regression coefficients given the data, from it we can compute
# the posterior means. 

# arguments:
# - y response
# - X design matrix
# - n_it number of iterations
# - b_0 starting value of the chain
# - prior list of prior mean and covariance matrix
# - S_prop covariance matrix of the proposal distribution (also a multivariate normal)
# - family either possion or bernoulli depending if we want to do
# a poisson bayesian regression or logistic bayesian regression

mcmc_glm <- function(y,X,n_it,b_0,prior,S_prop,family){
  
  samples <- matrix(NA, nrow = n_it, ncol = p)
  n_acc <- 0
  
  # Prior: beta ~ N(mu_prior, S_prior)
  mu_prior <- prior[[1]]
  S_prior <- prior[[2]]
  
  # Prior log-density 
  lprior <- function(b, mu, S){
    dmvnorm(as.numeric(b), mean = mu, sigma = S, log = TRUE)
  }
  
  if(family=="bernoulli"){
    loglik <- function(b,X,y){
      y <- as.matrix(y)
      X <- as.matrix(X)
      b <- matrix(b,ncol=1)
      eta <- as.vector(X%*%b)
      ll <- sum(y*eta-log(1+exp(eta)))
      return(ll)
    }
  }
  
  if(family=="poisson"){
    loglik <- function(b,X,y){
      y <- as.matrix(y)
      X <- as.matrix(X)
      b <- matrix(b,ncol=1)
      eta <- as.vector(X%*%b)
      ll <- sum(y*eta-log(factorial(y))-exp(eta))
      return(ll)
    }
  }
  
  # Posterior log-kernel
  lpost <- function(b, X, y, mu, S_prior){
    loglik(b, X, y) + lprior(b, mu, S_prior)
  }
  
  b_t <- b_0
  
  for(i in 1:n_it){
    
    # Proposal
    b_star <- as.numeric(rmvnorm(1, mean = b_t, sigma = S_prop))
    
    # Proposed posterior kernel
    lp_star <- lpost(b_star, X, y, mu=mu_prior, S_prior = S_prior)
    # previous iter posterior kernel
    lp_t <- lpost(b_t, X, y, mu=mu_prior, S_prior = S_prior)
    
    # log acceptance ratio
    log_ratio <- lp_star - lp_t
    log_acc <- min(0, log_ratio) # 0=log(1)
    
    if(log(runif(1)) <= log_acc){
      b_t <- b_star
      lp_t <- lp_star
      n_acc <- n_acc + 1
    }
    
    samples[i, ] <- b_t
  }
  
  acc_rate <- n_acc / n_it
  out <- list(samples,acc_rate)
  return(out)
}

# Bayesian Logistic regression

# Data loading / cleaning
data("CreditCard")
data <- CreditCard

# Response 0/1
y <- ifelse(data$card == "yes", 1, 0)

# Design matrix
X <- data[, !names(data) %in% c("share", "expenditure", "card")]

X$selfemp <- ifelse(data$selfemp == "yes", 1, 0)
X$owner   <- ifelse(data$owner   == "yes", 1, 0)

# Make numeric matrix
X <- as.matrix(X)

# Standardize predictors
X <- scale(X)

# Add intercept
n <- nrow(X)
X <- cbind(1, X)

# Dimensions
p <- ncol(X)

# Prior: beta ~ N(mu, S_prior)
mu <- rep(0, p)
prior_sd <- 1.0
S_prior <- (prior_sd^2) * diag(p)

prop_s <- 0.13  # tune this 
S_prop <- (prop_s^2) * diag(p)

n_it <- 100000
burn <- 20000

# Starting value (vector)
b_0 <- c(0, rep(0, p-1)) 

# MCMC
output <- mcmc_glm(y,X,n_it=100000,b_0=b_0,prior=list(mu,S_prior),S_prop=S_prop,family="bernoulli")

n_acc <- output[[2]]

# Basic diagnostics
samples <- output[[1]]
post <- samples[(burn+1):n_it,]
n_post <- nrow(post)

#install.packages("glue")
library(glue)
colnames(X)

par(mfrow=c(1,2))
colnam <- c("Intercept",colnames(X[,-1]))

for(j in 1:p){
  plot(samples[,j], type="l", main=glue("Trace: {colnam[j]}"), xlab="iter", ylab=glue("beta{j-1}"))
  abline(h = mean(samples[,j]),col="red",lwd=2)
  acf(post[,j], main=glue("ACF: {colnam[j]}"),lag.max=100)
}

# estimated posterior mean (with thinning)
post_est <- colMeans(post[seq(2,n_post,by=5),])

# Bayesian Poisson regression:

# data loading

library("COUNT")
data("azcabgptca")
data <- azcabgptca
y <- as.matrix(data$los)
X <- as.matrix(data[, which(names(data) != "los")])
X <- scale(X)

# Add intercept
n <- nrow(X)
X <- cbind(1, X)

# Dimensions
p <- ncol(X)

# Prior: beta ~ N(mu, S_prior)
mu <- rep(0, p)
prior_sd <- 1.0
S_prior <- (prior_sd^2) * diag(p)

prop_s <- 0.13  # tune this 
S_prop <- (prop_s^2) * diag(p)

n_it <- 100000
burn <- 20000

# Starting value (vector)
b_0 <- c(0, rep(0, p-1)) 

output <- mcmc_glm(y,X,n_it=100000,b_0=b_0,prior=list(mu,S_prior),S_prop=S_prop,family="poisson")

par(mfrow=c(1,2))
colnam <- c("Intercept",colnames(X[,-1]))

for(j in 1:p){
  plot(samples[,j], type="l", main=glue("Trace: {colnam[j]}"), xlab="iter", ylab=glue("beta{j-1}"))
  abline(h = mean(samples[,j]),col="red",lwd=2)
  acf(post[,j], main=glue("ACF: {colnam[j]}"),lag.max=100)
}
