#Build a log-posterior

logposterior_logistic <- function(beta, X, y, s0=5, s=2.5) {
  eta <- as.vector(X %*% beta)
  loglik <- sum(y * eta - log1p(exp(eta)))  # logistic log-likelihood
  logprior <- -0.5 * (beta[1]^2 / s0^2 + sum(beta[-1]^2) / s^2)  # Normal priors
  loglik + logprior
}
#Parameter definition
true_mu <- 1
CreditCard <- na.omit(CreditCard)
n <- nrow(CreditCard)
y <- ifelse(CreditCard$card == "yes", 1, 0)
X_data <- CreditCard[, !(names(CreditCard) %in% c("card", "share", "expenditure", "income", "age", "months"))]
X <- model.matrix(~ ., data = X_data)
X[, -1] <- scale(X[, -1])
p <- ncol(X)


iters <- 20000
proposal_sd <- 0.15   # tuning parameter (you may adjust)

beta_current <- rep(0, p)     # starting value
samples <- matrix(NA, iters, p)
accept <- 0
iters <- 20000
proposal_sd <- 1
#defining jumper standar deviation
jumper_stds <- rep(0.3, p)
samples <- matrix(NA, iters, p)
accept <- rep(0, p)
pilot_iters <- 2000

for (i in 1:iters) {
  update_order <- sample(1:p) 
  for (j in update_order) {
    
    
    beta_proposal <- beta_current
    beta_proposal[j] <- beta_current[j] + rnorm(1, 0, jumper_stds[j])
    
    log_ratio <- logposterior_logistic(beta_proposal, X, y) -
      logposterior_logistic(beta_current, X, y)
    
    if (log(runif(1)) < log_ratio) {
      beta_current <- beta_proposal
      accept[j] <- accept[j] + 1
    }
  }
  
  samples[i, ] <- beta_current
}

accept_rate <- accept / iters
accept_rate



accept_rate <- accept / pilot_iters

#computing Optimal Jump distribution
for (j in 1:p) {
  if (accept_rate[j] > 0.7)
    jumper_stds[j] <- jumper_stds[j] * 1.3
  
  if (accept_rate[j] < 0.3)
    jumper_stds[j] <- jumper_stds[j] * 0.7
}




#Final MCMC

iters <- 20000
samples <- matrix(NA, iters, p)
accept <- rep(0, p)

for (i in 1:iters) {
  
  update_order <- sample(1:p)
  
  for (j in update_order) {
    
    beta_proposal <- beta_current
    beta_proposal[j] <- beta_current[j] +
      rnorm(1, 0, jumper_stds[j])
    
    log_ratio <- logposterior_logistic(beta_proposal, X, y) -
      logposterior_logistic(beta_current, X, y)
    
    if (log(runif(1)) < log_ratio) {
      beta_current <- beta_proposal
      accept[j] <- accept[j] + 1
    }
  }
  
  samples[i, ] <- beta_current
}

accept_rate_final <- accept / iters
accept_rate_final


