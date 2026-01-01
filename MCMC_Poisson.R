#Build a log-posterior

logposterior_poisson <- function(beta, X, y, s = 2.5) {
  eta <- as.vector(X %*% beta)
  loglik <- sum(y * eta - exp(eta))
  logprior <- -0.5 * sum(beta^2 / s^2)
  loglik + logprior
}
#Parameter definition
true_mu <- 0
azcabgptca <- na.omit(azcabgptca)
n <- nrow(azcabgptca)
y <- azcabgptca$los
X_data <- azcabgptca[, !(names(azcabgptca) %in% c("los"))]
X <- model.matrix(~ . - 1, data = X_data)
X <- scale(X)
p <- ncol(X)



iters <- 20000
proposal_sd <- 0.15   # tuning parameter (you may adjust)

beta_current <- rep(0, p)     # starting value
samples <- matrix(NA, iters, p)
accept <- 0
iters <- 20000
proposal_sd <- 1
#defining jumper standard deviation
jumper_stds <- rep(0.1, p)
samples <- matrix(NA, iters, p)
accept <- rep(0, p)



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
    
    log_ratio <- logposterior_poisson(beta_proposal, X, y) -
      logposterior_poisson(beta_current, X, y)
    
    if (log(runif(1)) < log_ratio) {
      beta_current <- beta_proposal
      accept[j] <- accept[j] + 1
    }
  }
  
  samples[i, ] <- beta_current
}

accept_rate <- accept / iters
accept_rate
