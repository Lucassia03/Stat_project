#Log-posterior


log1pexp <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))

logposterior_logistic <- function(beta, X, y, s0 = 5, s = 2.5) {
  eta <- drop(X %*% beta)
  loglik <- sum(y * eta - log1pexp(eta))
  logprior <- -0.5 * (beta[1]^2 / s0^2 + sum(beta[-1]^2) / s^2)
  loglik + logprior
}
# Data preparation


CreditCard2 <- na.omit(CreditCard)

to01 <- function(x) as.integer(tolower(as.character(x)) %in% c("yes","1","true","t"))

CreditCard2$card    <- to01(CreditCard2$card)
CreditCard2$owner   <- to01(CreditCard2$owner)
CreditCard2$selfemp <- to01(CreditCard2$selfemp)

y <- CreditCard2$card

# variables to EXCLUDE

X_data <- subset(CreditCard2, select = -c(card, share, expenditure, months))

# continuous variables to scale
cont_vars <- c("reports", "age", "income", "dependents")

# store means & sds BEFORE scaling
x_means <- sapply(X_data[cont_vars], mean)
x_sds   <- sapply(X_data[cont_vars], sd)

# scale ONLY continuous predictors
X_data[cont_vars] <- scale(X_data[cont_vars])

# build design matrix
X <- model.matrix(~ ., data = X_data)
p <- ncol(X)


#MCMC (Metropolis-within-Gibbs)

set.seed(123)

iters <- 20000
burn  <- 5000

jumper_stds <- rep(0.3, p)

beta_current <- rep(0, p)
samples <- matrix(NA_real_, iters, p)
colnames(samples) <- colnames(X)

accept <- integer(p)
# Iterating across variable, two times as it has multiple variables.
for (i in 1:iters) {
  update_order <- sample(1:p)
  for (j in update_order) {
    beta_prop <- beta_current
    beta_prop[j] <- beta_prop[j] + rnorm(1, 0, jumper_stds[j])
    
    log_ratio <- logposterior_logistic(beta_prop, X, y) -
      logposterior_logistic(beta_current, X, y)
    
    if (log(runif(1)) < log_ratio) {
      beta_current <- beta_prop
      accept[j] <- accept[j] + 1L
    }
  }
  samples[i, ] <- beta_current
}

accept_rate <- accept / iters
print(accept_rate)

post_samples <- samples[(burn + 1):iters, , drop = FALSE]

post_rescaled <- post_samples

# rescale continuous slopes
for (v in cont_vars) {
  post_rescaled[, v] <- post_samples[, v] / x_sds[v]
}

post_rescaled[, "(Intercept)"] <-
  post_samples[, "(Intercept)"] -
  post_samples[, cont_vars] %*% (x_means / x_sds)

post_mean <- colMeans(post_rescaled)
post_sd   <- apply(post_rescaled, 2, sd)

bayes_summary <- data.frame(
  Coefficient = names(post_mean),
  Mean        = post_mean,
  SD          = post_sd,
  row.names   = NULL
)

print(bayes_summary)


