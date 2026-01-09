#Data initialization

az <- na.omit(azcabgptca)

# response (count)
y <- as.integer(az$los)

X_data <- subset(az, select = -c(los))

cont_vars <- c("age")

x_means <- sapply(X_data[cont_vars, drop = FALSE], mean)
x_sds   <- sapply(X_data[cont_vars, drop = FALSE], sd)

X_data[cont_vars] <- scale(X_data[cont_vars, drop = FALSE])

# design matrix intercept
X <- model.matrix(~ ., data = X_data)
p <- ncol(X)


#Log-posterior: Poisson + Normal priors

logposterior_poisson <- function(beta, X, y, s0 = 10, s = 2.5) {
  eta <- as.vector(X %*% beta)
  # prevents exp overflow
  if (any(eta > 50)) return(-Inf)
  loglik <- sum(y * eta - exp(eta))
  if (!is.finite(loglik)) return(-Inf)
  #priors
  logprior <- dnorm(beta[1], 0, s0, log = TRUE) + sum(dnorm(beta[-1], 0, s, log = TRUE))
  loglik + logprior
}

# Initialize at Poisson GLM
glm_fit <- glm(y ~ ., data = data.frame(y = y, X_data), family = poisson())
beta_current <- coef(glm_fit)[colnames(X)]  # align order with X columns




set.seed(123)
iters <- 20000
burn  <- 5000
jumper_stds <- rep(0.05, p)
names(jumper_stds) <- colnames(X)

samples <- matrix(NA_real_, iters, p, dimnames = list(NULL, colnames(X)))
accept  <- integer(p)

lp_current <- logposterior_poisson(beta_current, X, y)
#Loop the variables
for (i in 1:iters) {
  for (j in sample.int(p)) {
    
    beta_prop <- beta_current
    beta_prop[j] <- beta_prop[j] + rnorm(1, 0, jumper_stds[j])
    
    lp_prop <- logposterior_poisson(beta_prop, X, y)
    log_ratio <- lp_prop - lp_current
    
    if (log(runif(1)) < log_ratio) {
      beta_current <- beta_prop
      lp_current <- lp_prop
      accept[j] <- accept[j] + 1L
    }
  }
  samples[i, ] <- beta_current
}

accept_rate <- accept / iters
print(accept_rate)

post_samples <- samples[(burn + 1):iters, , drop = FALSE]

#Rescale coefficients back to original units

post_rescaled <- post_samples

# rescale continuous slopes
for (v in cont_vars) {
  post_rescaled[, v] <- post_samples[, v] / x_sds[v]
}

# adjust intercept for centering/scaling
post_rescaled[, "(Intercept)"] <-
  post_samples[, "(Intercept)"] -
  post_samples[, cont_vars, drop = FALSE] %*% (x_means / x_sds)

bayes_summary <- data.frame(
  Coefficient = colnames(post_rescaled),
  Mean = colMeans(post_rescaled),
  SD   = apply(post_rescaled, 2, sd),
  row.names = NULL
)

print(bayes_summary)


# MCMC diagnostics for coefficient


# choose coefficient, change name if you want other
coef_name <- "procedure"
coef_idx  <- which(colnames(X) == coef_name)

burnin <- 2000
coef_sample <- samples[(burnin + 1):iters, coef_idx]
par(mfrow = c(2, 2))

# Trace plot
plot(coef_sample,
     type = "l",
     main = paste("Trace plot –", coef_name),
     xlab = "Iteration",
     ylab = expression(beta))

#Autocorrelation 
acf(coef_sample,
    main = paste("ACF –", coef_name))

#Histogram
hist(coef_sample,
     breaks = 40,
     freq = FALSE,
     main = paste("Posterior distribution –", coef_name),
     xlab = expression(beta))

#Density 
lines(density(coef_sample),
      col = "red",
      lwd = 2)

#Posterior mean 
abline(v = mean(coef_sample), col = "blue", lwd = 2)
abline(v = quantile(coef_sample, c(0.025, 0.975)),
       col = "darkgreen", lty = 2, lwd = 2)


# Plot
par(mfrow = c(1, 1))

