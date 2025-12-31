# Fisher's scoring

# Author: Federico Boiocchi

# the following is a unique function that performs Fisher'scoring
# for binary regression (link logit and probit)
# and poisson regression (link log).

# In the case of canonical link, the Fisher's scoring algorithm boils down
# to a Newton-Raphson applied to the log-likelihood of the sample (whose mean is modeled
# with a GLM). Therefore no computations of the expected fisher information
# matrix which requires integration is required. In practice the Fisher's
# scoring has been implemented using its IRLS (iteratively reweighted least squares)
# formulation. (also in the probit non-canonical link case the computations are easy
# since the expectation of the negative hessian removes a quantity depending on y leading
# to a fisher info which is the same as the canonical link fisher info matrix)

# The IRLS algorithm is a Weighted Leasted Square (WLS) applied to a pseudoresponse
# z that is updated at each iteration using the current estimate of the betas
# and using the current matrix of weights W which depends also on beta


rm(list = ls())

# arguments:
# - y: response (binary or count)
# - X: design matrix (no intercept included and no response included)
# - family: string denoting the DE1 distribution to use
# - link: string denoting the link function used (implemented are: logit,probit,log)
# - trace: boolean, if TRUE the whole sequence of estimates is produced with the index of current iteration
# - intercept: boolean, if TRUE an intercept is included
# - eps: tolerance under which convergence is declared L2 norm of beta_new-beta_old < eps
# - maxit: maximum number of iterations after which convergence is declared even though the tolerance level is not reached


fs <- function(y, X, family, link, trace, intercept, eps, maxit) {
  # response variable
  y <- as.matrix(y)

  # design matrix (no intercept)
  X <- as.matrix(X)

  col_names <- colnames(X)

  # sample size
  n <- length(y)

  # number of variables not including the intercept
  p <- ncol(X)
  nc <- p

  # starting value
  b1 <- rep(0, nc)

  if (intercept) {
    X <- cbind(rep(1, n), X)
    b1 <- rep(0, p + 1)
    nc <- p + 1
  }

  # temporary initial distance (to enter in the while loop)
  dist <- eps + 1

  counter <- 1

  while (dist > eps) {
    if (counter > maxit) {
      break
    }

    # linear predictor
    eta <- as.vector(X %*% b1)

    if (family == "poisson") {
      mu <- exp(eta)
      w <- mu
      z <- eta + (y - mu) / w
    }

    if (family == "bernoulli") {
      if (link == "logit") {
        mu <- 1 / (1 + exp(-eta))
        w <- mu * (1 - mu)
        z <- eta + (y - mu) / w
      }
      if (link == "probit") {
        mu <- pnorm(eta, mean = 0, sd = 1)
        vmu <- mu * (1 - mu)
        w <- (dnorm(eta)^2) / vmu
        z <- eta + (y - mu) / dnorm(eta)
      }
    }

    # Weight matrix
    W <- diag(w)

    # updating IRLS equation with Tikhonov regularization to avoid sigularities
    b2 <- solve((t(X) %*% W %*% X + diag(rep(1e-15, nc))), t(X) %*% W %*% z)

    # L2 norm distance (used for the stopping rule)
    dist <- sqrt(sum((b2 - b1)^2))
    b1 <- b2

    if (trace) {
      out <- t(b1)
      out <- c(b1, as.integer(counter))
      if (intercept) {
        names(out) <- c("(intercept)", col_names, "n_iter")
        cat("\nCoefficients:\n")
        print(out)
      } else {
        names(out) <- c(col_names, "n_iter")
        cat("\nCoefficients:\n")
        print(out)
      }
    }
    counter <- counter + 1
  }
  out <- as.numeric(t(b1))
  return(out)
}


# Test with data

# azcabgptca data

library("COUNT")
data("azcabgptca")
data <- azcabgptca
mod <- glm(los ~ ., family = poisson(link = "log"), data = data)
y <- data$los
X <- data[, which(names(data) != "los")]

# poisson regression
fs(y, X, family = "poisson", trace = TRUE, intercept = TRUE, eps = 1e-15, maxit = 100)
coef(mod)

# CreditCard data
library(AER)
data("CreditCard")
data <- CreditCard
y <- data$card
X <- data[, !names(data) %in% c("share", "expenditure", "card")]

mod_logit <- glm(card ~ . - share - expenditure, family = binomial(link = "logit"), data = data)
mod_probit <- glm(card ~ . - share - expenditure, family = binomial(link = "probit"), data = data)

y <- ifelse(data$card == "yes", 1, 0)
X$selfemp <- ifelse(CreditCard$selfemp == "yes", 1, 0)
X$owner <- ifelse(CreditCard$owner == "yes", 1, 0)

# logistic regression
fs(y, X, family = "bernoulli", link = "logit", trace = TRUE, intercept = TRUE, eps = 1e-15, maxit = 100)
coef(mod_logit)

# probit regression
fs(y, X, family = "bernoulli", link = "probit", trace = TRUE, intercept = TRUE, eps = 1e-15, maxit = 100)
coef(mod_probit)
