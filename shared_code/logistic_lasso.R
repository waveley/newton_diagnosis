library(tidyverse)

# logistic function
logistic <- function(x) 1 / (1 + exp(-x))

# shrinkage function
S <- function(beta, gamma) {
  if(abs(beta) <= gamma) {
    0
  } else if(beta > 0) {
    beta - gamma
  } else {
    beta + gamma
  }
}

# probability adjustment function
p_adj <- function(p, epsilon) {
  if (p < epsilon) {
    0
  } else if(p > 1 - epsilon) {
    1
  } else {
    p
  }
}

# weight adjustment function
w_adj <- function(p, epsilon) {
  if ((p < epsilon) | (p > 1 - epsilon)) {
    epsilon
  } else {
    p * (1 - p)
  }
}

# executes logistic lasso regression via coordinate descent
logistic_lasso <- function(
  # a numeric design matrix or data frame with named columns
  inputs
  # a vector of outputs; we must have length(output) == nrow(inputs)
  , output
  # a vector of descending penalization factors, ideally on a logarithmic scale
  , lambda_vec
  # standardize inputs using scale
  , standardize   = TRUE
  # a buffer to prevent divergence when fitted probabilities approach 0 or 1
  , epsilon       = 10^-8
  # maximum number of updates to quadratic approximation of likelihood
  , outer_maxiter =  100
  # maximum number of cycles for coordinate descent given quadratic approximation
  , inner_maxiter = 1000
  # tolerance for convergence of coordinate descent
  , tolerance     = 10^-12
) {

  # standardize data unless otherwise specified
  if(standardize) {

    # format data
    X <- as.matrix(cbind(rep(1, nrow(inputs)), scale(inputs)))
    y <- output  

  } else {

    # format data
    X <- as.matrix(cbind(rep(1, nrow(inputs)), inputs))
    y <- output    

  }

  # initialize coefficients at origin
  beta    <- rep(0, ncol(X))
  beta_df <- NULL

  # begin lambda decrement
  for(lambda in lambda_vec) {

    outer_term <- 0
    outer_iter <- 1

    # update quadratic approximation, execute coordinate descent until convergence, repeat
    while(outer_term < 1) {

      # update quadratic approximation; i.e., taylor expand around current estimates  
      p <- map_dbl(logistic(X %*% beta), p_adj, epsilon)
      w <- map_dbl(p, w_adj, epsilon)
      z <- X %*% beta + (y - p) / w

      inner_term <- 0
      inner_iter <- 1

      # given current quadratic approximation, execute coordinate descent
      while(inner_term < 1) {

        beta_old <- beta

        # execute a complete cycle of coordinate descent
        for(k in 1:ncol(X)) {

          # un-penalized coefficient update
          b_k_temp <- sum(w * (z - X[ , -k] %*% beta[-k]) * X[ , k]) / sum(w * X[ , k]^2)
          # shrinkage update
          b_k      <- S(b_k_temp, (k > 1) * lambda / mean(w * X[ , k]^2))
          # update coefficient vector
          beta[k] <- b_k

        }

        inner_iter <- inner_iter + 1

        if(inner_iter == inner_maxiter | max(abs(beta - beta_old)) < tolerance) {

          inner_term <- 1

        }

      }

      outer_iter <- outer_iter + 1

      if(outer_iter == outer_maxiter | inner_iter == 2) {

        outer_term <- 1

      }

    }

    beta_df <- rbind(beta_df, t(beta))

  }

  # format data frame of coefficient estimates
  colnames(beta_df) <- c("intercept", names(inputs))
  beta_df <- as_tibble(beta_df)

  # extract number of variables selected for each lambda
  selected_vec <- apply(beta_df, 1, function(x) sum(x != 0) - 1)

  # output results
  list(lambda = lambda_vec, beta = beta_df, selected = selected_vec)
}