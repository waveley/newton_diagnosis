# ##########
#
# This file contains the lasso function.
# The lasso function takes in:
# > datx: dataset of covariates. all non-predictor variables
#         should already have been removed.
# > daty: vector of outcomes. should already be binary (0 or 1).
# > lambda_vec: vector of lambdas to iterate over
#
# The function returns:
# > beta_matrix: matrix of lambdas and beta coefficients.
#                first column of matrix are lambdas.
#
# ##########

epsilon <- 10^(-5)

lasso <- function(x, y, lambda_vec){
  
    X <- x %>% mutate(intercept = 1) %>% as.matrix()
    Y <- y
  
    source("shared_code/lasso_helper.R")  
  
    # initialize parameters
    beta <- rep(0, ncol(X))
    beta_matrix <- t(c(NA, beta))
    
    # pb <- progress_bar$new(
    #  format = " lasso-ing [:bar] :percent eta: :eta",
    #  total = length(lambda_vec), clear = FALSE, width = 60)
    
    for (lambda in lambda_vec) {
    #  pb$tick()
      outer_term <- 0
      outer <- 1
      while (outer_term < 1) {
        
        terminate <- 0
        iter <- 1
        while (terminate < 1) {
          
          beta_old <- beta
          
          p <- map_dbl(logistic(X %*% beta), p_adj, epsilon)
          w <- map_dbl(p, w_adj, epsilon)
          z <- X %*% beta + (Y - p) / w
          
          for (k in 1:ncol(X)) {
            x_k    <- X[ , k]
            x_notk <- X[ , -k]
            b_notk <- beta[-k]
            
            # un-penalized coefficient update
            b_k_temp <- sum(w * (z - x_notk %*% b_notk) * x_k) / sum(w * x_k^2)
            
            # shrinkage update
            b_k      <- S(b_k_temp, lambda * (k > 1) / mean(w * x_k^2))
            
            # update beta vector along with other parameters
            beta[k] <- b_k
          }
          
          iter <- iter + 1
          
          if (iter == 100 | max(abs(beta - beta_old)) < 10^(-10)) {
            terminate <- 1
          }
          
        }
        
        outer <- outer + 1
        
        # if outer == 100 or iter == 2, quit while loop and add to beta_matrix
        if (outer == 100 | iter == 2) {
          outer_term <- 1
        }
        
      }
      
      beta_matrix <- rbind(beta_matrix, t(c(lambda, beta)))
      
    }
    
    return(beta_matrix)

}

#### testing function ###


#n    <- nrow(data)
#X    <- scale(data[ , -c(1, 2)])
#X    <- data[ , -c(1, 2)]
#X    <- as.matrix(cbind(rep(1, n), X))
#y    <- data$diagnosis

#beta <- rep(0, ncol(X))
#lambda_vec <- exp(seq(log(0.4), log(0.4/1000), -0.2))
#lambda_vec <- rep(seq(10, 2, -1), 4) / c(rep(1, 9), rep(10, 9), rep(100, 9), rep(1000, 9))
#lambda_vec <- exp(seq(log(0.4), log(0.4/10), -0.1))
#lambda_vec <- seq(0.40, 0.20, -0.01)

#lasso(X, y, lambda_vec)

#beep()
























