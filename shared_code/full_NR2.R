#########################################################################
## This file contains the function to run a full Newton Raphson Model  ##
##
## The code in this script has been adapted from newton_fullmod_clean.Rmd
## This script should be used in conjuction with full_nr_lasso.Rmd
##
## Model Input: 
##  - dat : this is the data vector. Note that diagnosis need to be y 
##          and y needs to take values 0 and 1. 
##  - start : this is the inital values for the beta vectors 
##            please note that the values should be small 0.01 or less
##  - tol :  1e-8, 
##  - maxiter : 200
##
## Model Outputs: 
##  - "i" : what iteratioon we are on 
##  - "loglik" : logliklihood value 
##  - "step" : To see the half stepping work 
##  - Coef : These values are named by the variables put into the model. 
##           There will be as many covariates specified + one intercept 
##
## Using Function 
##   To get Just the coef see comment at the bottom 
## 



# Here we have the log-likelihood function that produces 
# the log-likelihood, gradient vector, and hessian matrix, given a dataset and beta vector
loglike_func <- function(dat, betavec){
  
  dat = bc_trn
  
  # x matrix
  dat_temp <-
    dat %>%
    mutate(intercept = 1) %>%
    select(-diagnosis) %>%
    relocate(intercept)
  
  dat_x <-  
    dat_temp %>%
    as.matrix() %>%
    unname() 
  
  # pi vector
  u <- dat_x %*% betavec
  pi <- exp(u) / (1 + exp(u))
  
  # loglikelihood
  loglik <- sum(dat[,1]*u - log(1 + exp(u)))
  
  #gradient 
  grad <- t(dat_x) %*% (dat[,1] - pi)
  
  # Hessian 
  W <- diag(nrow(pi))
  diag(W) <- pi*(1 - pi)
  hess <- -(t(dat_x) %*% W %*% (dat_x))
  
  return(list(loglik = loglik, grad = grad, hess = hess))
  
}

# Next, we have a Newton Raphson algorithm to determine beta coefficients that maximize the likelihood of the function.
NewtonRaphson <- function(dat,  start, tol = 1e-8, maxiter = 200){
  
  i <- 0
  cur <- start 
  stuff <- loglike_func(dat, cur)
  res <- c(i = 0, "loglik" = stuff$loglik,  "step" = 1, cur)
  
  prevloglik <- -Inf # to make sure it iterates
  
  while (i < maxiter && abs(stuff$loglik - prevloglik) > tol) {
    step <- 1
    i <- i + 1
    prevloglik <- stuff$loglik
    
    # check negative definite  
    eigen_vals <- eigen(stuff$hess)
    
    if (max(eigen_vals$values) <= 0 ) { # check neg def, if not change 
      hess <- stuff$hess
    } else { # if it is pos def then need to adjust 
      hess <- stuff$hess - (max(eigen_vals$values) + 0.1)*diag(nrow(stuff$hess))
    } 
    
    prev  <- cur
    cur   <- prev - step*(solve(stuff$hess) %*% stuff$grad)
    stuff <- loglike_func(dat, cur) # log-lik, gradient, Hessian
    
    # step halving
    while (stuff$loglik < prevloglik) {
      stuff <- loglike_func(dat, prev)
      step  <- step / 2 # this is where half steping happens 
      cur   <- prev - step*(solve(stuff$hess) %*% stuff$grad)
      stuff <- loglike_func(dat, cur)
    }
    # add current values to results matrix
    res <- rbind(res, c(i, stuff$loglik, step, cur))
  }
  
  colnames(res) <- c("i", "loglik",  "step", "intercept", names(dat[,-1]))
  return(res)
}

# getting our best estimates
betavec <- c(rep(0.01, ncol(bc_trn))) %>% as.matrix()
nr1 <- NewtonRaphson(bc_trn, betavec)
nr_beta_est <- nr1[nrow(nr1), -c(1:3)] %>% as.vector()