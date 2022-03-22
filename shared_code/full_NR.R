########################################################################
## This file contains the function to run a full Newton Rapson Model  ##
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

loglike_func <- function(dat, betavec){
  # setting up an intercept 
  dat_temp = dat %>%
    rename(intercept = y) %>%
    mutate(intercept = rep(1, nrow(dat) ))
  dat_x = unname(as.matrix(dat_temp)) # creating the x matrix 
  
  # finding the pi values 
  u = dat_x %*% betavec
  pi <- exp(u) / (1+exp(u))
  
  # loglikelihood
  loglik <- sum(dat$y*u - log(1 + exp(u)))
  
  #gradient 
  grad <- t(dat_x)%*%(dat$y - pi)
  
  # Hessian 
  W = matrix(0, nrow = dim(dat)[1], ncol = dim(dat)[1]) 
  diag(W)= pi*(1-pi) 
  hess = -(t(dat_x)%*% W %*% (dat_x))
  
  return(list(loglik = loglik, grad = grad, hess = hess))
  
}
#loglike_func(dat, betavec = c( rep(0.03, 31)))  # test! 


NewtonRaphson_a <- function(dat,  start, tol = 1e-8, maxiter = 200){
  i = 0
  cur = start 
  stuff = loglike_func(dat, cur)
  res <- c(i=0, "loglik" = stuff$loglik,  "step" = 1, cur)
  prevloglik <- -Inf # To make sure it iterates
  while(i < maxiter && abs(stuff$loglik - prevloglik) > tol) {
    step = 1
    i <- i + 1
    prevloglik <- stuff$loglik
    
    # checking negative definite  
    eigen_vals = eigen(stuff$hess)
    if(max(eigen_vals$values) <= 0 ){ # check neg def, if not change 
      hess = stuff$hess
    } else{ # if it is pos def then need to adjust 
      hess = stuff$hess - (max(eigen_vals$values) + 0.1)*diag(nrow(stuff$hess))
    } 
    
    
    prev  <- cur
    cur   <- prev - rep(step, length(prev))*(solve(stuff$hess) %*% stuff$grad)
    stuff <- loglike_func(dat, cur) # log-lik, gradient, Hessian
    
    # doing half stepping 
    while(stuff$loglik < prevloglik){
      stuff <- loglike_func(dat, prev)
      step  <- step / 2 # this is where half steping happens 
      cur   <- prev - step*(solve(stuff$hess) %*% stuff$grad)
      stuff <- loglike_func(dat, cur)
    }
    # Add current values to results matrix
    res <- rbind(res, c(i, stuff$loglik, step, cur))
  }
  
  colnames(res) <- c("i", "loglik",  "step", "intercept", names(dat[,-1]))
  return(res)
}



betavec = c(rep(0.01, ncol(trn_data))) # beta vec example
ans <- NewtonRaphson_a(trn_data, betavec) # running method using trn_data dataset

# get just the beta values
est_ans <- data.frame(ans) %>%
  select(-step, -loglik, -i) %>%
  filter(row_number() == n()) %>%
  pivot_longer(
    cols = everything(),
    names_to = "term",
    values_to = "vals"
  )


predicted <- predict(glmnet_fit, newx = tst_data %>% select(-y), type="response")

#calculate AUC
library(pROC)
auc(test$default, predicted)

