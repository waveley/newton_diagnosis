#################
#
# This function takes in the following parameters:
# > k: number of folds
# > training: training dataset, first column is outcome, remaining columns are predictors.
# > lambda_list: list containing the following:
#     - lambda_start: starting value of lambda vector, without lam_start_stop_func applied
#     - lambda_stop: stopping value of lambda vector, without lam_start_stop_func applied
#     - lambda_step: step size between numbers in linear sequence
#     - lambda_func: transformation of initial sequence of numbers
#     - lambda_start_stop_func: function applied to lambda_start/lambda_stop
#
#
# This function returns a list with two elements:
# > A matrix containing the results of the last CV:
#     - lambda: a column of lambdas
#     - mean_auc: mean AuC
#     - selected_vars: number of selected variables column (excluding intercept)
#     - num_dropped_Vars: number of variables deleted from previous lambda
# > A matrix containing the lambda vectors that were attempted.
#     - lambda_count: id of lambda vector/iteration number
#     - lambda_start: starting value of lambda vector, *without* the lam_start_stop_func
#       applied (i.e., in linear/untransformed units)
#     - lambda_stop: stopping value of lambda vector, *without* the lam_start_stop_func 
#       applied (i.e., in linear/untransformed units)
#     - lambda_step: step size between lambdas
#       if lambda_start > lambda_stop, this should be negative.
#
#################

identity <- function(x){
  return(x)
}

lambda_init <- function(start, stop, step, func = identity){
  lambda_vec <- func(seq(start, stop, step))
  return(lambda_vec)
}

cv_jt <- function(k = 5, training, func, lam_start_stop_func, lambda_list){
  
  lam_start <- lambda_list[[1]]
  lam_stop <- lambda_list[[2]]
  lam_step <- lambda_list[[3]]
  lam_func <- lambda_list[[4]]
  
  lam_list <- tibble(
    lam_count = 0,
    lam_start = 0,
    lam_stop = 0,
    lam_step = 0
  )
  
  lam_count <- 0
  del_too_many_var <- 1
  out_res <- list()
  res <- list()
  
  while (del_too_many_var > 0) {
    # print("working")
    
    lam_count <- lam_count + 1
    
    # saving lambda vector parameters
    cur_lam_list <- tibble(lam_count, lam_start, lam_stop, lam_step)
    
    lam_list <- bind_rows(lam_list, cur_lam_list) 
    
    new_lambda_vec <- lambda_init(lam_start, lam_stop, lam_step, lam_func)
    
    lasso_list <- list()
    auc_list <- list()
    
    pb <- progress_bar$new(
      format = " lasso-ing [:bar] :percent eta: :eta",
      total = 5, clear = FALSE, width = 60)
    
    for (i in 1:k) {
      
      pb$tick()
      
      # this will identify the training set as not i
      trn_set = 
        training %>% 
        filter(fold_id != i) %>% 
        select(-fold_id)
      
      # and this assigns i to be the test set
      tst_set =
        training %>% 
        filter(fold_id == i) %>% 
        select(-fold_id) %>%
        rename(y = diagnosis)
      
      # making matrices
      X_trn <- trn_set[,-1]
      Y_trn <- trn_set$diagnosis
      
      #print("about to lasso")
      
      # lasso_list    
      lasso_list <- func(inputs = X_trn, output = Y_trn, lambda_vec = new_lambda_vec)
      
      #print("done with lasso")
      
      lasso_lambda <- lasso_list[[1]]
      lasso_beta <- lasso_list[[2]]
      lasso_selected <- tibble(selected_num = lasso_list[[3]])
      
      lasso_lam_bet <- cbind(lasso_lambda, lasso_beta) %>% as.matrix()
      
      trn_roc <- auc_calc_lasso(lasso_lam_bet, tst_set)
      auc_list <- bind_rows(auc_list, trn_roc)
      
    }
    
    auc_res <- 
      auc_list %>% 
      group_by(lambda) %>%
      summarise(mean_auc = mean(auc_vals))
    
    res[[lam_count]] <- bind_cols(auc_res, lasso_selected)
    
    res[[lam_count]] <- res[[lam_count]] %>% 
      mutate(num_dropped_vars = selected_num - lag(selected_num, 1))
    
    del_too_many_var <- sum(na.omit(res$num_dropped_vars) > 1)
    
    # del_too_many_var <- 0
    
    if (del_too_many_var > 0) {
      max_auc_lam <- res %>% filter(mean_auc == max(mean_auc)) %>% pull(lambda) %>% mean()
      lam_start <- lam_start_stop_func(max_auc_lam) + 2*abs(lam_step)
      lam_stop <- lam_start_stop_func(max_auc_lam) - 2*abs(lam_step)
      lam_step <- sign(lam_step)*(lam_start - lam_stop)/length(new_lambda_vec)
    }  
  }
  
  # creating dataframe to show lambda values and corresponding mean AUC
  out_res[[1]] <- res
  out_res[[2]] <- lam_list
  
  return(out_res)
}