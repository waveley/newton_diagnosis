######################### AUC for LASSO #########################

#load("./test_matrix.RData")

library(pROC)

auc_calc_lasso <- function(beta_matrix, tst_data){
  # initalizing the storage of the probs 
  pred_list = list()
  
  # looping through all the lambda values 
  for (i in c(1:nrow(beta_matrix))) {
    # getting one lambda's beta values. excluding the lambda value
    beta = beta_matrix[i,-1] 
    
    # getting the x values, do not want the first y value 
    xvals <- cbind(1, tst_data[,-1]) 
    pred = as.matrix(xvals) %*% beta # cross product to get the linear function 
    pred_prob = exp(pred) / (1 + exp(pred)) # link function 
    
    pred_list[[i]] = pred_prob # saving the probabilities 
  }
  
  # putting the probabilities in a data.frame
  pred_tib <- tibble(lambda = beta_matrix[,1] ,pred_list = pred_list, ) 
  
  
  # trying to get AUC
  pred_tib = pred_tib %>% mutate(
    auc_vals = map_dbl( c(1:length(pred_list)), function(x) auc(tst_data$y, pred_list[[x]] %>% as.vector()))
  )
  return(pred_tib)
}


#auc_calc_lasso(beta_matrix, tst_data)
