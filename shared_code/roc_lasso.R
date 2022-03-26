#################
#
# This function takes in the following parameters:
# > beta estimates: a vector of beta coefficients
# > test data: the heretofore un-touched data for final comparisons
#
# This function returns a ROC object, and can be changed to return just the AUC.
# However, the AUC can be extracted from the ROC object by roc_object$auc. Additionally,
# the ROC plot can be easily plotted by plot(roc_object, legacy.axes = TRUE).
#
#################

roc_func <- function(beta_est, test_data){
  # we have this flexible in case we want to test fewer variables 
  # terms <- beta_est %>% pull(term) 
  # col.num <- which(colnames(test_data) %in% terms)
  # select the desired x values 
  test_data = bc_tst
  xvals <- test_data[,-1] %>% 
    mutate(
      intercept = 1 # create a intercept variable 
    ) %>%
    relocate(intercept) # move it to the front 
  pred <- as.matrix(xvals) %*% beta_est # get the cross product of the linear model
  logit_pred <- exp(pred) / (1 + exp(pred)) # link function to get probabilities
  
  #auc_val <- auc(test_data[,1], as.vector(logit_pred)) # calculating the AUC
  
  roc_res <- roc(test_data[,1], as.vector(logit_pred)) # output ROC
  
  return(roc_res)
}