## Creating CV Folds or Sets ##
# this function should be used on training data

cv_sets <- function(k = 5, training){
  # again generating a probability value
  fold_p = runif(nrow(training), min = 0, max = 1)
  
  split_list = list()
  # looping to establish cutoff values for each fold
  for (i in 1:k){
    split_list[[i]] = data.frame("cutoff" = i/k, "fold" = i)
  }
  split_df = as.data.frame(do.call(rbind, split_list))
  # this is some rather manual code to assign each observation to a fold
  # if this function needed to handle varying values of k, this would need to change
  fold_id = case_when(fold_p <= split_df[1,1] ~ split_df[1,2],
                      fold_p <= split_df[2,1] & fold_p > split_df[1,1] ~ split_df[2,2],
                      fold_p <= split_df[3,1] & fold_p > split_df[2,1] ~ split_df[3,2],
                      fold_p <= split_df[4,1] & fold_p > split_df[3,1] ~ split_df[4,2],
                      fold_p <= split_df[5,1] & fold_p > split_df[4,1] ~ split_df[5,2])
  
  data_new = cbind(training, fold_p, fold_id)
  
  return(data_new)
}