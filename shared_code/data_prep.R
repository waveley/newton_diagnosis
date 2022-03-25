#################
#
# This file provides a way to standardize data import/initial variable selection across programs.
# This code: 
# > imports the raw data
# > converts diagnosis to binary,
# > deselects id
# > cuts out some predictors that are correlated with others
# > creates training and test datasets
#
################

source("./shared_code/partition.R")

bc <- 
  read_csv("./data/breast-cancer.csv") %>% 
  mutate(diagnosis = 1 * (diagnosis == "M")) %>% 
  select(-id)

remove_bad_vars <- function(indat, bad_vars){
  
  outdat <-
    indat %>% 
    dplyr::select(-bad_vars)
    return(outdat)

}

bad <- c(
  "area_mean", "area_worst", "perimeter_mean", "perimeter_worst", "radius_mean"
  , "perimeter_se", "area_se"
  , "concave points_worst", "concavity_mean"
  , "texture_worst"
)


bc_trunc <- remove_bad_vars(bc, bad)

part_bc <- partition(p = 0.8, data = bc_trunc)

bc_trn <- 
  part_bc %>% 
  filter(part_id == "train") %>% 
  select(-part_id)

bc_tst <- 
  part_bc %>% 
  filter(part_id == "test") %>% 
  select(-part_id)





















