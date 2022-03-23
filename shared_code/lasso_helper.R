##################
#
# This file contains four lasso helper functions:
# > logistic: returns logit of input x
# > p_adj: probability adjustment function.
#          pulls input probability to 0 or 1 if 
#          big/small enough (determined by epsilon).
#          otherwise, returns input probability.
# > w_adj: weight adjustment function
#          sends weight to epsilon if input probability is 
#          big/small enough (determined by epsilon)
#          otherwise, return p*(1-p).
#
##################

# logistic function
logistic <- function(x) 1 / (1 + exp(-x))

# shrinkage function
S <- function(beta, gamma) {
  if (abs(beta) <= gamma) {
    0
  } else if (beta > 0) {
    beta - gamma
  } else {
    beta + gamma
  }
}

# probability adjustment function
p_adj <- function(p, epsilon) {
  if (p < epsilon) {
    0
  } else if (p > 1 - epsilon) {
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