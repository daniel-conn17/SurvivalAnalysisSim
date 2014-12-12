linearAFT <- function(n, p, sd_err, nonzero_beta, target_p_cens) {
  if(length(beta) > p){
    print("exit: beta is larger than p")
    return(NULL)
  }
  if(p > n){
    print("warning: p > n")
  }
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  support_size <- length(beta)
  beta <- c(nonzero_beta,rep(0,p - support_size))
  y <- X%*%beta + rnorm(n,mean=0, sd=sd_err)
  t <- exp(y)
  #generate data to figure censoring parameter
  n_cens <- 10000
  X_cens <- matrix(rnorm(n_cens*p), nrow=n_cens, ncol=p)
  y_cens <- X_cens%*%beta + rnorm(n_cens, mean=0, sd=sd_err)
  t_cens <- exp(y_cens)
  test_quantiles <- seq(.1,.995,.0025)
  test_quant_length <- length(test_quantiles)
  cens_prob <- rep(0, test_quant_length)
  for(i in 1:test_quant_length){
    current_quant <- quantile(t_cens, test_quantiles[i])
    cens <- rexp(n_cens, 1/current_quant)
    cens_prob[i] <- table(cens < t_cens)[2]/n_cens
  }
  diffs <- abs(cens_prob - target_p_cens)
  best_cens_ind <- which.min(diffs)
  cens_param <- quantile(t_cens, test_quantiles[best_cens_ind])
  c <- rexp(n, 1/cens_param)
  delta <- t < c
  delta <- as.numeric(delta)
  t <- apply(cbind(y,c), 1, min)
  return(list(X, y, c, t, delta))
}
