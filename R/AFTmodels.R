#' Generate data from linear log-normal AFT model.
#'
#' @param  n Sample size.
#' @param  p Number of covariates.
#' @param  sd_err Standard deviation of the error term, \eqn{\epsilon}.
#' @param  support_beta Support of \eqn{\beta} (non-zero coefficients).
#' @param  censoring_rate Censoring rate.
#' @return linearAFT returns an object of class "survsim".  An object of class
#'        "survsim" is a list containing the following components.
#'    \itemize{
#'    \item X: matrix of covariates
#'    \item y: unobserved log-survival times
#'    \item t: unobserved survival times
#'    \item c: unobserved censoring times.
#'    \item obs_times: observed time on study (min(t,c)).
#'    \item delta: vector indicating whether observation is
#'    censored (1=failure, 0=censored).
#'    }
#'
#' @examples linear_aft(1000, 5, .1, 1, .3)


linear_aft <- function(n, p, sd_err, support_beta, censoring_rate) {
  if(length(support_beta) > p){
    print("exit: support of beta is larger than p")
    return(NULL)
  }
  if(p > n){
    print("warning: p > n")
  }
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  support_size <- length(support_beta)
  beta <- c(support_beta,rep(0,p - support_size))
  y <- X%*%beta + rnorm(n,mean=0, sd=sd_err)
  t <- exp(y)
  #Censoring distribution will be exponential.
  #We will generate a larger data set to determine
  #what it's mean should be.
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
  diffs <- abs(cens_prob - censoring_rate)
  best_cens_ind <- which.min(diffs)
  cens_param <- quantile(t_cens, test_quantiles[best_cens_ind])
  c <- rexp(n, 1/cens_param)
  delta <- t < c
  delta <- as.numeric(delta)
  obs_times <- apply(cbind(t,c), 1, min)
  out <- list(X, y, t, c, obs_times, delta)
  class(out) <- "survsim"
  return(out)
}
