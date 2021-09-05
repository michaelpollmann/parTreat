#' Estimate treatment effects under additivity using the influence function
#'
#' Efficiently estimate the treatment effect under an additive model based on the influence function.
#' Finds the parameter \eqn{\theta} such that 
#'     \deqn{\int f_{0}'(x)/f_{0}(x) f_{0}(x) dx = \int f_{0}'(y-\theta)/f_{0}(y-\theta) f_{1}(y) d y (*)}
#' and estimates the standard error of \eqn{\hat{\theta}} (under correct specification).
#' Beware of warnings (about failures to find a solution) coming from pracma::fzero. These are currently not handled properly.
#' 
#'
#' @param X numeric vector, the outcomes of the control observations.
#' @param Y numeric vector, the outcomes of the treated observations.
#' @param xf numeric vector, points at which the density is given in d1_logf (and d2_logf).
#'          If xf is NULL but d1_logf is supplied, it is assumed that the density was estimated at each points in X.
#' @param d1_logf numeric vector, (estimates of) the first derivative of the log density of the control observations X.
#'          If d1_logf is null, the density (and its derivatives) of X is estimated using estimate_density_d_logs.
#' @param d2_logf numeric vector, (estimates of) the second derivative of the log density of the control observations X.
#'          Used to estimate standard errors.
#'          If d2_logf is null but calc_se is TRUE (default), the density (and its derivatives) of X is estimated using estimate_density_d_logs.
#' @param theta_init numeric, initial value to find the solution to (*). Defaults to the difference in medians.
#' @param calc_se logical, should standard errors be calculated? (default = TRUE)
#'          This is slower if d1_logf but not d2_logf is provided as then the density needs to be estimated for standard errors but not for the point estimate.
#' @param ... additional arguments passed on estimate_density_d_logs for density estimation.
#' @return list of two elements if calc_se = TRUE or one element otherwise: 
#'         \item{tau}{the point estimate of the treatment effect}
#'         \item{se}{(if calc_se = TRUE) the estimated standard error}
#' @examples
#' # draw a random sample with additive treatment effect
#' X <- rexp(n=1000, rate=2)
#' Y <- 0.5 + rexp(n=200, rate=2)
#' eif_additive(X,Y)
#' @export
eif_additive <- function(X, Y, xf=NULL, d1_logf=NULL, d2_logf=NULL, theta_init=NULL,calc_se=TRUE,...) {
  # input checks
  if (!is.numeric(X) || !is.numeric(Y)) { stop("X and Y must be numeric") }
  if (!is.logical(calc_se)) { stop("calc_se must be either TRUE or FALSE") }
  
  # compute default estimates of optional arguments if necessary
  if (is.null(theta_init)) {
    theta1 <- stats::median(Y) - stats::median(X)
  } else if (!is.numeric(theta_init)) {
    stop("theta_init must be NULL or a numeric")
  } else if (length(theta_init) > 1) {
    stop("theta_init must be length 1")
  }
  if (is.null(d1_logf)) {
    dens <- estimate_density_d_logs(X=X,estDerivs=T,...)
    xf <- dens$xf
    d1_logf <- dens$d1_logf
    d2_logf <- dens$d2_logf
  } else if (is.null(xf)) {
    if (length(xf)==length(X)) {
      xf <- X
    } else {
      stop("derivative of log density (d1_logf) supplied but not the points at which the density was evaluated")
    }
  }
  if (!(length(xf) == length(d1_logf))) {
    stop("lengths of estimated density derivatives (d1_logf) and points at which the density was estimated (xf) differ")
  }
  
  # integral over derivative of log density at X
  mean_control <- mean(stats::approx(x=xf,y=d1_logf,xout=X, rule=2)$y)
  
  # integral over derivative of log density at transformed Y
  fun <- function(theta) { 
    return(mean_control -
             mean(stats::approx(x=xf,y=d1_logf,xout=Y-theta, rule=2)$y))
  }
  
  theta1 <- pracma::fzero(fun, theta1)$x
  
  if (calc_se) {
    se <- eif_additive_se(X=X,p=length(Y)/(length(X)+length(Y)),
                          xf=xf,d2_logf=d2_logf)
    return(list(tau=theta1,
                se=se))
  } else {
    return(list(tau=theta1))
  }
}

#' Estimate the standard error of the efficient estimators under additivity
#'
#' The standard errors can be used for both eif_additive and waq under correct specification of the additive model.
#' 
#' @param X numeric vector, the outcomes of the control observations. (Or a separate sample of observations from the same distribution as the control observations.)
#' @param p numeric between 0 and 1, the fraction of the sample that is treated.
#' @param xf numeric vector, points at which the density is given in d2_logf.
#'          If xf is NULL but d2_logf is supplied, it is assumed that the density was estimated at each points in X.
#' @param d2_logf numeric vector, (estimates of) the second derivative of the log density of the control observations X.
#'          If d2_logf is null, the density (and its derivatives) of X is estimated using estimate_density_d_logs.
#' @param n integer, the sample size. Defaults to length(X)/(1-p).
#' @param ... additional arguments passed on estimate_density_d_logs for density estimation.
#' @return the estimated standard error of the efficient additive treatment effect estimators.
#' @examples
#' # draw a random sample of observations
#' X <- rexp(n=1000, rate=2)
#' # standard error assuming half of the estimation sample is control
#' eif_additive_se(X=X,p=0.5)
#' # standard error assuming there are 2000 control and 200 treated in the estimation sample
#' eif_additive_se(X=X,p=200/(2000+200),n=2000+200)
#' 
#' @export
eif_additive_se <- function(X,p,xf=NULL,d2_logf=NULL,n=NULL,...) {
  # input checks
  if (!is.numeric(X)) { stop("X must be numeric") }
  if (! ((p > 0) & (p < 1))) {
    stop("treatment probability p must be strictly between 0 and 1")
  }
  
  # compute default estimates of optional arguments if necessary
  if (is.null(d2_logf)) {
    dens <- estimate_density_d_logs(X=X,estDerivs=T,...)
    xf <- dens$xf
    d2_logf <- dens$d2_logf
  } else if (is.null(xf)) {
    if (length(xf)==length(X)) {
      xf <- X
    } else {
      stop("derivative of log density (d2_logf) supplied but not the points at which the density was evaluated")
    }
  }
  if (!(length(xf) == length(d2_logf))) {
    stop("lengths of estimated density derivatives (d2_logf) and points at which the density was estimated (xf) differ")
  }
  if (is.null(n)) {
    n <- length(X)/(1-p)
  }
  
  se <- sqrt(1/(p*(1-p)*(-1)*mean(stats::approx(x=xf,y=d2_logf,xout=X,rule=2)$y)))/sqrt(n)
  return(se)
}
