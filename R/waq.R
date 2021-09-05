#' Estimate treatment effects under additivity by weighting quantile differences
#' 
#' Efficiently estimate the treatment effect under an additive model based on taking a weighted average of the quantile differences:
#'     \deqn{\sum_{i} w(X_{(i)}) (Y_{(i)} - X_{(i)}) / \sum_{i} w(X_{(i)})}
#' where \eqn{w(X_{(i)})} is the (estimated) second derivative of the log density of the control, and subscripts (i) indicate that the vectors are sorted.
#' If X and Y are not of the same length, elements of the shorter vector are duplicated appropriately.
#'
#' @param X numeric vector, the outcomes of the control observations.
#' @param Y numeric vector, the outcomes of the treated observations.
#' @param xf numeric vector, points at which the density is given in d2_logf.
#'          If xf is NULL but d2_logf is supplied, it is assumed that the density was estimated at each points in X.
#' @param d2_logf numeric vector, (estimates of) the second derivative of the log density of the control observations X.
#'          If d2_logf is null, the density (and its derivatives) of X is estimated using estimate_density_d_logs.
#' @param truncate_wf logical, should weights be truncated if (normalized) w/f is too large? (default TRUE)
#' @param return_weights logical, should the weights be returned?
#' @param ... additional arguments passed on estimate_density_d_logs for density estimation.
#' @return list of two elements (three if return_weights is TRUE): \item{tau}{the point estimate of the treatment effect} \item{se}{the estimated standard error} \item{w}{(if return_weights is TRUE) weights, length equal to length of longer of X and Y, in order of sorted X and Y}
#' @examples
#' # draw a random sample with additive treatment effect
#' X <- rexp(n=1000, rate=2)
#' Y <- 0.5 + rexp(n=200, rate=2)
#' waq(X,Y)
#' @export
waq <- function(X,Y,xf=NULL,d2_logf=NULL,truncate_wf=TRUE,return_weights=TRUE,...) {
  # input checks
  if (!is.numeric(X) || !is.numeric(Y)) { stop("X and Y must be numeric") }
  if (!is.logical(truncate_wf)) { stop("truncate_wf must be either TRUE or FALSE") }
  if (!is.logical(return_weights)) { stop("return_weights must be either TRUE or FALSE") }
  
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
  if (!is.logical(truncate_wf)) {
    stop("truncate_wf must be either true or false")
  }
  
  m <- max(length(X),length(Y));  # larger group size
  q = (1:m)/(m+1);  # quantiles at which to evaluate
  X_i <- sort(X[ceiling(q*length(X))])
  Y_i <- sort(Y[ceiling(q*length(Y))])
  w <- stats::approx(x=xf,y=d2_logf,xout=X_i,rule=2)$y
  if (truncate_wf) {
    # normalize weights
    w <- w / sum(w)
    # density f standardized by mad
    f <- pmin(stats::mad(X_i) * (1/m) / diff(X_i), stats::mad(Y_i) * (1/m) / diff(Y_i))
    # threshold
    h <-  log(log(m)) * m^(1/4) / log(m)
    # truncation based on w/f
    w[c(FALSE, w[-1]/f >h)] <- 0
    w[c(w[-m]/f>h, FALSE)]  <- 0
  }
  tau <- sum(w*(Y_i-X_i))/sum(w)
  
  se <- eif_additive_se(X=X,p=length(Y)/(length(X)+length(Y)),
                        xf=xf,d2_logf=d2_logf)
  
  if (return_weights) {
    return(list(tau=tau,
                se=se,
                w=w))
  } else {
    return(list(tau=tau,
                se=se))
  }
}
