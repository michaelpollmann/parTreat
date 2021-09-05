#' Estimate the density (and its derivatives), for sorted input only
#'
#' Estimate the density and its first and second derivatives as input for efficient estimators of the treatment effect.
#' The function implements fast adaptive kernel density estimation by calling a C++ subroutine.
#' For fast computation, the input must be sorted.
#' 
#' @param x0 numeric vector, points at which the density should be evaluated (must be sorted from small to large)
#' @param dat numeric vector, random sample from the distribution for which the density should be estimated.
#' @param estDerivs logical, should the first and second derivative be estimated? (default = TRUE)
#'     The first and second derivative of the density are needed for efficient estimators of the density. But adaptive estimates of the density require an initial estimate of the density (but not its derivatives).
#'     Hence it can be convenient to be able to skip estimating the derivatives for these initial estimates.
#' @param kernel string, indicates which kernel should be used. Currently only "triweight" is implemented.
#' @param sd_dat string or positive numeric, how should the standard deviation be calculated for the Silverman rule-of-thumb bandwidth?
#'     Can be "sd" (the usual standard deviation) or "norm90" (the default) which takes the difference of the 0.95 and 0.05 quantiles and divides by the corresponding range of the normal distribution for an estimate of the standard deviation (if the data is normally distributed) that tends to be substantially less affected by outliers even if the data is not normal.
#'     Alternatively, can be a positive number giving the standard error to be plugged into Silverman's rule-of-thumb bandwidth to manually over / undersmooth.
#' @param adapt logical, should adaptive estimation be used? (default = TRUE) Intuitively, adaptive estimation recognizes that there is less data in the tails of the distribution and hence a larger bandwidth may be helpful in low-density regions than in high-density regions. However, it is computationally more cumbersome.
#' @param fdat numeric vector of the same length as dat, the density evaluated at the points in dat. Used to calculate the adaptive bandwidth if adapt is TRUE. If adapt is TRUE and fdat is not supplied, the function calls itself to estimate the density.
#' @return a matrix with length(x0) rows and 3 columns if estDerivs is TRUE, otherwise 1 column.
#'     The first column contains the estimated density at each point in x0.
#'     The second and third columns contain the estimated first and second derivatives of the density (not log density) at each point in x0, respectively.
#' @examples
#' # draw a random sample from the standard normal distribution and ensure data is sorted
#' X <- sort(rnorm(n=1000))
#' # estimate the density
#' fX <- estimate_density_sorted(X,X,estDerivs=FALSE)
#' # plot the familiar bell curve
#' plot(X,fX)
#' @export
estimate_density_sorted <- function(x0,dat,estDerivs=TRUE,kernel="triweight",sd_dat="norm90",adapt=TRUE,fdat=NULL) {
  # input checks
  if (!is.numeric(x0) || !is.numeric(dat)) { stop("x0 and dat must be numeric") }
  # check that input data is sorted
  if (is.unsorted(x0)) {
    stop("x0 is not sorted")
  }
  if (is.unsorted(dat)) {
    stop("dat is not sorted")
  }
  # check optional arguments
  if (! (kernel %in% c("triweight"))) {
    stop(paste("kernel",kernel,"not implemented"))
  }
  if (! ((sd_dat %in% c("sd","norm90")) | 
         (is.numeric(sd_dat) & (sd_dat > 0)))) {
    stop(paste("sd_data",sd_dat,"not implemented. Must be positive number or one of ",
               paste0("(","sd","norm90",")")))
  }
  if (! is.null(fdat)) {
    if (! (length(fdat)==length(dat))) {
      stop("fdat must have the same length as dat")
    }
    if (any(fdat <= 0)) {
      stop("density fdat must be strictly positive for all points in dat")
    }
  }
  
  # need density if adapt is TRUE
  if ((adapt == TRUE) & is.null(fdat)) {
    fdat <- estimate_density_sorted(dat,dat,estDerivs=FALSE,kernel=kernel,sd_dat=sd_dat,adapt=FALSE)
  }
  
  # Silverman's rule constants for different kernel and order
  # http://www.ssc.wisc.edu/~bhansen/718/NonParametrics1.pdf
  if (kernel == "triweight") {
    K <- K_tri
    silv_const_d0 <- 3.15
    silv_const_d1 <- 2.83
    silv_const_d2 <- 2.70
  }
  
  # estimate standard deviation of data
  if (sd_dat == "sd") {
    s <- stats::sd(dat)
  } else if (sd_dat == "norm90") {
    s <- (stats::quantile(dat,0.95)-stats::quantile(dat,0.05))/(2*stats::qnorm(0.95))
  } else {
    s <- sd_dat
  }
  
  # bandwidth 
  h_d0 <- s * silv_const_d0 / length(dat)^(1/5);
  h_d1 <- s * silv_const_d1 / length(dat)^(1/(5+2));
  h_d2 <- s * silv_const_d2 / length(dat)^(1/(5+4));
  
  # if the density at the data points is supplied, use an adaptive bandwidth
  if (adapt & !is.null(fdat)) {
    # adaptive kernel weights
    alpha <- 0.5
    lambda <- (fdat/exp(mean(log(fdat))))^(-alpha)
    # adjust bandwidth for each observation
    h_d0 <- h_d0 * lambda
    h_d1 <- h_d1 * lambda
    h_d2 <- h_d2 * lambda
  }
  
  
  # allocate output vector of estimated density
  fx <- matrix(nrow=length(x0),ncol=(1+2*estDerivs))
  
  # estimate kernel density
  if (estDerivs) {
    h_pow_2 <- h_d1^2
    h_pow_3 <- h_d2^3
  }
  for (i in seq(length(x0))) {
    # don't recompute density if next point is the same as previous point
    if (i>1) {
      if (x0[i] == x0[i-1]) {
        fx[i,] <- fx[i-1,]
      }
    }
    if (is.na(fx[i,1])) {
      fx[i,1] <- K(x0[i],dat,h_d0,h_d0,0);
      if (estDerivs) {
        fx[i,2] <- K(x0[i],dat,h_d1,h_pow_2,1);
        fx[i,3] <- K(x0[i],dat,h_d2,h_pow_3,2);
      }
    }
  }
  return(fx)
}


#' Estimate the density, its derivatives, and the first and second derivatives of the log density
#'
#' Wrapper around estimate_density_sorted to calculate the first and second derivative of the log of the density.
#' The function estimates the density at all unique points in the input data.
#' It can be useful to estimate the density outside of the treatment effect functions if one computes both the eif and the waq estimator to avoid estimating the density twice, or if one wishes to inspect the density itself.
#' 
#' @param X numeric vector, random sample from the distribution for which the density should be estimated.
#' @param ... additional arguments passed to estimate_density_sorted. DO NOT DISABLE ESTIMATION OF DERIVATIVES.
#' @return a list with four elements:
#'     \item{xf}{numerical vector, the points where the density was estimated (sorted unique points in X)}
#'     \item{f}{matrix with number of rows equal to length(xf) and three columns for the density, its first derivative, and its second derivative}
#'     \item{d1_logf}{(unless estDerivs is set to FALSE) the first derivative of the log density, calculated as f[,2]/f[,1]}
#'     \item{d2_logf}{(unless estDerivs is set to FALSE) the second derivative of the log density, calculated as (f[,1]*f[,3] - f[,2]^2)/f[,1]^2}
#' @examples
#' # draw a random sample with additive treatment effect
#' X <- rexp(n=1000, rate=2)
#' Y <- 0.5 + rexp(n=200, rate=2)
#' # estimate the density and its derivatives
#' dens <- estimate_density_d_logs(X)
#' # pass as argument to avoid computing density inside estimation functions
#' eif_additive(X,Y,xf=dens$xf,d1_logf=dens$d1_logf,d2_logf=dens$d2_logf)
#' waq(X,Y,xf=dens$xf,d2_logf=dens$d2_logf)
#' # calculate and inspect the weights used by waq
#' w_waq <- approx(x=dens$xf,y=dens$d2_logf,xout=sort(X),rule=2)$y
#' # normalize weights and get correct sign
#' w_waq <- w_waq/sum(w_waq)
#' # plot the weights
#' plot(sort(X),w_waq)
#' # plot the weights against quantiles
#' plot(seq(X)/length(X),w_waq)
#' @export
estimate_density_d_logs <- function(X,...) {
  # input checks
  if (!is.numeric(X)) { stop("X") }
  # points at which to estimate the density
  xf <- sort(unique(X))
  # density estimates
  f <- estimate_density_sorted(x0=xf,dat=sort(X),...)
  # only if derivatives were estimated
  if (ncol(f) > 1) {
    # first derivative of log density
    d1_logf <- f[,2]/f[,1]
    # second derivative of log density
    d2_logf <- (f[,1]*f[,3] - f[,2]^2)/f[,1]^2
    return(list(xf = xf,
                f = f,
                d1_logf = d1_logf,
                d2_logf = d2_logf))
  } else {
    return(list(xf = xf,
                f = f))
  }
}
