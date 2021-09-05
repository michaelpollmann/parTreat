#' Calculate the difference of asymmetric trimmed means between Y and X
#' 
#' Calculate the difference in asymmetric trimmed means between Y and X.
#' The percent trimmed from left (\eqn{\alpha}) and right ((\eqn{\beta})) are chosen adaptively using grid search.
#' If gridsize = 'sqrtn', the grid of \eqn{\alpha} and \eqn{\beta} values are spaced by 1/sqrt(n) where n is the size of the smaller group.
#' If gridsize = 'n', the grid values are spaced by 1/n.
#' Note that the two-dimensional grid search is then of size n^2 which can be prohibitively slow for large n, while sqrt(n) will be near optimal and feasible even for large samples.
#' Alternatively, specify any integer of at least 3, up to the size of the smaller group.
#' 
#' The parameters max_alpha and max_beta allow to set bounds for the amount of trimming from the left and right, respectively.
#' Note that the grid spacing is not affected by this; that is, fewer values of \eqn{\alpha} or \eqn{\beta} need to be searched over when changing the defaults.
#' To only consider trimming of up to the smallest 3% and largest 10% of the data, set max_alpha = 0.03 and max_beta = 0.1.
#' To only trim the right tail of the data, set max_alpha = 0; forcing the function to not trim the left tail may improve stability when appropriate, and should be much faster computationally.
#' The size of the sample after trimming is 1 - \eqn{\alpha} - \eqn{\beta} where \eqn{\alpha} and \eqn{\beta} are the selected trimming percentages.
#' 
#' 
#' @param X numeric vector, the outcomes of the control observations.
#' @param Y numeric vector, the outcomes of the treated observations.
#' @param gridsize either the string 'sqrtn' or 'n' or an integer giving the number of alpha and beta values to try 
#' @param max_alpha numeric between 0 and 1, the maximum fraction of observations to trim from the left (default 1)
#' @param max_beta numeric between 0 and 1, the maximum fraction of observations to trim from the right (default 1)
#' @param min_frac numeric between 0 and 1, the minimum fraction of observations to KEEP after trimming (default 0.25)
#' @param min_obs numeric, the minimum number of observations of each treatment state to KEEP after trimming (default 50)
#' @return list of four elements:
#'         \item{tau}{the point estimate of the treatment effect}
#'         \item{se}{the estimated standard error (if calc_se = TRUE)}
#'         \item{alpha}{the trimming fraction from the left}
#'         \item{beta}{the trimming fraction from the right}
#' @examples
#' # draw a random sample with additive treatment effect
#' X <- rnorm(n=1000)
#' Y <- rnorm(n=1000)
#' atm_diff(X,Y)
#' # for smaller samples, a n-grid may be computationally feasible
#' X <- rnorm(n=100)
#' Y <- rnorm(n=100)
#' atm_diff(X,Y, gridsize='n')
#' @export
atm_diff <- function(X, Y, gridsize="sqrtn", max_alpha = 1, max_beta = 1, min_frac = 0.25, min_obs = 50) {
  # input checks
  if (max_alpha < 0 || max_alpha > 1 || !is.numeric(max_alpha)) { stop("max_alpha must be between 0 and 1")}
  if (max_beta < 0 || max_beta > 1 || !is.numeric(max_beta)) { stop("max_beta must be between 0 and 1")}
  if (min_frac < 0 || min_frac > 1 || !is.numeric(min_frac)) { stop("min_frac must be between 0 and 1")}
  if (!is.numeric(min_obs)) { stop("min_obs must be an integer") }
  if (min_obs >= min(length(X),length(Y))) { stop("min_obs must be less than the number of observations of either group") }
  
  # sort observations
  xs <- if (is.unsorted(X)) { sort(X) } else X
  ys <- if (is.unsorted(Y)) { sort(Y) } else Y
  
  # use grid search
  n <- min(length(xs), length(ys))
  if (is.character(gridsize)) {
    if (gridsize == "sqrtn") {
      gridsize <- ceiling(sqrt(n))
    } else if (gridsize == "n") {
      gridsize <- n
    } else {
      stop("unknown character supplied for parameter gridsize, should be either one of c('sqrtn','n') or an integer")
    }
  } else if (is.numeric(gridsize)) {
    if (gridsize > n) { stop(paste("parameter gridsize cannot be larger than smaller group, which has size",n)) }
    if (gridsize < 3) { stop("gridsize of less than 3 points not allowed")}
    gridsize <- ceiling(gridsize)
  } else {
    stop("unknown type for parameter gridsize, should be either one of c('sqrtn','n') or an integer")
  }
  grids <- utils::combn(seq(0, 1, length.out = gridsize), 2)
  left_bounds <- grids[1,]
  right_bounds <- grids[2,]
  # impose minimum fraction of observations
  lr_allowed <- (right_bounds - left_bounds) >= min_frac
  left_bounds <- left_bounds[lr_allowed]
  right_bounds <- right_bounds[lr_allowed]
  # impose max_alpha
  right_bounds <- right_bounds[left_bounds <= max_alpha]
  left_bounds <- left_bounds[left_bounds <= max_alpha]
  # impose max_beta
  left_bounds <- left_bounds[1 - right_bounds <= max_beta]
  right_bounds <- right_bounds[1 - right_bounds <= max_beta]
  
  VecATrimmedMean <- function(sortedx, leftbs, rightbs) {
    # assume the x vector is sorted
    # returns the variances of asymmetric trimmed means for paired (left, right) bounds
    nx <- length(sortedx)
    L <- pmin(nx, 1+floor(nx * leftbs))  # location of left trimming
    R <- ceiling(nx  * rightbs)  # location of right trimming
    S <- cumsum(sortedx)   # first moment
    Q <- cumsum(sortedx^2)  # second  moment
    # trimmed mean
    vecmu_trimmed <- (S[R] - ifelse(L>1, S[L-1], 0)) / (R-L+1)
    # winsorized mean (for variance formula)
    vecmu <- ((L-1)*sortedx[L] + (nx-R) * sortedx[R]  +  (S[R] - ifelse(L>1, S[L-1], 0)))  / nx
    # variance formula
    vecsigma2  <-  ((L-1)*sortedx[L]^2 + (nx-R) * sortedx[R]^2  +
                      (Q[R] - ifelse(L>1,  Q[L-1], 0)) -  nx * vecmu^2)
    # don't allow L>=R
    vecsigma2[L>=R] <- Inf
    # enforce minimum number of observations
    vecsigma2[L>=R-ceiling(min_obs)+1] <- Inf
    return(list(mu = vecmu_trimmed,
                sigma2 = vecsigma2 / (R - L + 1)^2))
  }
  
  # calculate asymmetric trimmed mean at all grid points
  atm_x <- VecATrimmedMean(xs, left_bounds, right_bounds)
  atm_y <- VecATrimmedMean(ys, left_bounds, right_bounds)
  
  # pick variance minimizing trimming parameters
  sigma2 <-  atm_x$sigma2 + atm_y$sigma2
  kopt <- which.min(sigma2)
  alpha <- left_bounds[kopt]
  beta <- 1 - right_bounds[kopt]
  
  # trimmed mean with these parameters
  ymu <- atm_y$mu[kopt]
  xmu <- atm_x$mu[kopt]
  
  return(list(tau=ymu - xmu,
              se=sqrt(sigma2[kopt]),
              alpha=alpha,
              beta=beta
  ))
}