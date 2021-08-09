#' Translate estimated multiplicative (in logs) effects to level effects (+ standard errors)
#' 
#' After estimating treatment effects with the outcomes in logs, this helper function translates the effect back into a level effect.
#' It also calculates the appropriate standard error for the in-sample effect using the delta method.
#' 
#' If one assumes a multiplicative model \eqn{Y(1) = \theta Y(0)}, the specification is additive after taking the natural logarithm.
#' Hence, one can use the function eif_additive or waq applied to outcomes in logs to estimate log\eqn{(\theta)}.
#' However, one is often substantively interested in the average effect of the treatment rather than in the multiplicative parameter \eqn{\theta}.
#' This function therefore translates the effect into the sample average level effect for a sample of control observations X_pop and treated observations Y_pop.
#' 
#' In the multiplicative model, even if \eqn{\theta} is known exactly, estimates of the population average effect, under the multiplicative model
#'     \deqn{Y(1) - Y(0) = (\theta-1) Y(0)}
#' can be noisy because they require estimating the population mean of the control (e.g. by the sample mean).
#' The standard errors here are therefore calculate for a fixed sample of control X and treated Y.
#' Under the assumption of the multiplicative model, outcomes can then be considered fixed, and the only variation comes from estimating \eqn{\theta}
#'
#' @param est_log a list of (at least) the two elements tau and se, for instance as returned by eif_additive or waq. It is assumed that the estimate is the (untransformed) estimate of an additive specification with the outcome variables in logs.
#' @param X_pop numeric vector, the outcomes of the control observations in levels, or alternatively of a super-population of interest.
#' @param Y_pop numeric vector, the outcomes of the treated observations in levels, or alternatively of a super-population of interest.
#' @return list of four elements: 
#'         \item{tau_level}{the point estimate of the treatment effect in levels}
#'         \item{se_level}{the standard error for the treatment effect in levels}
#'         \item{tau_log}{the point estimate of the treatment effect with the outcome in logs (copied from the input list, est_log$tau)}
#'         \item{se_log}{the standard error for the treatment effect with the outcome in logs (copied from the input list, est_log$se)}
#' @examples
#' # draw a random sample with multiplicative treatment effect
#' X <- rexp(n=1000, rate=2)
#' Y <- 1.1 * rexp(n=200, rate=2)
#' # estimate the treatment effect with an additive specification in logs, either waq or eif
#' est_waq_log <- waq(log(X),log(Y))
#' est_eif_log <- eif_additive(log(X),log(Y))
#' # translate the effects into level effects
#' log_to_level(est_waq_log,X,Y)
#' log_to_level(est_eif_log,X,Y)
#' @export
log_to_level <- function(est_log, X_pop, Y_pop) {
  
  mu_X <- mean(X_pop)
  mu_Y <- mean(Y_pop)
  n0 <- length(X_pop)
  n1 <- length(Y_pop)
  
  tau_log <- est_log$tau
  theta <- exp(tau_log)
  se_log <- est_log$se
  
  # average treatment effect estimate
  tau_level <- (n0 * (theta * mu_X - mu_X) + n1 * (mu_Y - mu_Y/theta)) / (n0 + n1)
  
  # translate standard errors
  se_level <- sqrt((n0*theta*mu_X + n1*mu_Y/theta)/(n0+n1) * (se_log^2) * (n0*theta*mu_X + n1*mu_Y/theta)/(n0+n1))
  
  return(list(tau_level = tau_level,
              se_level = se_level,
              tau_log = tau_log,
              se_log = se_log))
}