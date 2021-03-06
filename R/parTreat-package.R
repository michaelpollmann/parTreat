#' parTreat: Efficiently Estimate Treatment Effects Based On A Parametric Model For Treatment Effects
#'
#' Efficiently estimate treatment effects imposing a parametric assumption on treatment effects 
#'   of the form \eqn{Y(1) = h(Y(0),\theta)}.
#' Currently, two estimators for the additive model \eqn{Y(1) = Y(0) + \theta} are implemented,
#'   as well as a helper function to use them to estimate treatment effects based on a multiplicative model.
#' 
#' @section Usage:
#' 
#' If X are the outcomes of the control observations, and Y are the outcomes of the treated observations, call
#' 
#' \code{eif_additive(X,Y)}
#' 
#' or
#' 
#' \code{waq(X,Y)}
#' 
#' to efficiently estimate the treatment effect assuming additivity.
#' 
#' For a multiplicative model, use
#' 
#' \code{est_eif_log <- eif_additive(log(X),log(Y))}
#' 
#' \code{log_to_level(est_eif_log, X, Y)}
#' 
#' or
#' 
#' \code{est_waq_log <- waq(log(X),log(Y))}
#' 
#' \code{log_to_level(est_waq_log, X, Y)}
#'
#' Under correct specification (of the additive or multiplicative model),
#'   the eif and waq estimators are asymptotically equivalent.
#' They are adaptive, in the sense that if the additive model is correctly specified,
#'   they are just as efficient as parametric estimators
#'   that use knowledge of the true distribution of the outcomes (up to shift).
#' 
#' @section Note:
#' 
#' The functions estimate the density based on the control observations.
#' This is convenient if the control group is larger than (or at least as large as) the treatment group.
#' If the treated group is larger, one can call the functions with reverse groups Y and X.
#' The treatment effect is then minus the estimated effect, and no adjustment is necessary for the standard error.
#' 
#' Generally
#' \code{waq(X,Y)$tau}
#' and
#' \code{-waq(Y,X)$tau}
#' are not the same.
#' Under correct specification, they are asymptotically equivalent,
#' but in finite samples they will differ based on which group is used to estimate the densities.
#' The same is also true when using \code{eif_additive}.
#' 
#' This package ports the Matlab code used in the simulations of Athey et al. (2021) to R.
#' Please use the Matlab replication code to exactly replicate the results
#'   and for implementations of other estimators studied by Athey et al. (2021).
#'
#' @section Reference:
#' 
#' Susan Athey, Peter Bickel, Aiyou Chen, Guido Imbens, and Michael Pollmann. \emph{Semiparametric Estimation of Treatment Effects in Randomized Experiments}. 2021.
#' 
#' @docType package
#' @name parTreat
#' @useDynLib parTreat, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL
#> NULL
