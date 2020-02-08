#' @title
#' Estimate the proportion of true null hypotheses using Storey's estimation
#'
#' @description
#' Estimates the proportion of true nulls in a vector of p-values using the estimation provided by Storey (2002) and Storey, Taylor & Seigmund (2004).
#' Can be set to use a correction for finite hypotheses.
#'
#' @param obs_scores_pvals Vector of p-values for the observed scores. A low p-value corresponds to more evidence against the null.
#' @param lambda The value of the tuning parameter in the estimate.
#' @param add Either 0 or 1. Determines what is added to the numerator in the calculation. 0 corresponds to Storey's asymptotic formulation, while 1 corresponds to the finite correction.
#'
#' @return
#' \item{pi_0_est}{The estimated value of \eqn{\pi_0}.}
#'
#' @author Kristen Emery
#'
#' @references
#' Storey JD. (2002) A direct approach to false discovery rates. Journal
#' of the Royal Statistical Society, Series B, 64: 479-498. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00346/abstract}
#'
#' Storey JD, Taylor JE, and Siegmund D. (2004) Strong control,
#' conservative point estimation, and simultaneous conservative
#' consistency of false discovery rates: A unified approach. Journal of
#' the Royal Statistical Society, Series B, 66: 187-205. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2004.00439.x/abstract}
#'
#' @export

pi_0_est_St <- function(obs_scores_pvals, lambda, add = c(0, 1)) {

  n <- length(obs_scores_pvals)

  pi_0_est <- min((sum(obs_scores_pvals > lambda) + 1 * add)/((1 - lambda) * n), 1)

  return(pi_0_est)

}
