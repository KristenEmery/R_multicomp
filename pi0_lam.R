#' @title
#' Estimate the proportion of true null hypotheses using Storey's estimation
#'
#' @description
#' Estimates the proportion of true nulls in a vector of p-values using the estimation provided by Storey (2002) and Storey, Taylor & Seigmund (2004). Uses either a pre-set value
#' of \eqn{\lambda} or selects one based off the data.
#' Can be set to use a correction for finite hypotheses.
#'
#' @param p Vector of p-values for the observed scores. A low p-value corresponds to more evidence against the null.
#' @param d The number of decoys.
#' @param add Either 0 or 1. Determines what is added to the numerator in the calculation. 0 corresponds to Storey's asymptotic formulation, while 1 corresponds to the finite correction.
#' @param lambda The value of the tuning parameter in the estimate. Either a number between 0 and 1 (non-inclusive) or \code{NULL} with the latter corresponding to using
#' the binomial test method to select \eqn{\lambda} according to the data.
#' @param lam_sig The significance level of the binomial test to choose \eqn{\lambda}.
#'
#' @return
#' \item{pi_0_est}{The estimated value of \eqn{\pi_0}.}
#' \item{lambda}{The value of lambda used to calculated \eqn{\pi_0}.}
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
#' Emery K, Keich U, Hasam S and Nobel W. (2019) Multiple competition based FDR control. arXiv:1907.01458\cr
#' \url{https://arxiv.org/abs/1907.01458}
#'
#' @export


pi0est_mc = function(p, d, add, lambda = NULL, lam_sig = 0.1){

  n <- length(p)
  lam_length <- length(lambda)

  if(lam_length > 1){
    stop("The length of lambda should be 1 (a number between 0 and 1) or 0 (empty).")
  }

  if(is.numeric(lambda)){

    if(lam <=0 | lam >= 1){
      stop("Lambda must be between 0 and 1 (non-inlcusive).")
    }

    pi_0 <- min((sum(p > lambda) + 1 * add)/((1 - lambda) * n), 1)

  } else {

    lambda <- lam_est(p, d, lam_sig)
    pi_0 <- min((sum(p > lambda) + 1 * add)/((1 - lambda) * n), 1)

  }
  if(pi_0 < 0){
    stop("pi_0 < 0")
  }

  return(list(pi_0 = pi_0, lambda = lambda))

}

