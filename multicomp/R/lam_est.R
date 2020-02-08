#' @title
#' Estimate the parameter \eqn{\lambda} for multiple decoy competion.
#'
#' @description
#' Determines the value of \eqn{\lambda} for use in multiple competition. This value is obtained by tabulating the p-values and
#' sequentially using a binomial test to find the point at which the remaining number of p-values to the right are evenly distributed around their mid-point.
#'
#' @param obs_scores_pvals Vector of p-values for the observed scores. Low p-values correspond to more evidence against the null.
#' @param d The number of competing scores for each hypothesis.
#' @param q The significance level we perform the binomial test at. This test determines whether we classify the p-values as 'evenly distributed'.

#' @return
#' \item{lam}{The selected value of \eqn{\lambda}.}
#'
#' @author Kristen Emery
#'
#' @references
#' Emery K, Keich U, Hasam S and Nobel W. (2019) Multiple competition based FDR control. arXiv:1907.01458\cr
#' \url{https://arxiv.org/abs/1907.01458}
#'
#' @export

lam_est <- function(obs_scores_pvals, d, q) {
  # Estimates lambda using a binomial test on the data.

  fin <- 0
  i <- 1
  while (fin == 0 && i < floor(0.95 * (d + 1))) {
    # Count how many p-values are equal to each possible value
    pval_count <- rep(0, d + 1)
    for (j in 1:(d + 1)) {
      pval_count[j] <- sum(obs_scores_pvals == j/(d + 1))
    }
    removed_count <- pval_count[-(1:i)]  #Remove current values

    if (length(removed_count)%%2 > 0) {
      # In the event of an uneven split we remove the middle value
      gpA <- sum(removed_count[1:floor(length(removed_count)/2)])  #Number in first group
      n_remain <- sum(removed_count[-ceiling(length(removed_count)/2)])  #Total remaining observations
    } else {
      gpA <- sum(removed_count[1:(length(removed_count)/2)])
      n_remain <- sum(removed_count)
    }

    if (n_remain <= 0) {
      fin <- 1
    } else {
      p <- binom.test(gpA, n_remain, p = 0.5, alternative = "greater")$p.value  #Test if the remaining scores are evenly distributed across the groups

      if (p > q) {
        fin <- 1  #Stop if binomial test is not siginificant
      } else {
        i <- i + 1
      }
    }
  }

  lam <- i/(d + 1)

  return(lam)
}
