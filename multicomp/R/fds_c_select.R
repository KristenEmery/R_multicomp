#' @title
#' Select the value of the parameter c for use in multiple decoy competition using the Finite Decoy Storey procedure
#'
#' @description
#' Determines the value of c according the Finite Decoy Storey (FDS) procedure. This procedure sets c equal to the resultant threshold given
#' by the procedure outlined in Storey, Taylor and Seigmund (2004).
#' Can return c for either FDS and FDS1.
#'
#' @param obs_scores_pvals A vector of p-values for the observed scores. Low p-values correspond to more evidence against the null.
#' @param d The number of competing scores for each hypothesis.
#' @param pi_0 The value of \eqn{\pi_0}, the proportion of true null hypotheses.
#' @param max_c The maximum value of c that can be obtained. Default is 0.95.
#' @param add Either 0 or 1. Determines if we add 1/(d+1) to the c chosen by the procedure. Set to 1 if using FDS1.
#' @param q The value of \eqn{\alpha}, the desired FDR threshold for the Storey test.
#'
#' @return
#' \item{c}{The selected value of c.}
#'
#' @author Kristen Emery
#'
#' @references
#' Storey JD, Taylor JE, and Siegmund D. (2004) Strong control,
#' conservative point estimation, and simultaneous conservative
#' consistency of false discovery rates: A unified approach. Journal of
#' the Royal Statistical Society, Series B, 66: 187-205. \cr
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2004.00439.x/abstract}
#'
#' @references
#' Emery K, Keich U, Hasam S and Nobel W. (2019) Multiple competition based FDR control. arXiv:1907.01458\cr
#' \url{https://arxiv.org/abs/1907.01458}
#'
#' @export

fds_c_select <- function(obs_scores_pvals, d, pi_0, max_c = 0.95, add = 0, q) {

  n <- length(obs_scores_pvals)
  possible_c <- seq(1/(d + 1), max_c, 1/(d + 1))

  FDR_est <- rep(0, d)
  for (i in 1:d) {
    # Run through each possible value of c (except 1)
    if (i %in% 1:length(possible_c)) {
      FDR_est[i] <- (i/(d + 1) * pi_0 * n)/sum(obs_scores_pvals <= i/(d + 1))  #Estimate FDR if we rejected all pvals less than i/(n_p + 1)
    } else {
      FDR_est[i] <- 1
    }
  }
  c <- suppressWarnings(max(which(FDR_est <= q))/(d + 1)) + add/(d + 1)
  if (c <= 0) {
    c <- 1/(d + 1)
  } else if (c > max_c) {
    c <- possible_c[max(can_c)]
  }

  return(c)
  # The final choice of c is the maximum c that satifies the Storey condition.
}
