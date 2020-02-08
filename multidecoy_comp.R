
#' @title
#' Multiple competition with parameter selection
#'
#' @description
#' Performs multiple competition on a given set of observed scores and corresponding competing scores and returns the discovery list.
#' Selects c and \eqn{\lambda} by a user input method. Can test a range of FDR significance levels at once.
#'
#' @details
#' This is a function that combines parameter selection with mulitple competition using a \code{\link{mirandom}} mapping scheme (which it calls internally once parameters have been selected) to obtain the discovery list.
#' Depending on the parameter method chosen (fds, fds1) \code{\link{lam_est}}, \code{\link{pi_0_est_St}} and \code{\link{fds_c_select}} are called internally to calculate \eqn{\lambda}, \eqn{\pi_0} and c repectively.
#' If lbm is the selected parameter method then \code{\link{LBM_method_select}} is called to determine the optimal method of selecting c and \eqn{\lambda} for each \eqn{\alpha}.
#' The procedure will then apply the chosen method to acquire the discoveries. The parameters \code{bs_runs}, \code{test_methods} and \code{fallback} are not used unless \code{c_lam_selection = 'lbm'} and can be left as default.
#'
#' @param obs_scores A vector corresponding to the observed scores.
#' @param decoy_scores A \eqn{d * m} matrix of competiting null scores. Each column corresponds to a specific hypothesis/observed score. Each row corresponds to a single competing set.
#' @param alpha_range A vector of \eqn{\alpha}, the FDR thresholds, that we wish to test at.
#' @param c_lam_selection Method of choosing the parameters c and \eqn{\lambda}.
#'  Possible options: \code{'mirror'} - mirror test (\eqn{c = \lambda = 0.5}), \code{'lf'} - Lei-Fithian choice (\eqn{c = \alpha, \lambda = 0.5}), \code{'fds'} - FDS, \code{'fds1'} - FDS1, \code{'lbm'} - LBM.
#' @param LBM_bs_runs LBM parameter: A value that sets the number of bootstrap iterations to perform.
#' @param LBM_test_methods LBM parameter: A vector of parameter selecting methods that we choose between (see \code{c_lam_selection} for options). The order indicates the priority when we break ties. Default is FDS > mirror > FDS1.
#' @param LBM_fallback LBM parameter: The name of the method that we fall back upon if the estimated maximum FDR is greater than \eqn{\alpha} for our bootstrap samples (see \code{c_lam_selection} for options). Default is FDS1.
#'
#' @return
#' Returns a list made up of:
#' \item{Discoveries}{A list showing the indicies of the discoveries. Each entry in the list corresponds to the discovery list for a specific value of \eqn{\alpha} in alpha_range. A response of 0 corresponds to no discoveries.}
#' \item{c}{A vector showing the value of c that was used for each \eqn{\alpha}.}
#' \item{lambda}{A vector showing the value of \eqn{\lambda} that was used for each \eqn{\alpha}.}
#'
#' @author Kristen Emery
#'
#' @references
#' Emery K, Keich U, Hasam S and Nobel W. (2019) Multiple competition based FDR control. arXiv:1907.01458\cr
#' \url{https://arxiv.org/abs/1907.01458}
#'
#' @examples
#' m = 200; k = 20; d = 5; alpha_range = c(0.01,0.05,0.1,0.2)
#' null_mu = rep(0, m); mu = null_mu; sig = rep(2, m)
#' alts = sample(m, k)
#' mu[alts] = 5
#' obs_scores = rnorm(m,mu,sig)
#' #Rows correspond to decoy sets, columns correspond to hypotheses.
#' decoy_scores = matrix(rnorm(m*d,null_mu,sig), ncol = m, nrow = d, byrow = TRUE)
#'
#' #Basic usage of multidecoy_comp with the LBM method:
#' results = multidecoy_comp(obs_scores, decoy_scores, alpha_range, c_lam_selection = 'lbm')
#' print(results)
#'
#' @export
#'

multidecoy_comp <- function(obs_scores, decoy_scores, alpha_range, c_lam_selection, LBM_bs_runs = 50, LBM_test_methods = c("fds", "mirror", "fds1"),
                            LBM_fallback = "fds1") {

  nqs <- length(alpha_range)

  scores <- rbind(obs_scores, decoy_scores)
  d <- nrow(scores) - 1

  rank_scores <- apply(scores, 2, rank, ties.method = "random")  #Rank each target-decoy set.
  obs_scores_pvals <- (d + 2 - rank_scores[1, ])/(d + 1)  #Generate empirical p-values

  c_choice <- rep(0, nqs)
  lam_choice <- rep(0, nqs)

  Discoveries <- as.list(rep(0, nqs))

  if (c_lam_selection == "mirror") {

    c_choice <- rep(0.5, nqs)
    lam_choice <- rep(0.5, nqs)

    for (iq in 1:nqs) {
      res_c_mirror <-mirandom(scores, 0.5, 0.5, alpha_range[iq], rank_scores = rank_scores)
      Discoveries[[iq]] <- res_c_mirror$Discoveries_ind
    }

  } else if (c_lam_selection == "lf") {

    # Uses c = min(alpha, 0.5), lambda = 0.5 If c = alpha is not exactly obtainable then selects max(c <= alpha)

    lat_vals <- (1:d)/(d + 1)  #Possible values for c and lambda

    lam_choice <- rep(0.5, nqs)

    for (iq in 1:nqs) {
      c_FL <- suppressWarnings(max(lat_vals[lat_vals - 1e-12 < alpha_range[iq]]))
      if (c_FL < min(lat_vals)) {
        c_FL <- min(lat_vals)
      }

      c_choice[iq] <- min(c_FL, 0.5)

      res_c_lf <-mirandom(scores, min(c_FL, 0.5), 0.5, alpha_range[iq], rank_scores = rank_scores)

      Discoveries[[iq]] <- res_c_lf$Discoveries_ind
    }

  } else if (c_lam_selection == "fds") {

    # Choose c as the maximum value of the emprical p-values of the observed scores rejected by Storey's FDR controlling procedure Use lambda =
    # estimated lambda If c > lambda. Set c = min(c, lambda)

    # Estimate lambda and pi_0
    lam <- lam_est(obs_scores_pvals, d, 0.1)
    pi_0_est1 <- pi_0_est_St(obs_scores_pvals, lam, add = 1)

    lam_choice <- rep(lam, nqs)

    for (iq in 1:nqs) {
      c_pvals <- fds_c_select(obs_scores_pvals, d, pi_0_est1, lam, add = 0, alpha_range[iq])

      c_choice[iq] <- min(c_pvals, lam)

      res_c_fds <-mirandom(scores, min(c_pvals, lam), lam, alpha_range[iq], rank_scores = rank_scores)

      Discoveries[[iq]] <- res_c_fds$Discoveries_ind
    }

  } else if (c_lam_selection == "fds1") {

    # Choose c as the maximum value of the emprical p-values of the observed scores rejected by Storey's FDR controlling procedure + 1/(d + 1)
    # Use lambda = estimated lambda If c > lambda. Set lambda = max(c, lambda)

    # Estimate lambda and pi_0
    lam <- lam_est(obs_scores_pvals, d, 0.1)
    pi_0_est <- pi_0_est_St(obs_scores_pvals, lam, add = 0)

    for (iq in 1:nqs) {
      c_pvals <- fds_c_select(obs_scores_pvals, d, pi_0_est, 0.95, add = 1, alpha_range[iq])

      c_choice[iq] <- c_pvals
      lam_choice[iq] <- max(c_pvals, lam)

      res_c_fds1 <-mirandom(scores, c_pvals, max(c_pvals, lam), alpha_range[iq], rank_scores = rank_scores)

      Discoveries[[iq]] <- res_c_fds1$Discoveries_ind
    }


  } else if (c_lam_selection == "lbm") {

    lam <- lam_est(obs_scores_pvals, d, 0.1)

    pi_0_est <- pi_0_est_St(obs_scores_pvals, lam, add = 0)
    pi_0_est1 <- pi_0_est_St(obs_scores_pvals, lam, add = 1)

    # Perform bootstrap method choice
    res_bs <- LBM_method_select(scores, alpha_range, bs_runs = LBM_bs_runs, test_methods = LBM_test_methods, fallback = LBM_fallback, input_c = lam,
                                rank_scores = rank_scores, pi_0 = pi_0_est)

    method_choice <- res_bs$m_choice

    for (iq in 1:nqs) {
      if (method_choice[iq] == "mirror") {

        lam_choice[iq] <- 0.5
        c_choice[iq] <- 0.5

        res_c_mirror <-mirandom(scores, 0.5, 0.5, alpha_range[iq], rank_scores = rank_scores)
        Discoveries[[iq]] <- res_c_mirror$Discoveries_ind

      } else if (method_choice[iq] == "lf") {

        lat_vals <- (1:d)/(d + 1)  #Possible values for c and lambda

        lam_choice[iq] <- 0.5

        c_FL <- suppressWarnings(max(lat_vals[lat_vals < alpha_range[iq]]))
        if (c_FL < min(lat_vals)) {
          c_FL <- min(lat_vals)
        }

        c_choice[iq] <- min(c_FL, 0.5)

        res_c_lf <-mirandom(scores, min(c_FL, 0.5), 0.5, alpha_range[iq], rank_scores = rank_scores)
        Discoveries[[iq]] <- res_c_lf$Discoveries_ind

      } else if (method_choice[iq] == "fds") {

        c_pvals <- fds_c_select(obs_scores_pvals, d, pi_0_est1, lam, add = 0, alpha_range[iq])

        lam_choice[iq] <- lam
        c_choice[iq] <- min(c_pvals, lam)

        res_c_fds <-mirandom(scores, min(c_pvals, lam), lam, alpha_range[iq], rank_scores = rank_scores)

        Discoveries[[iq]] <- res_c_fds$Discoveries_ind

      } else if (method_choice[iq] == "fds1") {

        c_pvals <- fds_c_select(obs_scores_pvals, d, pi_0_est, 0.95, add = 1, alpha_range[iq])

        lam_choice[iq] <- max(c_pvals, lam)
        c_choice[iq] <- c_pvals

        res_c_fds1 <-mirandom(scores, c_pvals, max(c_pvals, lam), alpha_range[iq], rank_scores = rank_scores)

        Discoveries[[iq]] <- res_c_fds1$Discoveries_ind
      }

      # Monotonicity check
      if (iq > 1) {
        if (length(Discoveries[[iq]]) < length(Discoveries[[iq - 1]])) {
          method_choice[iq] <- method_choice[iq - 1]
          if (method_choice[iq] == "mirror") {

            lam_choice[iq] <- 0.5
            c_choice[iq] <- 0.5

            res_c_mirror <-mirandom(scores, 0.5, 0.5, alpha_range[iq], rank_scores = rank_scores)
            Discoveries[[iq]] <- res_c_mirror$Discoveries_ind

          } else if (method_choice[iq] == "lf") {

            lat_vals <- (1:d)/(d + 1)  #Possible values for c and lambda

            lam_choice[iq] <- 0.5

            c_FL <- suppressWarnings(max(lat_vals[lat_vals < alpha_range[iq]]))
            if (c_FL < min(lat_vals)) {
              c_FL <- min(lat_vals)
            }

            c_choice[iq] <- min(c_FL, 0.5)

            res_c_lf <-mirandom(scores, min(c_FL, 0.5), 0.5, alpha_range[iq], rank_scores = rank_scores)
            Discoveries[[iq]] <- res_c_lf$Discoveries_ind

          } else if (method_choice[iq] == "fds") {

            c_pvals <- fds_c_select(obs_scores_pvals, d, pi_0_est1, lam, add = 0, alpha_range[iq])

            lam_choice[iq] <- lam
            c_choice[iq] <- min(c_pvals, lam)

            res_c_fds <-mirandom(scores, min(c_pvals, lam), lam, alpha_range[iq], rank_scores = rank_scores)

            Discoveries[[iq]] <- res_c_fds$Discoveries_ind

          } else if (method_choice[iq] == "fds1") {

            c_pvals <- fds_c_select(obs_scores_pvals, d, pi_0_est, 0.95, add = 1, alpha_range[iq])

            lam_choice[iq] <- max(c_pvals, lam)
            c_choice[iq] <- c_pvals

            res_c_fds1 <-mirandom(scores, c_pvals, max(c_pvals, lam), alpha_range[iq], rank_scores = rank_scores)

            Discoveries[[iq]] <- res_c_fds1$Discoveries_ind

          }
        }
      }
    }
  }

  return(list(Discoveries = Discoveries, c = c_choice, lambda = lam_choice))

}
