#' @title
#' Labeled Bootstrap Maximisation
#'
#' @description
#' Performs Labeled Bootstrap Maximisation (LBM) on a given set of observed and competing scores to select to optimal method of selecting the parameters c and
#' \eqn{\lambda} for use in multiple competition. Can test a range of FDR significance levels.
#'
#' @param scores A \eqn{(d+1) * m} matrix of combined observed and competing scores. Each column corresponds to a indivdual hypothesis while each row corresponds to an observed or competing score set. The first row is the observed scores while the subsequent rows correspond to specific competing scores.
#' @param alpha_range A vector of \eqn{\alpha}, the desired FDR significance levels.
#' @param bs_runs The number of bootstrap samples to generate.
#' @param test_methods A vector of the names of the parameter selecting methods that we choose between. The order indicates the priority when we break ties.
#' Possible options: \code{'mirror'} - mirror test (\eqn{c = \lambda = 0.5}), \code{'lf'} - Lei-Fithian choice (\eqn{c = \alpha, \lambda = 0.5}), \code{'fds'} - FDS, \code{'fds1'} - FDS1, \code{'lbm'} - LBM. Default is FDS > mirror > FDS1.
#' @param fallback The name of the method that we fall back upon if the estimated maximum FDR is greater than \eqn{\alpha}. Default is FDS1.
#' @param input_c Input value of the parameter c that is used in the labeling heuristic. Determines how we select a target win. Lower value results in a stricter condition.
#' @param max_c The maximum allowed value of c in the LF, FDS and FDS1 methods. By default this is set of 0.95.
#' @param rank_scores A \eqn{(d+1) * m} matrix of ranked observed and competing scores. Same structure as the \code{scores} parameter but with the score ranked across observed-competing sets. Default is NULL which will result in the function caluclating the ranks, but can be manually input.
#' @param pi_0 An estimate of pi_0. NULL by default in which case it will be caluclated in the function, but can otherwise be manually set.
#'
#' @details This function internally calls \code{\link{label_heuristic}} to estimate null/alternative labels for each bootstrap samples as well as \code{\link{mirandom}}
#' to perform multiple-decoy competition. Depending on the chosen \code{test_methods}
#' this function will call \code{\link{lam_est}}, \code{\link{pi_0_est_St}} and \code{\link{fds_c_select}} to determine c and \eqn{\lambda}.
#'
#' @return
#' Returns a list:
#' \item{m_choice}{A vector indicating the selected method for each specific \eqn{\alpha}.}
#' \item{Fallback}{A vector indicating if the fallback was triggered for each specific \eqn{\alpha}.}
#'
#' @author Kristen Emery
#' @references
#' Emery K, Keich U, Hasam S and Nobel W. (2019) Multiple competition based FDR control. arXiv:1907.01458\cr
#' \url{https://arxiv.org/abs/1907.01458}
#'
#' @export


LBM_method_select <- function(scores, alpha_range, bs_runs = 50, test_methods = c("fds", "mirror", "fds1"), fallback = "fds1", input_c = 0.5,
                              max_c = 0.95, rank_scores = NULL, pi_0 = NULL) {
  # Hybrid bootstrap to attain a vector of c for each alpha given in alpha_range. We assign the scores as null or alt and resample the indicies
  # of the observed scores then the decoys of the nulls. Then for each bootstrap run we acquire the FDP of the best run across all the c If the
  # average of the FDR is over alpha then we select c = flk, otherwise we select c equal to the best c over the average of the bootstrap runs.

  n <- ncol(scores)  #Number of hypotheses
  d <- nrow(scores) - 1  #Number of decoys
  n_m <- length(test_methods)  #Number of methods to select between
  nqs <- length(alpha_range)  #Number of alphas to test over.
  possible_c <- seq(1/(d + 1), max_c, 1/(d + 1))

  if (is.null(rank_scores)) {
    rank_scores <- apply(scores, 2, rank, ties.method = "random")  #Rank each target-decoy set.
  }

  if (is.null(pi_0)) {
    obs_scores_pvals <- (d + 2 - rank_scores[1, ])/(d + 1)
    lam <- lam_est(obs_scores_pvals, d, 0.1)
    pi_0 <- pi_0_est_St(obs_scores_pvals, lam, add = 0)
  }

  n_target_discs <- array(0, dim = c(bs_runs, length(alpha_range), length(test_methods)))
  FDR_est <- n_target_discs

  # We estimate nulls and alternatives using a segmenting heuristic. Resample indicies of the observed scores then resample the decoys of the
  # nulls.
  for (ibs in 1:bs_runs) {

    estTrueNullInds <- label_heuristic(scores, input_c, rank_scores = rank_scores)  #Estimate null labels
    resample_indicies <- sample(1:n, n, replace = TRUE)
    est_nulls_bs <- which(estTrueNullInds[resample_indicies] == 1)  #Which of the bootstrap resampled scores are null/alt

    scores_bs <- scores[, resample_indicies]

    # Permute null labelled scores and their decoys
    if (length(est_nulls_bs) > 0) {
      null_score_bs <- as.matrix(scores_bs[, est_nulls_bs])
      null_score_bs_perm <- null_score_bs
      for (ip in 1:length(est_nulls_bs)) {
        permute_decoy <- sample(1:(d + 1), (d + 1), replace = FALSE)
        null_score_bs_perm[, ip] <- null_score_bs[permute_decoy, ip]
      }

      scores_bs[, est_nulls_bs] <- null_score_bs_perm
    }

    rank_scores_bs <- apply(scores_bs, 2, rank, ties.method = "random")
    bs_scores_pvals <- (d + 2 - rank_scores_bs[1, ])/(d + 1)


    # Mirror c Selection MDT Test
    if ("mirror" %in% test_methods) {
      loc <- which(test_methods == "mirror")

      # c = 0.5 lambda = 0.5

      for (iq in 1:nqs) {
        res_c_mirror <- mirandom(scores_bs, 0.5, 0.5, alpha_range[iq], rank_scores = rank_scores_bs)

        S <- res_c_mirror$Discoveries_ind
        FDR_est[ibs, iq, loc] <- sum(estTrueNullInds[resample_indicies][S])/max(length(S), 1)
        if (any(S == 0)) {
          n_target_discs[ibs, iq, loc] <- 0
        } else {
          n_target_discs[ibs, iq, loc] <- length(S)
        }
      }
    }

    if ("lf" %in% test_methods) {

      loc <- which(test_methods == "lf")

      for (iq in 1:nqs) {
        c_FL <- suppressWarnings(max(lat_vals[lat_vals - 1e-12 < alpha_range[iq]]))
        if (c_FL < min(lat_vals)) {
          c_FL <- min(lat_vals)
        }

        res_c_lf <- mirandom(scores_bs, min(c_FL, 0.5), 0.5, alpha_range[iq], rank_scores = rank_scores_bs)

        S <- res_c_lf$Discoveries_ind
        FDR_est[ibs, iq, loc] <- sum(estTrueNullInds[resample_indicies][S])/max(length(S), 1)
        if (any(S == 0)) {
          n_target_discs[ibs, iq, loc] <- 0
        } else {
          n_target_discs[ibs, iq, loc] <- length(S)
        }
      }

    }

    if ("fds" %in% test_methods) {

      # Choose c as the maximum value of the emprical p-values of the observed scores rejected by Storey's FDR controlling procedure Use lambda =
      # estimated lambda Only search in (0, lambda] for c

      loc <- which(test_methods == "fds")

      lam <- lam_est(bs_scores_pvals, d, 0.1)
      pi_0_est1 <- pi_0_est_St(bs_scores_pvals, lam, add = 1)

      for (iq in 1:nqs) {
        c_pvals <- fds_c_select(bs_scores_pvals, d, pi_0_est1, lam, add = 0, alpha_range[iq])

        res_c_fds <- mirandom(scores_bs, min(c_pvals, lam), lam, alpha_range[iq], rank_scores = rank_scores_bs)

        S <- res_c_fds$Discoveries_ind
        FDR_est[ibs, iq, loc] <- sum(estTrueNullInds[resample_indicies][S])/max(length(S), 1)
        if (any(S == 0)) {
          n_target_discs[ibs, iq, loc] <- 0
        } else {
          n_target_discs[ibs, iq, loc] <- length(S)
        }
      }
    }

    if ("fds1" %in% test_methods) {

      # Choose c as the maximum value of the emprical p-values of the observed scores rejected by Storey's FDR controlling procedure + 1/(d + 1)
      # Use lambda = estimated lambda If c > lambda. Set lambda = max(c, lambda)

      loc <- which(test_methods == "fds1")

      lam <- lam_est(bs_scores_pvals, d, 0.1)
      pi_0_est <- pi_0_est_St(bs_scores_pvals, lam, add = 0)

      for (iq in 1:nqs) {
        c_pvals <- fds_c_select(bs_scores_pvals, d, pi_0_est, 0.95, add = 1, alpha_range[iq])

        res_c_fds1 <- mirandom(scores_bs, c_pvals, max(c_pvals, lam), alpha_range[iq], rank_scores = rank_scores_bs)

        S <- res_c_fds1$Discoveries_ind
        FDR_est[ibs, iq, loc] <- sum(estTrueNullInds[resample_indicies][S])/max(length(S), 1)
        if (any(S == 0)) {
          n_target_discs[ibs, iq, loc] <- 0
        } else {
          n_target_discs[ibs, iq, loc] <- length(S)
        }
      }
    }

  }


  # For each bootstrap run and alpha we pick the best method with respect to the target discoveries and record the FDP. Then consider the mean
  # and sd of this 'maximum' estimated FDP for that alpha.

  FDR <- rep(0, nqs)
  SD <- FDR
  for (iq in 1:nqs) {
    best_c_bs_pos <- apply(n_target_discs[, iq, ], 1, which.max)
    FDR_bs <- rep(0, bs_runs)
    for (ibs in 1:bs_runs) {
      FDR_bs[ibs] <- FDR_est[ibs, iq, best_c_bs_pos[ibs]]
    }
    SD[iq] <- sd(FDR_bs)
    FDR[iq] <- mean(FDR_bs)
  }


  m_choice <- rep(0, nqs)
  Fallback <- m_choice

  # For each alpha if FDR_est < FDR_thresh select the best method Use FDR_thres = alpha + pi_1 *4 * SE
  FDR_thres <- alpha_range + (1 - pi_0) * 4 * SD/sqrt(bs_runs)
  # Take method with best average rank over all the bootstrap runs If two methods are equally the best in average rank. Split using the order
  # of test_methods.
  for (iq in 1:nqs) {
    r_target_discs <- t(apply(n_target_discs[, iq, ], 1, rank, ties.method = "last"))
    r_target_discs_bsavg <- apply(r_target_discs, 2, mean)
    if (FDR[iq] > FDR_thres[iq]) {
      Fallback[iq] <- TRUE
      m_choice[iq] <- fallback
    } else {
      Fallback[iq] <- FALSE
      m_choice[iq] <- test_methods[which.max(r_target_discs_bsavg)]
    }
  }

  return(list(m_choice = m_choice, Fallback = Fallback))
}
