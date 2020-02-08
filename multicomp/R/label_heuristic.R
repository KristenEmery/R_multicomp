#' @title
#' Estimate null and alternative labels for a set of hypotheses.
#'
#' @description
#' Estimates null and alternative labels for use in segmented resampling.
#' Calculates a W and L score for each hypothesis according to the mirandom mapping scheme with c pre-selected and \eqn{\lambda = c}.
#' Labels are determined by estimating the number of false discoveries in each segment and randomly assigning labels to hypotheses in that segment weighted by the empirical p-values.
#'
#' @param scores A \eqn{(d+1) * m} matrix of combined observed and competing scores. Each column corresponds to a indivdual hypothesis while each row corresponds to an observed or competing score set. The first row is the observed scores while the subsequent rows correspond to specific competing scores.
#' @param c Input value of the parameter c. Determines how we select a observed win. Lower value results in a stricter condition.
#' @param rank_scores A \eqn{(d+1) * m} matrix of ranked observed and competing scores. Same structure as the \code{scores} parameter but with the score ranked across observed-competing sets. Default is NULL which will result in the function caluclating the ranks, but can be manually input.
#'
#' @return
#' \item{estTrueNullInds}{A vector of labels for each hypothesis indicating whether a specific hypothesis was estimated as truly null. 1 = Null, 0 = Alternative.}
#'
#' @author Kristen Emery & Uri Keich
#'
#' @references
#' Emery K, Keich U, Hasam S and Nobel W. (2019) Multiple competition based FDR control. arXiv:1907.01458\cr
#' \url{https://arxiv.org/abs/1907.01458}
#'
#' @export

label_heuristic <- function(scores, c, rank_scores = NULL) {
  # Heuristic that assigns to each score a Null label (TRUE) or alternative label (FALSE)

  if (is.null(rank_scores)) {
    rank_scores <- apply(scores, 2, rank, ties.method = "random")  #Rank each target-decoy set.
  }

  n <- ncol(scores)
  d <- nrow(scores) - 1

  lat_vals <- 1:(d)/(d + 1)

  c <- max(lat_vals[lat_vals <= c])

  obs_scores_pvals <- (d + 2 - rank_scores[1, ])/(d + 1)
  rev_pvals <- 1 - obs_scores_pvals

  ori_select <- rep(FALSE, n)
  decoy_select <- rep(FALSE, n)

  ori_select[which(obs_scores_pvals - 1e-12 <= c)] <- TRUE
  decoy_select[which(obs_scores_pvals - 1e-12 > c)] <- TRUE

  W <- rep(0, n)
  L <- rep(0, n)

  W[ori_select] <- scores[1, ori_select]  #Take original score
  decoy_select_loc <- which(decoy_select)
  L[decoy_select] <- 1

  if (sum(decoy_select > 0)) {
    if (d > 1) {
      n_decoys_in_draw <- floor(c * (d + 1))
      n_obs_ranks <- d + 1 - n_decoys_in_draw
      max_mapped_decoy_rank <- matrix(0L, n_obs_ranks, 1)
      max_mapped_decoy_prob <- max_mapped_decoy_rank
      current_decoy_rank <- d
      current_decoy_coverage <- 0
      uni_decoy_coverage <- n_obs_ranks/n_decoys_in_draw

      for (i in 1:n_obs_ranks) {
        max_mapped_decoy_rank[i] <- current_decoy_rank
        if (current_decoy_coverage + 1 > uni_decoy_coverage) {
          max_mapped_decoy_prob[i] <- uni_decoy_coverage - current_decoy_coverage
          remainder_obs_coverage_prob <- 1 - max_mapped_decoy_prob[i]
          current_decoy_coverage <- uni_decoy_coverage
        } else {
          max_mapped_decoy_prob[i] <- 1
          current_decoy_coverage <- 1 + current_decoy_coverage
          remainder_obs_coverage_prob <- 0
        }
        if (current_decoy_coverage >= uni_decoy_coverage - 1e-10) {
          current_decoy_rank <- current_decoy_rank - floor(remainder_obs_coverage_prob/uni_decoy_coverage + 1e-12) - 1
          current_decoy_coverage <- remainder_obs_coverage_prob - floor(remainder_obs_coverage_prob/uni_decoy_coverage + 1e-12) * uni_decoy_coverage
        }
      }
      if (n_decoys_in_draw > 0 && (current_decoy_coverage > 1e-10 || current_decoy_rank != n_obs_ranks - 1)) {
        stop("Error in mapping s_i - the winning ranks.")
      }
      # next, apply the mapping to every observation that lost
      obs_ranks <- rank_scores[1, decoy_select]
      rands <- runif(sum(decoy_select))
      mapped_obs_ranks <- max_mapped_decoy_rank[obs_ranks] - ceiling((rands - max_mapped_decoy_prob[obs_ranks])/uni_decoy_coverage - 1e-12) *
        (rands > max_mapped_decoy_prob[obs_ranks])

      rank_decoy_scores <- as.matrix(rank_scores[-1, decoy_select_loc])
      rank_obs_scores <- rank_scores[1, decoy_select_loc]
      for (i in 1:sum(decoy_select)) {
        rank_decoy_scores[, i][rank_decoy_scores[, i] > rank_obs_scores[i]] <- rank_decoy_scores[, i][rank_decoy_scores[, i] > rank_obs_scores[i]] -
          1
      }

      # Rank only the decoy location

      for (i in 1:sum(decoy_select)) {
        dec_choice <- which.min(abs(rank_decoy_scores[, i] - mapped_obs_ranks[i]))
        W[decoy_select_loc[i]] <- scores[dec_choice + 1, decoy_select_loc[i]]
      }

    } else {
      W[decoy_select_loc] <- scores[2, decoy_select_loc]
    }
  }

  W_order <- order(W, decreasing = TRUE)  #Randomly break ties in W
  W_sort <- W[W_order]
  L_sort <- L[W_order]  #Order scores and corresponding labels in descending order
  pvals_sort <- rev_pvals[W_order]
  # Split the scores up into segments
  estTrueNullInds <- rep(TRUE, n)

  prev_seg_end <- 0
  seg_end <- 1
  prev_rej <- 0
  c_false_null <- NULL
  min_seg_length <- 1

  # Calculate the cumulative FDR in each segment and use that to estimate the number of null hypotheses in the segment. Then randomly select
  # that many hypotheses from the segment and assign them as nulls.

  while (seg_end <= n) {
    seg_select <- (prev_seg_end + 1):(seg_end)
    prev_seg_end <- seg_end

    c_seg_select <- 1:seg_end
    L_seg <- L_sort[seg_select]
    L_c_seg <- L_sort[c_seg_select]
    n_false_null <- max((length(L_c_seg) - sum(L_c_seg)) - round(sum(L_c_seg) * c/(1 - c)) - sum(c_false_null), 0)
    c_false_null <- c(c_false_null, n_false_null)

    if (n_false_null > 0) {
      extended_seg <- c((prev_rej + 1):seg_end)
      est_fn_inds <- sample(1:length(extended_seg), n_false_null, replace = FALSE, prob = pvals_sort[extended_seg])
      est_fn <- extended_seg[est_fn_inds]
      prev_rej <- max(est_fn)
      estTrueNullInds[W_order][est_fn] <- FALSE
    }

    if (n_false_null == 0) {
      min_seg_length <- min_seg_length + 1
    }

    if (seg_end == n) {
      seg_end <- n + 1
    } else {
      seg_end <- min(seg_end + min_seg_length, n)
    }
  }

  return(estTrueNullInds)
}
