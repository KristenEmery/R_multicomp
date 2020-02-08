#' @title
#' FDR control using multiple competition.
#'
#' @description
#' Performs multiple competition with a preselected value of c and \eqn{\lambda}. Uses mirandom mapping scheme to determine W and L.
#'
#' @param scores A \eqn{(d+1) * m} matrix of combined observed and competing scores. Each column corresponds to a indivdual hypothesis while each row corresponds to a observed or competing score set. The first row is the observed scores while the subsequent rows correspond to specific competing scores.
#' @param c Input value of the parameter c. Determines how we select a observed win. Lower value results in a stricter condition.
#' @param lambda Input value of the parameter \eqn{\lambda}. Determines how we select a competing win. Higher value results in a stricter condition.
#' @param alpha The desired FDR significance level.
#' @param rank_scores A \eqn{(d+1) * m} matrix of ranked observed and competing scores. Same structure as the \code{scores} parameter but with the score ranked across observed-competing score sets. Default is NULL which will result in the function calculating the ranks, but can be manually input.
#'
#' @details \code{rank_scores} can be pre-calculated and input so we don't have perform m iterations of ranking every time we consider a new \eqn{\alpha}
#'
#' @return
#'
#' Returns a list:
#'
#' \item{Discoveries}{The W values of the selected discoveries.}
#' \item{Discoveries_ind}{The indicies (location) of the selected discoveries. A response of 0 corresponds to no discoveries.}
#'
#' @author Kristen Emery & Uri Keich
#'
#' @references
#' Emery K, Keich U, Hasam S and Nobel W. (2019) Multiple competition based FDR control. arXiv:1907.01458\cr
#' \url{https://arxiv.org/abs/1907.01458}
#'
#' @export

mirandom <- function(scores, c, lambda, alpha, rank_scores = NULL){

  # Computes MDT with d decoys with both rand and mirror selection for either calibrated or uncalibrated scores using the new parameter c. Uses
  # only a single alpha to incorperate a c that varies over alpha.

  if (is.null(rank_scores)) {
    # Rank each target-decoy set.
    rank_scores <- apply(scores, 2, rank, ties.method = "random")
  }

  n <- ncol(scores)
  d <- nrow(scores) - 1

  lat_vals <- 1:(d)/(d + 1)

  # Take c and lambda of the form i/(d+1)
  lambda <- min(lat_vals[lat_vals >= lambda])
  c <- max(lat_vals[lat_vals <= c])

  if (c + 1e-10 < min(lat_vals)) {
    warning("c is equivalent to 0, you will make no discoveries")
  } else {
    c <- max(lat_vals[lat_vals - 1e-10 <= c])
  }

  if (lambda - 1e-10 > max(lat_vals)) {
    warning("lambda is equivalent to 1")
  } else {
    lambda <- min(lat_vals[lat_vals + 1e-10 >= lambda])
  }

  if (c > lambda) {
    stop("c > lambda, please set c <= lambda")
  }

  ori_select <- rep(FALSE, n)
  decoy_select <- rep(FALSE, n)

  obs_scores_pvals <- (d + 2 - rank_scores[1, ])/(d + 1)

  ori_select[which(obs_scores_pvals - 1e-12 <= c)] <- TRUE
  decoy_select[which(obs_scores_pvals - 1e-12 > lambda)] <- TRUE

  W <- rep(0, n)

  W[ori_select] <- scores[1, ori_select]
  decoy_select_loc <- which(decoy_select)

  if(sum(decoy_select) > 0){
    if(d > 1){
      n_decoys_in_draw <- floor(c*(d + 1) + 1e-12)   #the number of decoys the losing observations will be mapped to
      n_obs_ranks <- d + 1 - ceiling((d + 1) * lambda - 1e-12)  #the number of ranks which result in a decoy selection
      max_mapped_decoy_rank <- matrix(0L, n_obs_ranks,1)   #the highest decoy rank each losing observation is mapped to
      max_mapped_decoy_prob <- max_mapped_decoy_rank  #the probability of mapping to that maximal decoy rank
      current_decoy_rank <- d     #the decoy rank currently being mapped into, starting from the highest rank
      current_decoy_coverage <- 0   #how much of that decoy rank did we cover so far
      uni_decoy_coverage <- n_obs_ranks / n_decoys_in_draw  #how much coverage should each mapped-to-decoy get

      for (i in 1 : n_obs_ranks){
        max_mapped_decoy_rank[i] <- current_decoy_rank;  #assign the current highest decoy rank available
        if (current_decoy_coverage + 1 > uni_decoy_coverage){  #mapping the current observation (i) to this decoy overfills its quota
          max_mapped_decoy_prob[i] <- uni_decoy_coverage - current_decoy_coverage;  #so the probability is (quota) - (coverage so far)
          remainder_obs_coverage_prob <- 1 - max_mapped_decoy_prob[i];  #portion of the current observation we still need to map
          current_decoy_coverage <- uni_decoy_coverage;  #current decoy is saturated
        } else {                                        #current observation can be fully mapped to current decoy
          max_mapped_decoy_prob[i] <- 1;                 #so the probability is 1
          current_decoy_coverage <- 1 + current_decoy_coverage;  #Update current decoy coverage
          remainder_obs_coverage_prob = 0;              #there is no "change", the current observation was fully mapped
        }
        if (current_decoy_coverage >= uni_decoy_coverage - 1e-10){ #is current decoy saturated (up to a roundoff error)?
          # floor(...) below is the number of additional decoy ranks the current observation can saturate
          current_decoy_rank <- current_decoy_rank - floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) - 1;	# allow for roundoff errors
          current_decoy_coverage <- remainder_obs_coverage_prob - floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) * uni_decoy_coverage;
          #and the remainder of the observation is allocated to the new decoy
        }
      }
      if (n_decoys_in_draw > 0 && (current_decoy_coverage > 1e-10 || d - current_decoy_rank != n_decoys_in_draw)){  #sanity check - increase 1e-12 for larger m
        stop('Error in mapping s_i - the winning ranks.')
      }

      obs_ranks <- rank_scores[1,decoy_select]
      rands <- runif(sum(decoy_select))
      mapped_obs_ranks <- max_mapped_decoy_rank[obs_ranks] - ceiling((rands - max_mapped_decoy_prob[obs_ranks]) / uni_decoy_coverage - 1e-12) * (rands > max_mapped_decoy_prob[obs_ranks]);

      rank_decoy_scores <- as.matrix(rank_scores[-1,decoy_select_loc])
      rank_obs_scores <- rank_scores[1,decoy_select_loc]
      for(i in 1:sum(decoy_select)){
        rank_decoy_scores[,i][rank_decoy_scores[,i] > rank_obs_scores[i]] <- rank_decoy_scores[,i][rank_decoy_scores[,i] > rank_obs_scores[i]] - 1
      }

      for(i in 1:sum(decoy_select)){
        dec_choice <- which.min(abs(rank_decoy_scores[,i]  - mapped_obs_ranks[i]))
        W[decoy_select_loc[i]] <- scores[dec_choice + 1, decoy_select_loc[i]]
      }

    } else {
      W[decoy_select_loc] <- scores[2,decoy_select_loc]
    }
  }

  counted_hyp <- as.logical(ori_select | decoy_select)
  n_count <- sum(counted_hyp)
  permW <- sample(1:n_count, n_count, replace = FALSE)
  # Randomly break ties in W
  W_order <- order(W[counted_hyp][permW], decreasing = TRUE)
  W_sort <- W[counted_hyp][permW][W_order]

  ori_sort <- ori_select[counted_hyp][permW][W_order]
  decoy_sort <- decoy_select[counted_hyp][permW][W_order]

  # At each point calculate the number of observed discoveries and competing null discoveries if we were to reject every score above that
  # point.
  Obs <- cumsum(ori_sort)
  Null <- cumsum(decoy_sort) + 1

  FDP_est <- pmin(Null/Obs * c/(1 - lambda), 1)

  res_alpha <- which(FDP_est - 1e-12 <= alpha)
  if (length(res_alpha >= 1)) {
    Last_Disc <- max(res_alpha)
    Disc <- which(ori_sort[1:Last_Disc])
    score_ind <- W_order[Disc]
    Discoveries_ind <- c(1:n)[counted_hyp][permW[score_ind]]
    Discoveries <- W[Discoveries_ind]
  } else {
    Discoveries_ind <- 0
    Discoveries <- NULL
  }

  return(list(Discoveries = Discoveries, Discoveries_ind = Discoveries_ind))
}
