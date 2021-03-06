#' Multiple competition based FDR control.
#'
#'
#' This package implements multiple competition procedures to control the false discovery rate (FDR) in multiple hypotheses testing.
#' This procedure forgoes the use of p-values and instead uses direct competition between the observed scores and a set of pre-generated null scores for each hypothesis.
#'
#'
#' @section Outline:
#' Multiple competition allows for p-value free multiple hypothesis testing. By (externally) randomly generating sets of
#' null scores for each hypothesis and directly comparing these to the original observed scores we are able to estimate and
#' control the false discovery rate (FDR) for a given list of discoveries. These discoveries are generally those hypotheses that have original scores that are
#' "better" than their corresponding competing null scores.
#'
#' The procedure is an extension of the competition framework present in target-decoy competition and knockoff+, designed to utilize
#' multiple competing scores instead of just one. This allows for significant power gain without sacrificing FDR control.
#'
#'
#' @section Functions:
#' The functions can be separated into two levels; a primary level that is suitable for most users who are simply interested in application
#' and a secondary level for more advanced users looking for more control.
#'
#' Primary level:
#' \itemize{
#' \item{\code{\link{multidecoy_comp}}: The main function of this package. Takes a set of observed scores and
#' competing scores and returns a list of discoveries. Most users will use this function only.}
#' }
#'
#' Secondary level:
#'
#' \itemize{
#' \item{\code{\link{mirandom}}: Once the testing method has been selected (see accompanying paper for
#' more details) performs multiple decoy competition.}
#' \item{\code{\link{fds_c_select}}: Calculates \eqn{c}, how to determine an original score win, according the
#' Finite decoy Storey formulation.}
#' \item{\code{\link{lam_est}}: Calculates \eqn{\lambda}, how to determine a competing score win. }
#' \item{\code{\link{LBM_method_select}}: Selects the method of choosing \eqn{c} and \eqn{\lambda} using a bootstrap approach.}
#' \item{\code{\link{label_heuristic}}: Estimates null and alternative labels for the data for construction of
#' the bootstrap samples in LBM.}
#' \item{\code{\link{pi_0_est_St}}: Estimates \eqn{\pi_0} according to Storey, Taylor and Tibshirani's (2004) formulation.}
#'
#' }
#'
#'
#' @references
#' Emery K, Keich U, Hasam S and Nobel W. (2019) Multiple competition based FDR control. arXiv:1907.01458\cr
#' \url{https://arxiv.org/abs/1907.01458}
#'
#' @docType package
#' @name multicomp
#'
#'
NULL
