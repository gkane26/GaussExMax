#' get_ind_df
#'
#' Get data frame of individual parameters from EM fit
#'
#' @param em_fit list; output from GaussExMax::em
#'
#' @return returns a data frame with the following columns:
#' \item{Subject}{Subject ID}
#' \item{Parameters}{a separate column for each parameter}
#' \item{lik}{The likelihood for the individual subject}
#' \item{bic}{The group-level bic for the model (only if)}
#'
#' @export
get_ind_df <- function(em_fit) {
  res_df <- data.frame(Subject=colnames(em_fit$x),
                       t(em_fit$x),
                       lik = em_fit$l,
                       row.names = NULL)
  res_df$bic = em_fit$bic
  res_df
}
