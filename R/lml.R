#' lml
#'
#' Get log marginal likelihood from an em fit
#'
#' @param em_fit list; output from GaussExMax::em
#'
#' @return numeric; the log marginal likelihood across subjects
#'
#' @export
lml <- function(em_fit) {
  # laplace approximation to the log marginal likelihood
  nparam <- nrow(em_fit$x)
  nsub <- ncol(em_fit$x)

  log_det_hess <- apply(em_fit$h, 3, function(x) log(det(x)))
  -nparam / 2 * log(2 * pi) * nsub + sum(em_fit$l) - sum(log_det_hess) / 2
}
