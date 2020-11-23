#' ibic
#'
#' Get integrated bayes information criterion from em fit
#'
#' @param em_fit list; output from GaussExMax::em
#'
#' @return numeric; integrated bayes information criterion
#'
#' @export
ibic <- function(em_fit, ndata) {
  nHypPar <- length(c(em_fit$mus, Matrix::diag(em_fit$sigma)))
  lml(em_fit) + nHypPar / 2 * log(ndata)
}
