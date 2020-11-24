#' em
#'
#' Run expectation maximization algorithm to fit parameters to each subject of an experiment
#'
#' @param dat data.frame; data
#' @param mus vector; starting mean parameters
#' @param sigma vector; starting sd in parameter
#' @param likfun function; the objective function
#' @param bic logical; if TRUE, calculate the group-level bic after fitting model
#' @param ndata integer; number of observations. only needed to calculate bic
#' @param emtol numeric; tolerance for convergence
#' @param rel_and_abs logical; use both relative and absolute tolerance to evaluate convergence
#' @param maxiter integer; maximum number of iterations. Will stop if maxiter reached without convergence
#' @param lower numeric vector; lower bounds for parameters
#' @param upper numeric vector; upper bounds for parameters
#' @param perturb_start numeric; coefficient of variation to perturb starting optimization parameters
#' @param parallel logical; if TRUE, will optimize subjects in parallel (within an iteration)
#' @param return_ind_df logical; if TRUE, will only return a data frame of individual results (parameters and likelihood)
#' @param verbose logical; if TRUE, will print progress at each iteration
#' @param ... additional arguments passed to modelfitr::fit_model
#'
#' @return if return_ind_df = FALSE, returns a list containing the following fields:
#' \item{mus}{vector of mean parameters}
#' \item{sigma}{vector of standard deviation of each parameter}
#' \item{x}{npar x nsubjects matrix with individual parameters}
#' \item{l}{vector of likelihood for each subject}
#' \item{h}{npar x npar x nsubjects array. each slice contains the hessian matrix for each subject at that subjects individual parameters}
#'
#' if return_ind_df = TRUE, returns a data frame with the following columns:
#' \item{Subject}{Subject ID}
#' \item{Parameters}{a separate column for each parameter}
#' \item{lik}{The likelihood for the individual subject}
#' \item{bic}{The group-level bic for the model (only if)}
#'
#' @export
em <- function(dat, mus, sigma, likfun, bic = F, ndata = NULL, emtol = 1e-5, rel_and_abs = TRUE, maxiter = 100L, lower = -Inf, upper = Inf, perturb_start = 0, parallel = F, return_ind_df = F, verbose = T, ...) {
  nparam <- length(mus)
  pnames <- names(mus)
  subs <- unique(dat$Subject)

  ### check the length of mus and sigmas
  if (length(mus) != length(sigma)) stop("mus and sigma must be the same length!")

  ### set up design matrix and sigma matrix ###
  dm <- designmatrix(nparam, length(subs))
  sigma <- Matrix::Diagonal(length(sigma), sigma)

  ### first iteration ###
  oldparams <- c(mus, Matrix::diag(sigma))
  fit <- estep(dat, subs, mus, sigma, likfun, lower = lower, upper = upper, perturb_start = perturb_start, parallel = parallel, verbose = verbose, ...)
  par_list <- mstep(fit$x, dm, fit$h, sigma)
  mus <- par_list$mus
  sigma <- par_list$sigma
  newparams <- c(mus, Matrix::diag(sigma))
  pardiff <- calc_par_diff(oldparams, newparams, rel_and_abs)

  ### repeat until convergence ###
  iter <- 1
  while (max(pardiff) > emtol & iter < maxiter) {
    iter <- iter + 1
    oldparams <- newparams

    fit <- estep(dat, subs, mus, sigma, likfun, startx = fit$x, lower = lower, upper = upper, perturb_start = perturb_start, parallel = parallel, verbose = verbose, ...)
    par_list <- mstep(fit$x, dm, fit$h, sigma)
    mus <- par_list$mus
    sigma <- par_list$sigma
    newparams <- c(mus, Matrix::diag(sigma))
    pardiff <- calc_par_diff(oldparams, newparams, rel_and_abs)

    if (verbose) {
      cat("\n")
      cat("iter:", iter, "\n")
      cat("mus:", round(mus, 5), "\n")
      cat("sigma:", round(Matrix::diag(sigma), 5), "\n")
      cat("change:", round(pardiff, 5), "\n")
      cat("max:", round(max(pardiff), 5), "\n")
      cat("\n")
    }
  }

  dimnames(sigma) <- list(pnames, pnames)
  dimnames(fit$x) <- list(pnames, subs)
  names(fit$l) <- subs
  dimnames(fit$h) <- list(pnames, pnames, subs)
  res <- list(mus = mus, sigma = sigma, x = fit$x, l = fit$l, h = fit$h)

  if (bic) {
    if (is.null(ndata)) {
      warning("ndata needed to calculate bic!")
      res$bic <- NA
    } else {
      res$bic <- ibic(res, ndata)
    }
  }

  if (return_ind_df) {
    return(get_ind_df(res))
  } else {
    return(res)
  }

}

calc_par_diff <- function(oldparams, newparams, rel_and_abs = TRUE) {
  absdiff <- abs(newparams - oldparams)
  reldiff <- abs((newparams - oldparams) / oldparams)
  if (rel_and_abs) {
    absdiff <- abs(newparams - oldparams)
  } else {
    absdiff = rep(Inf, length(reldiff))
  }
  mapply(function(x, y) min(x, y), x=reldiff, y=absdiff)
}

estep <- function(dat, subs, mus, sigma, likfun, startx = NULL, lower = -Inf, upper = Inf, perturb_start = 0, parallel = F, verbose = T, ...) {
  nsub <- length(subs)
  nparam <- length(mus)
  if (is.null(startx)) startx <- matrix(rep(mus, nsub), ncol = nsub)
  startx = matrix(rnorm(length(startx), startx, abs(startx * perturb_start)), ncol=nsub)

  h <- array(0, c(nparam, nparam, nsub))
  l <- numeric(nsub)
  x <- matrix(0, nrow = nparam, ncol = nsub)

  if (!parallel) {
    for (i in 1:nsub) {
      if (verbose) cat(i, "..", sep = "")
      sub_fit <- opt_sub(gaussian_prior, startx[, i], lower = lower, upper = upper, mus = mus, var = sigma, dat = dat[dat$Subject == subs[i], ], likfun = likfun, ...)
      x[, i] <- sub_fit$par
      l[i] <- sub_fit$lik
      h[, , i] <- sub_fit$hess
    }
    if (verbose) cat("\n")
  } else {
    dat_list <- split(dat, dat$Subject)
    startx_list <- simplify2array(apply(startx, 2, list))
    opt_sub_parallel <- function(d_sub, stx, ...) {
      opt_sub(gaussian_prior, stx, lower = lower, upper = upper, dat = d_sub, ...)
    }
    fit_list <- parallel::mcmapply(opt_sub_parallel,
      d_sub = dat_list, stx = startx_list,
      MoreArgs = c(list(mus = mus), var = sigma, likfun = likfun, list(...)),
      mc.cores = parallel::detectCores() - 1
    )

    if (length(fit_list) == 0) browser()
    x <- simplify2array(fit_list[1, ])
    l <- as.numeric(fit_list[2, ])
    h <- simplify2array(fit_list[3, ])
    if (class(x) == "numeric") {
      x <- matrix(x, nrow=1)
      h <- array(h, dim=c(1, 1, length(x)))
    }
  }

  list(x = x, l = l, h = h)
}

mstep <- function(x, dm, h, sigma, full = F) {
  nsub <- ncol(x)
  isigma <- solve(sigma)

  ### find betas ###
  m1 <- apply(dm, 3, function(x) t(x) %*% isigma %*% x)
  if (class(isigma) != "matrix") {
    mean_m1 <- numeric(nrow(m1[[1]]))
    for (i in 1:nsub) {
      mean_m1 <- mean_m1 + diag(m1[[i]])
    }
    mean_m1 <- mean_m1 / nsub
    mean_m1 <- Matrix::Diagonal(length(mean_m1), mean_m1)
  } else {
    m1 <- array(m1, c(nrow(sigma), nrow(sigma), nsub))
    mean_m1 <- apply(m1, 1:2, mean)
  }

  mean_m2 <- numeric(nrow(x))
  for (i in 1:nsub) {
    mean_m2 <- mean_m2 + t(dm[, , i]) %*% isigma %*% x[, i]
  }
  mean_m2 <- mean_m2 / nsub

  mus <- as.numeric(solve(mean_m1) %*% mean_m2)

  ### find sigma ###

  sigma <- (x - mus) %*% t(x - mus) / nsub + apply(h, 1:2, mean, na.rm = T)

  if (!full) {
    sigma <- Matrix::Diagonal(nrow(sigma), Matrix::diag(sigma))
  }

  list(mus = mus, sigma = sigma)
}

opt_sub <- function(obj, start, lower = -Inf, upper = Inf, package = "optimx", method = "BFGS", opt_args = list(), obj_args = list(), ...) {
  it <- 1
  fit <- modelfitr::fit_model(obj,
                              start,
                              lower = lower,
                              upper = upper,
                              package = package,
                              method = method,
                              hessian = TRUE,
                              opt_args = opt_args,
                              obj_args = obj_args,
                              ...)

  return(list(par = fit$pars, lik = fit$value, hess = solve(fit$hess)))
}

gaussian_prior <- function(param_values, mus, var, dat, likfun, ...) {
  d <- length(param_values)
  lp <- -d / 2 * log(2 * pi) - 1 / 2 * log(Matrix::det(var)) - 1 / 2 * t(param_values - mus) %*% solve(var) %*% (param_values - mus)
  nll <- likfun(param_values, dat, ...)
  map <- nll - lp[1]
  if (is.infinite(map) | is.na(map)) map <- 1e10
  if (is.nan(map)) browser()
  map
}

designmatrix <- function(npar, nsub) array(diag(npar), c(npar, npar, nsub))

# ibic <- function(x, l, h, mus, sigma, ndata) {
#   nHypPar <- length(c(mus, Matrix::diag(sigma)))
#   lml(x, l, h) + nHypPar / 2 * log(ndata)
# }

# lml <- function(x, l, h) {
#   # laplace approximation to the log marginal likelihood
#   nparam <- nrow(x)
#   nsub <- ncol(x)
#
#   log_det_hess <- apply(h, 3, function(x) log(det(x)))
#   -nparam / 2 * log(2 * pi) * nsub + sum(l) - sum(log_det_hess) / 2
# }

fit_all <- function(dat, subs, obj, start, lower = -Inf, upper = Inf, parallel = F, ...) {
  nparam <- length(start)
  nsub <- length(subs)
  h <- array(0, c(nparam, nparam, nsub))
  l <- numeric(nsub)
  x <- matrix(0, nrow = nparam, ncol = nsub)

  if (!parallel) {
    for (i in 1:nsub) {
      cat(i, "..", sep = "")
      sub_fit <- opt_sub(obj, start, lower = lower, upper = upper, dat = dat[dat$Subject == subs[i], ], ...)
      x[, i] <- sub_fit$par
      l[i] <- sub_fit$lik
      h[, , i] <- sub_fit$hess
    }
  } else {
    dat_list <- split(dat, dat$Subject)
    opt_sub_parallel <- function(d_sub, stx, ...) {
      opt_sub(gaussian_prior, stx, lower = lower, upper = upper, dat = d_sub, ...)
    }
    fit_list <- parallel::mclapply(dat_list, opt_sub_parallel,
      MoreArgs = c(stx = start, list(...)),
      mc.cores = parallel::detectCores()
    )

    x <- simplify2array(fit_list[1, ])
    l <- as.numeric(fit_list[2, ])
    h <- simplify2array(fit_list[3, ])
  }

  list(x = x, l = l, h = h)
}
