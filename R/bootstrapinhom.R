#' # Bootstrap Bias/CIs for inhomogeneous version
#'
#' @param N Number of bootstrap draws.
#' @param n Number of replications in each sample.
#' @param estimate Vector of coefficient estimates obtained by a fitting procedure.
#' @param window Observation window.
#' @param lambda Lambda derived from the estimate.
#' @param lambda_max Upper bound to lambda.
#' @param covariates List of covariates
#' @param ndummy Number of dummy points used in the calls to the logistic regression.
#' @param alpha Probability to use in the bootstrap CIs.
#' @param nthreads Number of CPU threads to use.
#' @param force_nu Force nu to a given value?
#' @param parallel Use faster parallel version?
#'
#' @return List containing bootstrap bias, and CIs.
#'
#' @importFrom parallel clusterExport makeCluster parLapply stopCluster
#' @importFrom spatstat owin
#' @export
bootstrapinhom <- function(N,
                           n,
                           estimate,
                           window = owin(),
                           lambda,
                           lambda_max,
                           covariates,
                           ndummy = 1e4,
                           alpha = 0.05,
                           nthreads = 1,
                           force_nu,
                           parallel = TRUE) {
  if(!missing(force_nu)) {
    nu <- force_nu
  } else if(is.na(match('nu', names(estimate)))) {
    nu <- 1
  } else {
    nu <- estimate[match('nu', names(estimate))]
  }
  if(missing(force_nu)) {
    force_nu <- NA
  }
  cat("Starting sampling of bootstrap configurations.\n")
  samples <- lapply(seq_len(N), function(i) {
    cat(paste0(i, "...\n"))
    rcompppinhom(n = n, lambda = lambda, lambda_max = lambda_max, nu = nu, window = window)
  })
  cat("Done sampling.\n")
  cat(paste0("Mean points in sampled configurations: ", paste0(sapply(samples, function(sample) mean(sapply(sample, function(conf) length(conf$x)))), collapse = ", "), ".\n"))

  cat("Starting fitting procedure.\n")
  if(parallel) {
    libpath <- .libPaths()
    varlist <- c("rcomfitlogit", "samples", "covariates", "force_nu", "libpath", "ndummy")
    cl <- makeCluster(nthreads)
    clusterExport(cl, varlist = varlist, envir = environment())
    coefs <- parLapply(cl = cl, seq_len(N), function(i) {
      .libPaths(libpath)

      if(is.na(force_nu)) {
        rcomfitlogit(samples[[i]], covariates = covariates, ndummy = ndummy)$coef
      } else {
        rcomfitlogit(samples[[i]], covariates = covariates, ndummy = ndummy, force_nu = force_nu)$coef
      }
    })
    stopCluster(cl)
  } else {
    coefs <- lapply(seq_len(N), function(i) {
      cat(paste0(i, "...\n"))
      if(is.na(force_nu)) {
        rcomfitlogit(samples[[i]], covariates = covariates, ndummy = ndummy)$coef
      } else {
        rcomfitlogit(samples[[i]], covariates = covariates, force_nu = force_nu, ndummy = ndummy)$coef
      }
    })
  }
  cat("Done with the fitting.\n")

  parameters <- matrix(NA, nrow = N, ncol = length(estimate))
  colnames(parameters) <- names(estimate)
  for(i in seq_len(nrow(parameters))) {
    for(j in seq_len(ncol(parameters))) {
      stopifnot(length(coefs[[i]]) == ncol(estimate))
      parameters[i, j] <- coefs[[i]][which(names(estimate)[j] == names(coefs[[i]]))]
    }
  }

  bias <- colMeans(parameters) - estimate
  sorted_parameters <- apply(parameters, 2, sort)
  empirical_cdf <- seq_len(N) / N
  lower_parameters <- apply(sorted_parameters, 2, function(a) a[which(empirical_cdf >= alpha / 2)[1]])
  upper_parameters <- apply(sorted_parameters, 2, function(a) a[which(empirical_cdf >= 1 - alpha / 2)[1]])

  list(bias = bias,
       lower_parameters = lower_parameters,
       upper_parameters = upper_parameters)
}
