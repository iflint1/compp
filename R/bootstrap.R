#' # Bootstrap Bias/CIs for homogeneous version
#'
#' @param N Number of bootstrap draws.
#' @param n Number of replications in each sample.
#' @param estimate Vector of coefficient estimates obtained by a fitting procedure.
#' @param window Observation window.
#' @param ndummy Number of dummy points used in the calls to the logistic regression.
#' @param force_nu Force nu to a given value?
#' @param nthreads Number of CPU threads to use.
#' @param alpha Probability to use in the bootstrap CIs.
#'
#' @return List containing bootstrap bias, and CIs.
#'
#' @importFrom parallel clusterExport makeCluster parLapply stopCluster
#' @importFrom spatstat owin
#' @export
bootstrap <- function(N,
                      n,
                      estimate,
                      window = owin(),
                      ndummy = 1e4,
                      force_nu,
                      nthreads = 1,
                      alpha = 0.05) {
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

  lambda  <- exp(estimate[match("(Intercept)", names(estimate))])

  # Parallel section
  varlist <- c("n", "libpath")
  libpath <- .libPaths()
  cl <- makeCluster(nthreads)
  clusterExport(cl, varlist = varlist, envir = environment())
  coefs <- parLapply(cl = cl, seq_len(N), function(i) {
    .libPaths(libpath)

    require(COMPoissonReg)

    samples <- rcomppp(n = n,
                       lambda = lambda,
                       nu = nu,
                       window = window)
    if(is.na(force_nu)) {
      rcomfitlogit(samples, ndummy = ndummy)$coef
    } else {
      rcomfitlogit(samples, ndummy = ndummy, force_nu = force_nu)$coef
    }

  })
  stopCluster(cl)

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
  empirical_cdf <- 1:N / N
  lower_parameters <- apply(sorted_parameters, 2, function(a) a[which(empirical_cdf >= alpha / 2)[1]])
  upper_parameters <- apply(sorted_parameters, 2, function(a) a[which(empirical_cdf >= 1 - alpha / 2)[1]])

  list(bias = bias,
       lower_parameters = lower_parameters,
       upper_parameters = upper_parameters)
}
