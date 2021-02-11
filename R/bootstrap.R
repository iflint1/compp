#' # Bootstrap Bias/CIs for homogeneous version
#'
#' @param N Number of bootstrap draws.
#' @param n Number of replications in each sample.
#' @param estimate Vector of coefficient estimates obtained by a fitting procedure.
#' @param window Observation window.
#' @param ndummy Number of dummy points used in the calls to the logistic regression.
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
                      nthreads = 1,
                      alpha = 0.05) {
  lambda  <- exp(estimate[match("(Intercept)", names(estimate))])
  samples <- rcomppp(n = n * N,
                     lambda = lambda,
                     nu = estimate[match("nu", names(estimate))],
                     window = window)

  # Parallel section
  cl <- makeCluster(nthreads)
  clusterExport(cl, varlist = c("rcomfitlogit", "samples"))
  coefs <- parLapply(cl = cl, seq_len(N), function(n) {
    rcomfitlogit(samples[(1 + (i - 1) * n):(i * n)], ndummy = ndummy)$coef
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
