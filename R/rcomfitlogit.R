#' Fit a COM-Poisson point process
#'
#' @param configurations List of `spatstat::ppp` configurations.
#' @param force_nu For nu to a given value?
#' @param covariates List of covariates.
#' @param ndummy Number of dummy points.
#'
#' @return List with two elements: the actual `spatstat` fit, and the corresponding coefficient estimates.
#'
#' @import rgdal
#' @importFrom methods is
#' @importFrom raster extract
#' @importFrom spatstat.geom quadscheme.logi superimpose
#' @export
rcomfitlogit <- function(configurations, force_nu, covariates = list(), ndummy = 1e4) {
  N <- length(configurations)
  stopifnot(N > 0)

  # Convert configurations to a list
  if(!is.list(configurations)) {
    configurations <- list(configurations)
  }

  # Are the covariates a list of time-independent covariates,
  # or a list of lists of covariates, each one corresponding
  # to a configuration?
  is_time_dependent_covariates <- ifelse(length(covariates) == 0,
                                         FALSE,
                                         is(covariates[[1]], "list") &
                                           length(covariates[[1]]) == length(configurations))

  # Combine the configurations into a single one
  complete_configuration <- Reduce(superimpose, configurations)

  # Create the quadrature scheme and dummy points
  quadrature_scheme <- quadscheme.logi(data = complete_configuration,
                                       nd = c(ndummy, N),
                                       dummytype = "binomial")
  dummy <- quadrature_scheme$dummy

  # This is not equal to ndummy if the window is an `im`
  nd <- as.integer(length(dummy$x) / N)

  # Construct nu vector
  log_n <- lapply(configurations, function(conf) rep(-log(length(conf$x)), length(conf$x)))
  log_1pn <- lapply(seq_len(length(configurations)), function(i) {
    if(i < length(configurations)) {
      rep(-log(1 + length(configurations[[i]]$x)), nd)
    } else {
      rep(-log(1 + length(configurations[[i]]$x)), length(dummy$x) - (N - 1) * as.integer(length(dummy$x) / N) )
    }
  })
  nu <- c(Reduce(c, log_n), Reduce(c, log_1pn))

  # Construct offset vector
  os <- log(N) - nu
  if(!missing(force_nu)) {
    os <- os + force_nu * nu
  }

  # Function to extract from covariate depending on whether or not it's a raster
  extract_from_covariate <- function(covariate, x, y) {
    if(is(covariate, "RasterLayer")) {
      z <- terra::extract(covariate, data.frame(x = x, y = y))
    } else {
      z <- as.function(covariate)(x, y)
    }
    if(length(z) != length(x) | length(x) != length(y)) {
      stop("One of the covariates does not output a vector of the same length as the vector to which it has been applied.")
    }
    z
  }

  # Construct covariates matrix
  covariate_matrix <- if(is_time_dependent_covariates) {
    do.call("cbind", lapply(covariates, function(covariate) {
      presence <- Reduce(c, lapply(seq_len(length(covariate)), function(i) {
        extract_from_covariate(covariate = covariate[[i]],
                               x = configurations[[i]]$x,
                               y = configurations[[i]]$y)
      }))
      absence <- Reduce(c, lapply(seq_len(length(covariate)), function(i) {
        indices <- if(i < length(covariate)) {
          (1 + (i - 1) * nd):(i * nd)
        } else { # If the window is irregular, there might be some leftover quadrature points,
          # put them in the last bin.
          (1 + (i - 1) * nd):length(dummy$x)
        }
        extract_from_covariate(covariate = covariate[[i]],
                               x = dummy$x[indices],
                               y = dummy$y[indices])
      }))
      c(presence, absence)
    }))
  } else {
    do.call("cbind", lapply(covariates, function(covariate) {
      extract_from_covariate(covariate = covariate,
                             x = c(complete_configuration$x, dummy$x),
                             y = c(complete_configuration$y, dummy$y))
    }))
  }

  # Add colnames to covariate matrix
  if(is.null(names(covariates)) & length(covariates) > 0) {
    colnames(covariate_matrix) <- paste0("Covariate_", seq_len(ncol(covariate_matrix)))
  } else {
    colnames(covariate_matrix) <- names(covariates)
  }

  # Fit the COM-Poisson model
  dat <- data.frame(nu = nu, os = os)
  if(length(covariates) > 0) {
    dat <- cbind(dat, covariate_matrix)
  }
  trend <- stats::as.formula(paste0("~ 1 + offset(os)",
                                    ifelse(missing(force_nu), " + nu", ""),
                                    ifelse(length(covariates) > 0, " + ", ""),
                                    paste0(colnames(covariate_matrix), collapse = " + ")))

  fit_spatstat <- ppm(Q = quadrature_scheme,
                      trend = trend,
                      data = dat,
                      correction = "none",
                      method = "logi")
  cf <- coef(fit_spatstat)
  nu <- cf[names(cf) == 'nu']
  if(length(nu) > 0) {
    if(nu < 0) {
      return(rcomfitlogit(configurations = configurations,
                          force_nu = 0,
                          covariates = covariates,
                          ndummy = ndummy))
    }
  } else {
    cf <- c(cf, nu = force_nu)
  }
  list(fit = fit_spatstat, coef = cf)
}
