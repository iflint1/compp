draw_points <- function(total_points,
                        draw_points_lambda,
                        lambda_max) {
  points <- matrix(NA, nrow = total_points, ncol = 2)

  lambda_draws <- draw_points_lambda()
  index <- 1

  simultaneous_draws <- nrow(lambda_draws)

  if(total_points > 0) {
    for(i in seq_len(total_points)) {
      while(lambda_draws[index, 3] < stats::runif(1) * lambda_max) {
        if(index >= simultaneous_draws) {
          lambda_draws <- draw_points_lambda()
          index <- 1
        } else {
          index <- index + 1
        }
      }
      points[i, 1] <- lambda_draws[index, 1]
      points[i, 2] <- lambda_draws[index, 2]
      if(index >= simultaneous_draws) {
        lambda_draws <- draw_points_lambda()
        index <- 1
      } else {
        index <- index + 1
      }
    }
  }
  points
}

#' Inhomogeneous COM-Poisson point process simulation
#'
#' @param n Number of replications.
#' @param lambda Quasi-intensity.
#' @param lambda_max Upper bound to quasi-intensity.
#' @param nu Dispersion parameter
#' @param window Observation window.
#' @param simultaneous_draws Number of simultaneous draws to do when running the rejection sampling.
#' @param nquad Number of quadrature points to use to approximate the integral of the quasi-intensity.
#' @param drop Logical. If n = 1 and drop = TRUE (the default), the result will be a point pattern,
#' rather than a list containing a single point pattern.
#' @param use_cpp Use C++ version of the code?
#'
#' @return List of `spatstat::ppp` samples of an inhomogeneous COM-Poisson point process.
#'
#' @import COMPoissonReg
#' @importFrom spatstat.geom area.owin owin ppp
#' @importFrom spatstat.random runifpoint
#' @export
rcompppinhom <- function(n = 1,
                         lambda,
                         lambda_max,
                         nu = 1,
                         window = owin(),
                         simultaneous_draws = 1e5,
                         nquad = 1e5,
                         drop = TRUE,
                         use_cpp = TRUE) {
  quad <- runifpoint(n = nquad, win = window)
  area <- area.owin(window)

  result <- lapply(seq_len(n), function(i) {
    f <- function() {
      draws <- runifpoint(n = simultaneous_draws, win = window)
      lambda_draws <- lambda(draws$x, draws$y, i)
      cbind(draws$x, draws$y, lambda_draws)
    }

    integral_lambda <- mean(lambda(quad$x, quad$y, i), na.rm = TRUE) * area
    total_points <- COMPoissonReg::rcmp(n = 1, lambda = integral_lambda, nu = nu)
    if(use_cpp) {
      xyvector <- draw_points_cpp(total_points, f, lambda_max[i])
    } else {
      xyvector <- draw_points(total_points, f, lambda_max[i])
    }
    ppp(xyvector[, 1], xyvector[, 2], window = window)
  })
  if(length(result) == 1 && drop == TRUE) {
    result <- result[[1]]
  }
  result
}
