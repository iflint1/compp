#' Homogeneous COM-Poisson point process simulation
#'
#' @param n Number of replications.
#' @param lambda Quasi-intensity.
#' @param nu Dispersion parameter
#' @param window Observation window.
#' @param drop Logical. If n = 1 and drop = TRUE (the default), the result will be a point pattern,
#' rather than a list containing a single point pattern.
#'
#' @return List of `spatstat::ppp` samples of a homogeneous COM-Poisson point process.
#'
#' @import COMPoissonReg
#' @importFrom spatstat.geom area.owin owin
#' @importFrom spatstat.random runifpoint
#' @export
rcomppp <- function(n = 1,
                    lambda = 1,
                    nu = 1,
                    window = owin(),
                    drop = TRUE) {
  area <- area.owin(window)
  total_points <- COMPoissonReg::rcmp(n = n, lambda = lambda * area, nu = nu)
  result <- lapply(total_points, function(npoints) {
    runifpoint(n = npoints, win = window)
  })
  if(length(result) == 1 && drop == TRUE) {
    result <- result[[1]]
  }
  result
}
