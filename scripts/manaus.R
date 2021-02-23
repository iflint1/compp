library(spatstat)

manaus <- read.csv("data-raw/manaus.csv")

range <- c(1903, 1992)

crude_flood_years <- c(1904, 1908, 1909, 1913, 1918,
                       1920, 1921, 1922, 1944, 1953,
                       1955, 1971, 1972, 1973, 1975,
                       1976, 1982, 1989)


flood_years <- sapply(crude_flood_years, function(year) {
  index <- which.max(manaus$value[manaus$time >= year & manaus$time < year + 1])
  manaus$time[manaus$time >= year & manaus$time < year + 1][index]
})

plot(flood_years, rep(1, length(flood_years)))

configurations <- lapply(seq(from = range[1], to = range[2], by = 1), function(year) {
  if(year %in% crude_flood_years) {
    spatstat::ppp(x = flood_years[flood_years >= year & flood_years < year + 1] - year,
                  y = runif(1),
                  window = spatstat::owin())
  } else {
    spatstat::ppp(x = c(), y = c(), window = spatstat::owin())
  }
})

fit <- rcomfitlogit(configurations, covariates = list(), ndummy = 1e4)
