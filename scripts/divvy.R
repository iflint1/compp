library(compp)
library(raster)
library(spatstat)

lambda_from_fit <- function(fit, indices, covariates) {
  fit <- fit$fit
  lambda <- lapply(indices, function(n) {
    cov <- list(nu = function(x, y) rep(0, length(x)),
                os = function(x, y) rep(0, length(x)))
    if(length(covariates) > 0) {
      for(i in seq_len(length(covariates))) {
        current_covariate <- covariates[[i]][[n]]
        cov[[length(cov) + 1]] <- if(is(current_covariate, "RasterLayer")) {
          maptools::as.im.RasterLayer(current_covariate)
        } else {
          current_covariate
        }
        names(cov)[length(cov)] <- names(covariates)[i]
      }
    }

    predict.ppm(fit, covariates = cov, ngrid = c(64, 64), type = "lambda")
  })

  lambda
}

divvy_path <- "~/R/scripts/data/Divvy_Trips.csv"
elevation_path <- "~/R/scripts/data/Chicago_dem.tif"
temperature_path <- "~/R/scripts/data/chii2h2019.txt"

# Divvy spatio-temporal analysis
seed <- 1
set.seed(seed)
n_bootstrap <- 1000
nthreads <- 1
ndummy <- 1e3

divvy <- read.csv(divvy_path)
# This corresponds to the days that we consider.
days <- c("02/01", "02/02", "02/03", "02/04", "02/05", "02/06", "02/07", "02/08", "02/09",
          "02/10", "02/11", "02/12", "02/13", "02/14", "02/15", "02/16", "02/17", "02/18",
          "02/19", "02/20", "02/21", "02/22", "02/23", "02/24", "02/25", "02/26", "02/27",
          "02/28",
          "03/01", "03/02", "03/03", "03/04", "03/05", "03/06", "03/07", "03/08", "03/09",
          "03/10", "03/11", "03/12", "03/13", "03/14", "03/15", "03/16", "03/17", "03/18",
          "03/19", "03/20", "03/21", "03/22", "03/23", "03/24", "03/25", "03/26", "03/27",
          "03/28", "03/29", "03/30", "03/31",
          "04/01", "04/02", "04/03", "04/04", "04/05", "04/06", "04/07", "04/08", "04/09",
          "04/10", "04/11", "04/12", "04/13", "04/14", "04/15", "04/16", "04/17", "04/18",
          "04/19", "04/20", "04/21", "04/22", "04/23", "04/24", "04/25", "04/26", "04/27",
          "04/28", "04/29", "04/30",
          "05/01", "05/02", "05/03", "05/04", "05/05", "05/06", "05/07", "05/08", "05/09",
          "05/10", "05/11", "05/12", "05/13", "05/14", "05/15", "05/16", "05/17", "05/18",
          "05/19", "05/20", "05/21", "05/22", # The dates "05/23", "05/24", "05/25", "05/26", "05/27" contain 999 values for the temp
          "05/28", "05/29", "05/30", "05/31")

# Remove rows with NA values
divvy <- divvy[!is.na(divvy$FROM.LONGITUDE) & !is.na(divvy$FROM.LATITUDE), ]

# Small jitter to avoid duplicate points
divvy$FROM.LONGITUDE <- rnorm(length(divvy$FROM.LONGITUDE), mean = divvy$FROM.LONGITUDE, sd = 1e-6)
divvy$FROM.LATITUDE <- rnorm(length(divvy$FROM.LATITUDE), mean = divvy$FROM.LATITUDE, sd = 1e-6)

xrange <- c(min(divvy$FROM.LONGITUDE, na.rm = TRUE), max(divvy$FROM.LONGITUDE, na.rm = TRUE))
yrange <- c(min(divvy$FROM.LATITUDE, na.rm = TRUE), max(divvy$FROM.LATITUDE, na.rm = TRUE))

new_xrange <- xrange
new_yrange <- yrange
window <- owin(xrange, yrange)
configurations <- vector(mode = "list", length = length(days))
for(i in seq_len(length(days))) {
  df <- divvy[startsWith(as.vector(divvy$START.TIME), paste0(days[i], "/2019 11:00")) & endsWith(as.vector(divvy$START.TIME), "AM"), ]
  configurations[[i]] <- as.ppp(unique.ppp(ppp(df$FROM.LONGITUDE, df$FROM.LATITUDE, window = window, check = FALSE)))
}

remove(list = c("divvy", "df"))
gc()

elevation <- raster::raster(elevation_path)

contour(elevation, xlim = c(-87.7, -87.6), ylim = c(41.75, 42.1), nlevels = 10)
plot(configurations[[107]], add = TRUE, col = "red", pch = 20)

meteo <- read.table(temperature_path, header = TRUE)
temp_by_day <- sapply(seq_len(length(days)), function(i) {
  meteo[meteo$MM == as.integer(gsub("(.+)/.*", "\\1", days[i])) & meteo$DD == as.integer(gsub(".*/(.+)", "\\1", days[i])) & meteo$hh == 11 & meteo$mm == 0, ]$ATMP
})

nconfigurations <- length(configurations)

temp_by_day_raster <- lapply(temp_by_day, function(cov) {
  as.im(function(x, y) rep(cov, length(x)), W = window)
})
elevation_by_day_raster <- lapply(seq_len(nconfigurations), function(i) {
  elevation
})

covariates <- list(elevation = elevation_by_day_raster,
                   temperature = temp_by_day_raster)

set.seed(seed)
fit <- rcomfitlogit(configurations, covariates = covariates, ndummy = ndummy)
print(fit)

nu <- fit$coef[names(fit$coef) == 'nu']

lambdas_from_fit <- lambda_from_fit(fit = fit,
                                    #indices = seq_len(length(configurations)),
                                    indices = 107,
                                    covariates = covariates)

max_lambdas <- sapply(lambdas_from_fit, function(l) max(l))
lambdas_function <- lapply(lambdas_from_fit, function(l) as.function(l))
lambdas <- function(x, y, t) lambdas_function[[t]](x, y)

tm <- Sys.time()
set.seed(seed)
divvy_bootstrap_nu <- bootstrapinhom(N = n_bootstrap,
                                     n = length(configurations),
                                     estimate = fit$coef,
                                     window = window,
                                     lambda = lambdas,
                                     lambda_max = max_lambdas,
                                     covariates = covariates,
                                     nthreads = nthreads,
                                     ndummy = ndummy,
                                     parallel = FALSE)
Sys.time() - tm

print(divvy_bootstrap_nu)

# mppm version
H <- hyperframe(Y = configurations)
for(i in seq_len(length(covariates))) {
  H <- cbind(H, lapply(covariates[[i]], function(cov) if(is(cov, "RasterLayer")) {
    maptools::as.im.RasterLayer(cov)
  } else {
    cov
  }))
}
names(H) <- c("Y", names(covariates))

form <- as.formula(paste0("Y ~ 1 + ", paste0(names(covariates), collapse = '+')))

set.seed(seed)
fit_mppm <- spatstat::mppm(form, data = H)
summary_fit <- summary(fit_mppm)$Fit$FIT$coefficients

# Print MPPM fit results
print(summary_fit[, 1])
print(summary_fit[, 1] - 1.96 * summary_fit[, 2])
print(summary_fit[, 1] + 1.96 * summary_fit[, 2])

# Print predicted intensity + all points
# pred <- predict(fit_mppm, type = "trend", ngrid = c(1024, 1024))$trend
pred <- lambdas_from_fit
index <- 7
plot(pred[[index]])
points(Reduce(superimpose, configurations), col = 'white')
