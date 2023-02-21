library(compp)
library(dplyr)
library(geodata)
library(rgdal)
library(rnaturalearth)
library(sf)
library(spatstat)
library(terra)


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

    predict.ppm(fit, covariates = cov, ngrid = c(512, 512), type = "lambda")
  })

  lambda
}

seed <- 1
set.seed(seed)
use_pixellated_window <- FALSE
use_yearly <- TRUE
ndummy <- 1e4
nthreads <- 1
n_bootstrap <- 1000
#window <- spatstat::owin(c(143, 147), c(-39.5, -36))
# window <- spatstat::owin(c(134, 150), c(-45, -30))
window <- owin(c(145.5, 150), c(-39.5, -34))
begin <- 2012
end <- 2018

path_locations <- "~/pCloudDrive/UoMSave/R/CWdata_AU_SA.csv"
path_bias <- "~/pCloudDrive/UoMSave/Downloads/bias.tif"
path_phosphorus <- "~/pCloudDrive/UoMSave/Downloads/phosphorus.tif"
path_pH <- "~/pCloudDrive/UoMSave/Downloads/pH.tif"
path_tmax <- "~/Downloads/wc2.1_2.5m_tmax_2010-2018/wc2.1_2.5m_tmax_"
path_prec <- "~/pCloudDrive/UoMSave/Downloads/wc2.1_2.5m_prec_2010-2018/au_wc2.1_2.5m_prec_"
path_save <- "~/pCloudDrive/UoMSave/Downloads/fit.RData"
dat <- read.csv(path_locations)

dat <- dat[dat$year >= begin & dat$year <= end & dat$country == "E_AUS", ]
dat$date <- paste0(dat$month, "/", dat$year)

month_from_date <- function(date) {
  as.numeric(gsub("(1?[0-9])/[0-9]+", "\\1", date))
}
year_from_date <- function(date) {
  as.numeric(gsub("1?[0-9]/([0-9]+)", "\\1", date))
}

if(use_yearly) {
  dates <- unique(year_from_date(dat$date))
} else {
  dates <- unique(dat$date)
}

# Small jitter to avoid duplicate points
dat$lon <- rnorm(length(dat$lon), mean = dat$lon, sd = 1e-6)
dat$lat <- rnorm(length(dat$lat), mean = dat$lat, sd = 1e-6)

# Not ideal, but the covariates have a given range of values, so get rid of outside points.
dat <- dat[!(dat$lon < 120 | dat$lon > 150 | dat$lat < -60 | dat$lat > -30), ]

# Also get rid of points outside our window
dat <- dat[!(dat$lon < window$xrange[1] | dat$lon > window$xrange[2] | dat$lat < window$yrange[1] | dat$lat > window$yrange[2]), ]

# Remove duplicates in data
dat <- dat[!duplicated(dat), ]

dat_by_date <- lapply(dates, function(date) {
  if(use_yearly) {
    dat[dat$year == date, ]
  } else {
    dat[dat$date == date, ]
  }
})
names(dat_by_date) <- dates

coordinates_by_date <- lapply(dat_by_date, function(dd) {
  d <- data.frame(lon = dd$lon, lat = dd$lat)
  if(nrow(d) > 0) {
    lon_limits <- c(min(d$lon), max(d$lon))
    lat_limits <- c(min(d$lat), max(d$lat))

    coordinates(d) <- c("lon", "lat")
    proj4string(d) <- CRS("+init=epsg:4326")
    #d_proj <- sp::spTransform(d, CRS("+init=epsg:3111"))
    d
  } else {
    d
  }
})

full_coordinates <- data.frame(lon = Reduce(c, sapply(coordinates_by_date, function(coord) coord$lon)),
                               lat = Reduce(c, sapply(coordinates_by_date, function(coord) coord$lat)))

# window <- spatstat::owin(c(min(full_coordinates$lon), max(full_coordinates$lon)),
#                                c(min(full_coordinates$lat), max(full_coordinates$lat)))

if(use_pixellated_window) {
  m <- spatstat::pixellate(spatstat::ppp(x = full_coordinates$lon,
                                         y = full_coordinates$lat,
                                         window = window), dimyx = c(13, 13)) > 0
  m$v[!m$v] <- NA
  window <- spatstat::as.owin(m)
}

# Load some of the covariates to construct window properly
phosphorus <- raster::raster(path_phosphorus)
crs(phosphorus) <- CRS("+init=epsg:4326")

pH <- raster::raster(path_pH)
crs(pH) <- CRS("+init=epsg:4326")

bias <- raster::raster(path_bias)
crs(bias) <- CRS("+init=epsg:4326")

construct_mask_window <- function(r) {
  pixel <- !is.na(r)
  msk <- maptools::as.im.RasterLayer(pixel)
  msk$v[!msk$v] <- NA
  as.owin(msk)
}

window <- intersect.owin(window,
                         construct_mask_window(phosphorus),
                         construct_mask_window(pH),
                         construct_mask_window(bias))

configurations <- lapply(coordinates_by_date, function(coordinates) {
  ppp(x = coordinates$lon, y = coordinates$lat, window = window)
})

australia <- sf::st_crop(rnaturalearth::ne_countries(country = "Australia",
                                         scale = 50,
                                         returnclass = "sf"), c(xmin = 110, xmax = 155, ymin = -45, ymax = -5))


states <- c("South Australia",
            "New South Wales",
            "Victoria",
            "Tasmania")
states <- rnaturalearth::ne_states(country = "Australia",
                                   returnclass = "sf") %>%
  dplyr::filter(name %in% koala_states) %>%
  sf::st_union() %>%
  terra::vect()

tmax <- geodata::worldclim_global(var = "tmax",
                                  res = 10,
                                  path = "data-raw") %>%
  terra::crop(states) %>%
  terra::mask(states)


names(tmax) <- substr(names(tmax), 11, 100)

plot(tmax)

# Covariates
covariates <- list(pH = lapply(seq_len(length(configurations)), function(i) pH),
                   #phosphorus = lapply(seq_len(length(configurations)), function(i) phosphorus),
                   # is_spring = lapply(seq_len(length(configurations)), function(i) {
                   #   month <- month_from_date(names(configurations)[i])
                   #   val <- (1 + cos(2 * pi * (month - 10) / 11)) / 2
                   #   function(x, y) rep(val, length(x))
                   # }),
                   bias = lapply(seq_len(length(configurations)), function(i) bias),
                   # prec = lapply(seq_len(length(configurations)), function(i) {
                   #   if(use_yearly) {
                   #     month <- 1
                   #     year <- names(configurations)[i]
                   #   } else {
                   #     month <- month_from_date(names(configurations)[i])
                   #     year <- year_from_date(names(configurations)[i])
                   #   }
                   #   file_path <- paste0(path_prec,
                   #                       year,
                   #                       "-",
                   #                       formatC(month, width = 2, format = "d", flag = "0"),
                   #                       ".tif")
                   #   raster::raster(file_path)
                   # }),
                   tmax = lapply(seq_len(length(configurations)), function(i) {
                       if(use_yearly) {
                         month <- 10
                         year <- names(configurations)[i]
                       } else {
                         month <- month_from_date(names(configurations)[i])
                         year <- year_from_date(names(configurations)[i])
                       }
                     file_path <- paste0(path_tmax,
                                         year,
                                         "-",
                                         formatC(month, width = 2, format = "d", flag = "0"),
                                         ".tif")
                     raster::raster(file_path)
                   })
)

a1 <- sapply(configurations, function(configuration) sum(configuration$y <= mean(window$yrange)))
a2 <- sapply(configurations, function(configuration) sum(configuration$y > mean(window$yrange)))
print(cor.test(a1, a2))

set.seed(seed)
fit <- rcomfitlogit(configurations, covariates = covariates, ndummy = ndummy)
print(fit)

nu <- fit$coef[names(fit$coef) == 'nu']

lambdas_from_fit <- lambda_from_fit(fit = fit,
                                    indices = seq_len(length(configurations)),
                                    covariates = covariates)

# Plot for manuscript
plot(lambdas_from_fit[[6]])
points(configurations[[6]], col = 'green')
points(configurations[[1]])
points(configurations[[2]])
points(configurations[[3]])
points(configurations[[4]])
points(configurations[[5]])
points(configurations[[7]])


# We need to slightly shrink the window to where the lambdas are non NA
msk <- lambdas_from_fit[[1]]
window <- intersect.owin(window, spatstat::as.owin(msk))

max_lambdas <- sapply(lambdas_from_fit, function(l) max(l))
lambdas_function <- lapply(lambdas_from_fit, function(l) as.function(l))
lambdas <- function(x, y, t) lambdas_function[[t]](x, y)

# draw <- rcompppinhom(n = length(configurations),
#                      lambda = lambdas,
#                      lambda_max = max_lambdas,
#                      nu = nu,
#                      window = window)

tm <- Sys.time()
set.seed(seed)
capeweed_bootstrap_nu <- bootstrapinhom(N = n_bootstrap,
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

print(capeweed_bootstrap_nu)

# mppm version
H <- hyperframe(Y = configurations)
for(i in seq_len(length(covariates))) {
  H <- cbind(H, lapply(covariates[[i]], function(cov) maptools::as.im.RasterLayer(cov)))
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
