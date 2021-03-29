library(compp)
library(ggplot2)
library(scales)
library(spatstat)

set.seed(1)

manaus <- read.csv("data-raw/manaus.csv")

range <- c(1903, 1992)

flood_years <- c(1904, 1908, 1909, 1913, 1918,
                 1920, 1921, 1922, 1944, 1953,
                 1955, 1971, 1972, 1973, 1975,
                 1976, 1982, 1989)


flood_months <- sapply(flood_years, function(year) {
  index <- which.max(manaus$value[manaus$time >= year & manaus$time < year + 1])
  as.integer((manaus$time[manaus$time >= year & manaus$time < year + 1][index] - year) * 12)
})

df <- data.frame(year = rep(seq(from = range[1], to = range[2], by = 1), each = 12),
                 month = rep(1:12, range[2] - range[1] + 1))
df$flood <- sapply(seq_len(nrow(df)), function(n) {
  index <- which(df$year[n] == flood_years)
  if(length(index) > 0) {
    df$month[n] == flood_months[index]
  } else {
    FALSE
  }
})

df$month <- factor(df$month)
levels(df$month) <- month.abb

df$flood <- factor(df$flood)
levels(df$flood) <- c("No flood", "Flood")

ggplot(df, aes(x = year, y = month, fill = flood)) +
  geom_tile() +
  scale_y_discrete(limits = month.abb[length(month.abb):1]) +
  scale_x_continuous(breaks = seq(from = range[1], to = range[2], by = 4)) +
  scale_fill_manual(name = "", values = c('grey', 'blue')) +
  xlab("Year") +
  ylab("Month") +
  ggtitle("Flooding of Rio Negro river by month") +
  theme_minimal()


configurations <- lapply(seq(from = range[1], to = range[2], by = 1), function(year) {
  if(year %in% flood_years) {
    index <- which.max(manaus$value[manaus$time >= year & manaus$time < year + 1])
    spatstat.geom::ppp(x = manaus$time[manaus$time >= year & manaus$time < year + 1][index] - year,
                  y = runif(1),
                  window = spatstat.geom::owin())
  } else {
    spatstat.geom::ppp(x = c(), y = c(), window = spatstat.geom::owin())
  }
})

set.seed(1)
fit <- rcomfitlogit(configurations, covariates = list(), ndummy = 1e4)
b <- bootstrap(N = 10000, n = length(configurations), estimate = fit$coef, nthreads = 4)

set.seed(1)
fit_ppp <- rcomfitlogit(configurations, covariates = list(), force_nu = 1, ndummy = 1e4)
b_ppp <- bootstrap(N = 10000, n = length(configurations), estimate = fit_ppp$coef, nthreads = 4, force_nu = 1)
