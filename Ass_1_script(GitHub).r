# Assignment 1 - EC7413
## Group: Anas Haouat, Jonathan Holl & Andrea Jeppsson


#--**--##--**--##--**--##--**--##--**--##--**--#
#--**--##--**--##--**--##--**--##--**--##--**--#


### Problem 1:

# Start by loading the relevant libraries
library(dynlm)
library(portes) # to perform Ljung-Box test
library(tseries) # to perform Jarque-Bera test


# Load dataset (IncomeUK) from the Ecdat package (w/ downloading all of it)
data("IncomeUK", package = "Ecdat")

## Q1 - Estimate AR(1), AR(2), AR(3), AR(4), AR(5):
d_income <- diff(IncomeUK[,"income"]) # calculate the differences
d_income_m <- embed(d_income, dimension = 6) # create a matrix to estimate AR()
d_income_m <- ts(d_income_m, start = c(1971,2), frequency = 4)

# Run a loop to get the 5 AR models:
for (i in 1:5) {
  assign( 
    sprintf("AR%d", i), 
    dynlm( d_income_m[,1] ~ d_income_m[, 2:(i+1)] ) 
  )
}

# Summaries
models <- list(AR1, AR2, AR3, AR4, AR5)
for (mod in models) {
  print(summary(mod))
}


#--**--##--**--##--**--##--**--##--**--##--**--#


## Q2 - Compute AIC for the 5 models:

# Compute the Akaike’s information criterion:
AIC(AR1, AR2, AR3, AR4, AR5)

# Answer - we choose AR(4), given that it has the smallest AIC.

## Q2 - Tests
# Ljung-Box:
for (m in models) {
  # compute lag‐choice + parameter count -
  # We use the rule of thumb where q = 0.75*n^(1/3), rounded down to the nearest integer:
  q     <- floor(0.75 * nobs(m)^(1/3)) + length(m$coefficients)
  param <- length(m$coefficients)
  
  # run and print the Ljung–Box test
  print( LjungBox(m$residuals, lags = q, fitdf = param) )
}

# Comment - we fail to reject H0 for AR4 and AR5, meaning that there is no
          # evidence of auto correlation in either, unlike AR1, AR2 and AR3.


# Jarque-Bera
for(mod in list(AR1, AR2, AR3, AR4, AR5)) print(jarque.bera.test(mod$residuals))

# Comment - we only reject the H0 with AR1, meaning that the other 4 models show
          # no evidence against normality and residuals are thus drawn from a 
          # normal distribution.

# Comment - Ljunx-Box and AIC:
          # Using both tests, we decide that AR4 is the most optimal model for
          # data. AR4 is not autocorrelated (according to LB) and presents the
          # best fit.


#--**--##--**--##--**--##--**--##--**--##--**--#
#--**--##--**--##--**--##--**--##--**--##--**--#


### Problem 2:

library(pxweb) # pxweb_get()
library(forecast) # auto.arima(), forecast()
library(dynlm) # dynlm()
library(sandwich) # NeweyWest()
library(lmtest) # coeftest()


#--**--##--**--##--**--##--**--##--**--##--**--#


# Q1 - download KPIF data from SCB using their API

# Use the pxweb package to download the KPIF data from SCB
pxweb_query_list <- 
  list("ContentsCode"=c("000005HS"),
       "Tid"=c("*"))

px_data <- 
  pxweb_get(url = "https://api.scb.se/OV0104/v1/doris/sv/ssd/PR/PR0101/PR0101G/KPIF",
            query = pxweb_query_list)

# Convert the data to a time series object
KPIF <- as.data.frame(px_data)
KPIF_column <- as.numeric(KPIF$"KPIF, månadsförändring, 1987=100")

# Make a time series from 1987M01 onward
full_ts  <- ts(KPIF_column,
               start     = c(1987,1),
               frequency = 12)

# Subset a time series that covers exactly 2002M12–2020M11 and gives 216 values:
KPIF_ts <- window(full_ts,
                  start = c(2002,12),
                  end   = c(2020,11))
rm(KPIF_column)
rm(full_ts)


#--**--##--**--##--**--##--**--##--**--##--**--#


# Q2 - forecast twice, using ARIMA and AR1 models:

# We have 216 quarterly observations between 2002:M12 and 2020:M11. We want to forecast for the 
# horizons h = 1,12,24 so we have 36 forecast origins.
origin_1 <- c(2002,12) # first forecast origin - 2002:M12
horizons <- c(1,12,24) # the 3 forecast horizons we are using
hmax <- max(horizons) # 24 steps ahead
n_origins <- length(window(KPIF_ts, start = origin_1))
freq <- 12 # monthly


# Create empty lists for both models:
ARIMA_models <- vector("list", n_origins)
ARIMA_forecasts <- vector("list", n_origins)

AR1_models <- vector("list", n_origins)
AR1_forecasts <- vector("list", n_origins)

# Create an empty matrix where we will store our forecasts
ARIMA_results <- ts(matrix(NA, nrow = n_origins, ncol = length(horizons)), 
                    start = origin_1,
                    frequency = frequency(KPIF_ts))

AR1_results <- ts(matrix(NA, nrow = n_origins, ncol = length(horizons)), 
                    start = origin_1,
                    frequency = frequency(KPIF_ts))

# Change the column names to match our 3 horizons:
colnames(ARIMA_results) <- paste0("h: ", horizons)
colnames(AR1_results) <- paste0("h: ", horizons)

# Forecast using ARIMA:
for (i in 1:n_origins) {
  origin_ARIMA <- origin_1 + (i-1)/freq
  y_ARIMA <- na.omit(window(KPIF_ts, end = origin_ARIMA))
  m_ARIMA <- auto.arima(y_ARIMA, seasonal = FALSE)
  f_ARIMA <- forecast(m_ARIMA, h = hmax)    # single integer
  ARIMA_models   [[i]] <- m_ARIMA
  ARIMA_forecasts[[i]] <- f_ARIMA
  ARIMA_results[i,]   <- f_ARIMA$mean[horizons]
}

# Forecast using AR1:
for (i in 2:n_origins) {
  origin_AR1   <- origin_1 + (i - 1)/freq
  y_AR1    <- na.omit(window(KPIF_ts, end = origin_AR1))
  m_AR1    <- Arima(y_AR1, order = c(1,0,0), method = "ML")
  f_AR1    <- forecast(m_AR1, h = hmax)
  AR1_models    [[i]] <- m_AR1
  AR1_forecasts [[i]] <- f_AR1
  AR1_results   [i, ] <- f_AR1$mean[horizons]
}

# Get rid of any redundant values:
rm(i,origin_ARIMA,origin_AR1,y_ARIMA,y_AR1,m_ARIMA,m_AR1,f_ARIMA,f_AR1)


#--**--##--**--##--**--##--**--##--**--##--**--#


#Q3 - Test the 2 forecasts for bias using the Newey-West estimator

# Our bias test involves regressing the forecast error on a constant and checking
# whether the estimated coefficient is significant.

# We already have:
ARIMA_results   # a time series of size n_origins × 3 with h = 1,12,24
AR1_results     # same dimensions

## Now we need to compute each models' errors since we need them for the bias:

# we start by creating an empty ts:
outcome <- ts(
  matrix(NA, nrow = n_origins, ncol = length(horizons)),
  start     = origin_1,
  frequency = frequency(KPIF_ts)
)
colnames(outcome) <- paste0("h = ", horizons)

# now we fill in via direct indexing:
for (i in 1:n_origins) {
  for (j in seq_along(horizons)) {
    h <- horizons[j]
    outcome[i, j] <- as.numeric(KPIF_ts[i + h])
  }
}

ARIMA_errors <- outcome - ARIMA_results
AR1_errors   <- outcome - AR1_results

# label the columns
colnames(ARIMA_errors) <- paste0("h = ", horizons)
colnames(AR1_errors)   <- paste0("h = ", horizons)

# done
rm(outcome)

## Compute bias

# Define a function to test bias using the regressing method:
bias_test <- function(errors, h) {
  model <- dynlm(errors[,h] ~ 1)
  test <- coeftest(model, vcov. = NeweyWest(model, lag = h -1 ))
  pval <- test["(Intercept)", "Pr(>|t|)"]
  return(pval)
}

# Create an empty dataframe to store the 6 p-values with adequate col names:
bias_pvals <- data.frame(c(1,12,24), ARIMA = NA, AR1 = NA)
colnames(bias_pvals) <- c("Horizons", "ARIMA", "AR1")


# Run a loop to fill the matrix with our p-values.
for (h in 1:3) {
  bias_pvals$ARIMA[h] <- bias_test(ARIMA_errors,h)
  bias_pvals$AR1[h] <- bias_test(AR1_errors,h)
}
rm(h)

# Results!

bias_pvals

# There’s no evidence that the 6 forecasts systematically are biased.
## DOUBLE CHECK!


#--**--##--**--##--**--##--**--##--**--##--**--#


#Q4 - Diebold-Mariano-test:
## Compare the two models for each horizon (3 comparisons in total):


# We create an empty dataframe where we will store the p-valuesfrom the
# Diebold-Mariano test for the 3 different forecast horizons and two models.

DM_pvals <- data.frame(horizons, pval = NA)

# Fill the dataframe with p-values from the test

for (i in seq_along(horizons)) {
  h    <- horizons[i]
  test <- dm.test(
    ARIMA_errors[, i],   # 1st col = h=1, 2nd = h=12, 3rd = h=24
    AR1_errors[,   i],
    h    = h
  )
  DM_pvals$pval[i] <- test$p.value
}

rm(i,j)

DM_pvals

# The null cannot only be rejected for horizon 1 at a 90% confidence level.

#DONE!