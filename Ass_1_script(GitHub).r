# Assignment 1 - EC7413
## Group: Anas Haouat, Jonathan Holl & Andrea Jeppsson

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

