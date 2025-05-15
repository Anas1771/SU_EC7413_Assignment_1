# Assignment 1 - EC7413
## Group: Anas Haouat, Jonathan Holl & Andrea Jeppsson

# Start by loading the relevant libraries
library(dynlm)

# Load dataset (IncomeUK) from the Ecdat package (w/ downloading all of it)
data("IncomeUK", package = "Ecdat")

## Q1 - Estimate AR(1), AR(2), AR(3), AR(4), AR(5):

d_income <- diff(IncomeUK[,"income"]) # calculate the differences

d_income_m <- embed(d_income, dimension = 6) # create a matrix to estimate AR()

d_income_m <- ts(d_income_m, start = c(1971,2), frequency = 4)

# Run a loop:
AR1 <- dynlm(d_income_m[,1] ~ d_income_m[,2:2])
AR2 <- dynlm(d_income_m[,1] ~ d_income_m[,2:3])
AR3 <- dynlm(d_income_m[,1] ~ d_income_m[,2:4])
AR4 <- dynlm(d_income_m[,1] ~ d_income_m[,2:5])
AR5 <- dynlm(d_income_m[,1] ~ d_income_m[,2:6])

# Summaries
summary(AR1)
summary(AR2)
summary(AR3)
summary(AR4)
summary(AR5)


## Q2 -

