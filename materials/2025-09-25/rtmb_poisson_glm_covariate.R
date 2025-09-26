## Load TMB TRUE
library(RTMB)
library(tidyverse)

# read in the data
salmon <- read_csv("../2025-09-23/data/2025-ballard-sockeye-counts.csv") %>% 
  janitor::clean_names() %>% 
  mutate(date = as_date(date,format = "%d-%m"),
         month = month(date))
year(salmon$date) <- 2025
salmon

parameters <- list(log_beta=0)

# function to evaluate the objective function
f <- function(parms) {
  getAll(salmon, parms, warn=FALSE)
  ## Optional (enables extra RTMB features)
  count <- OBS(count)
  ## Data
  count %~% dpois(exp(log_beta))
  # ## Get predicted  uncertainties
  beta <- exp(log_beta)
  ADREPORT(beta)
  # ## Return
  # nll
}

#process the objective function
obj <- MakeADFun(f, parameters)

#fit the model using nlminb
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit
#TMBAIC(fit)

# get variance estimates for model parameters
sdr <- sdreport(obj)
sdr

summary(sdr) #gives the ADREPORT variables too.

# look at the residuals
#goodness of fit using One step ahead residuals
osa <- oneStepPredict(obj, 
                      method="oneStepGeneric",
                      discrete=TRUE,
                      range=c(0,Inf))
qqnorm(osa$res); abline(0,1)

## version with month covariate
parameters <- list(log_beta = 1:4) #rep(0,4))

# salmon$count
# salmon$month
# parameters$log_beta[salmon$month-5]

# function to evaluate the objective function
f2 <- function(parms) {
  getAll(salmon, parms, warn=FALSE)
  ## Optional (enables extra RTMB features)
  count <- OBS(count)
  ## Data
  # nll <- 0
  # for (i in 1:length(count))
  #   nll <- nll - dpois(count[i], exp(log_beta[as.integer(month[i]-5)]))
  # for (i in unique(month))
  #   nll <- nll - sum(dpois(count[month==i], exp(log_beta[i-5])))
  count %~% dpois(exp(log_beta[month-5]))
  # ## Get predicted  uncertainties
  beta <- exp(log_beta)
  ADREPORT(beta)
  # ## Return
  #nll
}

#process the objective function
obj <- MakeADFun(f2, parameters)
#fit the versionb with the common intercept 
# obj <- MakeADFun(f2, parameters, 
#                  map = list(log_beta=factor(c(1, 1, 1, 1))))

#fit the model using nlminb
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit
#TMBAIC(fit)

# get variance estimates for model parameters
sdr <- sdreport(obj)
sdr

summary(sdr) #gives the ADREPORT variables too.

# look at the residuals
#goodness of fit using One step ahead residuals
osa <- oneStepPredict(obj, 
                      method="oneStepGeneric",
                      discrete=TRUE,
                      range=c(0,Inf))
qqnorm(osa$res); abline(0,1)
