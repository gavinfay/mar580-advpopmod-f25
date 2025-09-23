## Load TMB TRUE
library(RTMB)
library(tidyverse)

# read in the data
salmon <- read_csv("data/2025-ballard-sockeye-counts.csv") %>% 
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


