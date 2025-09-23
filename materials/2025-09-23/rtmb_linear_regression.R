# load RTMB

#install.packages('RTMB')
library(RTMB)
library(tidyverse) 
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
#library(TMBhelper)

#linear regression
# Generate some data
set.seed(8675309)
data <- tibble(x = 1:20,
               y = 0.5 + 2*x + rnorm(20,0,2))
# view the data set
data
ggplot(data) +
  geom_point(aes(x=x,y=y)) +
  geom_smooth(aes(x=x,y=y),method="lm") +
  theme_minimal()

#  yhat(i) = b0 + b1*x(i)
lm1 <- lm(y~x, data = data)
lm1

# model parameters
parameters <- list(b0=0, b1=0, logSigma=0)
print(parameters)

# function to evaluate the objective function
f <- function(parms) {
  getAll(data, parms, warn=FALSE)
  ## Optional (enables extra RTMB features)
  y <- OBS(y)
  ## Initialize negative log likelihood
  nll <- 0
  ## Data
  y_pred <- b0 + b1 * x
  nll <- nll - sum(dnorm(y, y_pred, sd=exp(logSigma), log=TRUE))
  ## Get predicted weight uncertainties
  ADREPORT(y_pred)
  ## Return
  nll
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
#print(summary(sdr, p.value = TRUE))
#print(summary(sdr, select = "fixed"))

#can also extract results as list objects
as.list(sdr, "Est") ## parameter estimates
as.list(sdr, "Std") ## parameter uncertainties

as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
as.list(sdr, "Std", report=TRUE) ## ADREPORT uncertainties

#goodness of fit using One step ahead residuals
osa <- oneStepPredict(obj, method="fullGaussian", discrete=FALSE)
qqnorm(osa$res); abline(0,1)

#plot the predictions
ggplot(data) +
  geom_point(aes(x=x,y=y)) +
  geom_abline(intercept = fit$par["b0"],
              slope = fit$par["b1"]) +
  theme_minimal()


# probabilistic syntax example
f2 <- function(parms) {
  getAll(data, parms, warn=FALSE)
  ## Optional (enables extra RTMB features)
  y <- OBS(y)
  # ## Initialize negative log likelihood
  # nll <- 0
  ## Data
  y_pred <- b0 + b1 * x
  # nll <- nll - sum(dnorm(y, y_pred, sd=exp(logSigma), log=TRUE))
  y %~% dnorm(y_pred, sd=exp(logSigma))
  ## Get predicted weight uncertainties
  ADREPORT(y_pred)
  ## Return
  # nll
}

#process the objective function
obj <- MakeADFun(f2, parameters)

#fit the model using nlminb
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit
