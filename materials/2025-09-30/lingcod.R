# schnute lingcod example

## Load TMB TRUE
library(RTMB)
library(tidyverse)


## Data(and parameters
lingcod <- read_table("lingcod.dat",skip=7)
data <-list(obs=log(lingcod$`#data`))

parameters <- list(logr=0, logit_beta=0, log_sigma_proc=0, log_sigma_obs=0,
                   eta=rep(0,nrow(lingcod)-1))
## Model function
lingcod <- function(parms) {
  #model goes here
}


## Make a function object
obj <- MakeADFun(lingcod, parameters, random = "eta") #need the random effects

## Call function minimizer
opt <-nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
summary(sdr)

