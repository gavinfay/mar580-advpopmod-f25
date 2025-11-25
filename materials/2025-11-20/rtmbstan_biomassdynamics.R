## surplus production model

## Load TMB TRUE
library(RTMB)
library(tidyverse)
library(readxl)
library(tmbstan)
library(shinystan)

## Data and parameters

#catch <- read_xlsx("data/flounder_index_data.xlsx","catch")
#index <- read_xlsx("data/flounder_index_data.xlsx","survey")

#index_years <- index$Year - min(catch$Year) + 1

source("wildebeest.R")

Ca[Ca==0] <- 0.0001
data <- list(catches = Ca,
             index = Nhat,
             index_years = y)
gnu_data <- data

parameters <- list(log_n0 = log(0.5),
                   log_k = log(6), 
                   log_r = log(0.2), 
                   log_m_est = log(1),  #log(1) here is schaefer model
                   log_sigma_proc = log(0.015),
                   log_q = 0,
                   log_sigma_c = log(0.001), 
                   log_sigma_i = log(0.2),
                   #logq = -2, 
                   log_b = rep(1,length(data$catches)),
                   logit_u = rep(-3, length(data$catches)))
parameters



## 
prod_model_fun <- function(data, parameters) {
  pella_t <- function(b1, r, K, m, u_i) {
    u_use = exp(u_i)/(1.0+exp(u_i))
    b2 = b1 + b1*r*(1/(m-1))*(1.0-(b1/K)^(m-1.0)) - u_use*b1
    log_b2 = log(b2)
  }
  require(RTMB)
  "[<-" <- ADoverload("[<-")
  
  #make things visible
  getAll(data, parameters, warn=FALSE)
  ## 
  index <- OBS(index)
  catches <- OBS(catches)  
  
  nll <- 0
  
  #keep m above 1
  m <- 1. + exp(log_m_est)
  #print(c("m",m))
    
  log_bpred <- pella_t(exp(log_n0), exp(log_r), exp(log_k), m, logit_u[1])
  #log_bpred <- pella_t(exp(log_n0), exp(log_r), exp(log_k), m, catches[1])
  nll <- nll - dnorm(log_b[1], log_bpred, exp(log_sigma_proc), TRUE)
  #print(c(1, nll, log_bpred))
  for (year in 2:(length(catches))) {
    log_bpred <- pella_t(exp(log_b[year-1]), exp(log_r), exp(log_k), m, logit_u[year])
    #log_bpred <- pella_t(exp(log_b[year-1]), exp(log_r), exp(log_k), m, catches[year])
    nll <- nll - dnorm(log_b[year], log_bpred, exp(log_sigma_proc), TRUE)
   # print(c(year, nll, log_bpred))
  }
  
  biomass <- rep(0,length(catches)+1)
  biomass[1] <- exp(log_n0)
  biomass[2:(length(catches)+1)] <- exp(log_b)
  # for (year in 2:(length(catches)+1))
  #   biomass[year] <- exp(log_b[year-1])
  
  #observation model
  #predicted survey index, estimate survey q analytically
  # pred_bio <- biomass[index_years]
  # log_q <- 0
  # for (i in 1:length(index))
  #  log_q <- log_q + log(index[i]/pred_bio[i])
  # log_q <- log_q/length(index)
  # index_pred = exp(log_q)*pred_bio
  index_pred = exp(log_q)*biomass[index_years]
  
  ##predicted catches
  #catch_pred <- AD(rep(0, length(catches)))
  catch_pred <- biomass[-(length(catches)+1)]*exp(logit_u)/(1.0+exp(logit_u))
  # for (year in 1:length(catches))
  #   catch_pred[year] <- biomass[year]*exp(logit_u[year])/(1.0+exp(logit_u[year]))
  
  
  # Likelihood
  # ## CATCHES
  # for (i in 1:length(catches))
  #  nll <- nll - dnorm(log(catches[i]),log(catch_pred[i]),exp(log_sigma_c), TRUE)
  # ##SURVEY
  # for (i in 1:length(index))
  #  nll <- nll - dnorm(log(index[i]),log(index_pred[i]),exp(log_sigma_i), TRUE)
  ## CATCHES
  nll <- nll - sum(dnorm(log(catches),log(catch_pred),exp(log_sigma_c), TRUE))
  ##SURVEY
  nll <- nll - sum(dnorm(log(index),log(index_pred),exp(log_sigma_i), TRUE))
  
    
  q <- exp(log_q)

  ADREPORT(biomass)
  ADREPORT(q)
  ADREPORT(m)
  
  nll
}

prod_model_fun(gnu_data, parameters)


cmb <- function(f, d) function(p) f(p, d)

# ## Make a function object
# obj <- MakeADFun(prod_model_fun, parameters, 
#                  random = c("log_b"),
#                  #map = list(log_m_est=factor(NA), log_sigma_c = factor(NA),
#                  map = list(log_sigma_c = factor(NA),
#                             log_sigma_proc = factor(NA))) #,
#                  #control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Make a function object
obj <- MakeADFun(cmb(prod_model_fun, gnu_data), parameters, 
                 random = c("log_b"),
                 map = list(log_m_est=factor(NA), log_sigma_c = factor(NA),
                            log_q = factor(NA),
                 #map = list(log_sigma_c = factor(NA),
                            log_sigma_proc = factor(NA))) #,
#control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr, control= list(eval.max = 10000,
                                                     iter.max = 10000))

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
summary(sdr)

#obj$simulate()

# set.seed(1) 
# chk <- checkConsistency(obj, estimate=TRUE)
# summary(chk)


stanfit <- tmbstan(obj, chains=1) #, init = opt$par)

cores <- parallel::detectCores()-1
options(mc.cores = cores)
# init.fn <- function()
#   list(u=rnorm(114), beta=rnorm(2), logsdu=runif(1,0,10), logsd0=runif(1,0,1))
stanfit <- tmbstan(obj, chains=3, open_progress=FALSE) #, init=list(fit$par, fit$par, fit$par))

## To explore the fit use shinystan
launch_shinystan(stanfit)


