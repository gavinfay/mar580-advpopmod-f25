# schnute lingcod example

## Load TMB TRUE
library(RTMB)
library(tidyverse)


## Data(and parameters
lingcod <- read_table("data/lingcod.dat",skip=7)
data <-list(obs=log(lingcod$`#data`))

parameters <- list(logr=0, logit_beta=0, log_sigma_proc=0, log_sigma_obs=0,
                   log_b=rep(0,nrow(lingcod)-1))
## Model function
lingcod_model <- function(parms) {
  #make things visible
  getAll(data, parms, warn=FALSE)
  ## 
  obs <- OBS(obs)  
  
  #transform parameters
  sigma_proc = exp(log_sigma_proc)
  sigma_obs = exp(log_sigma_obs)
  beta = -2.0 + 2.0*exp(logit_beta)/(1.0+exp(logit_beta))
  r = exp(logr)
  
  #initialize section (storage variables, etc.)
  nll <- 0
  #process model
  #set initial state 
  # loop over time
  # predict expected biomass 
  # do process error model - compare prediction to state
    m = -r/beta #initial state
    nll = nll - dnorm(log_b[1], m, sigma_proc, TRUE)
    for(i in 2:(length(obs)-1)){
      m = r + (1+beta) * log_b[i-1] #// Gompertz
      nll <- nll - dnorm(log_b[i], m, sigma_proc, TRUE)
    }
    
    # #storage vector
    log_bpred <- rep(0, length(obs))
    print(class(log_bpred))
    log_bpred[1] = -r/beta
    print(class(log_bpred))
    for(i in 2:length(obs)){
      log_bpred[i] = log_b[i-1]
    }
    print(class(log_bpred))
    #// observation model
    nll <- nll - sum(dnorm(obs, log_bpred, sigma_obs, TRUE))

    #report_section
    ADREPORT(log_bpred)
    ADREPORT(r)
    ADREPORT(beta)
    ADREPORT(sigma_proc)
    ADREPORT(sigma_obs)
    nll
}


## Make a function object
obj <- MakeADFun(lingcod_model, parameters, random = "log_b") #need the random effects

## Call function minimizer
opt <-nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
summary(sdr)


biomass_pred <- tibble(mu = as.list(sdr, "Est", report = TRUE)$log_bpred,
                       sd = as.list(sdr, "Std", report = TRUE)$log_bpred) %>% 
  mutate(lo = qlnorm(0.025, mu, sdlog = sd),
         hi = qlnorm(0.975, mu, sdlog = sd),
         t = 1:33,
         y = exp(data$obs))
biomass_pred %>% 
  ggplot() +
  aes(x = t) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = gray(0.3), alpha = 0.2) +
  geom_point(aes(x = t, y=y), col = "blue") +
  geom_line(aes(y = exp(mu))) +
  theme_bw() + 
  labs(x = "time",
       y = "biomass",
       title = "lingcod state-space model") +
  ylim(0,1200)


# obj$simulate()
# 
set.seed(1)
chk <- checkConsistency(obj)
chk

chk <- checkConsistency(obj, estimate=TRUE)
summary(chk)
1-mean(summary(chk)$convergence). #very low convergence rate
GGally::ggpairs(slice(summary(chk)$estimate$par, which(summary(chk)$convergence==0)))
# poorly defined parameters
# this is not unexpected, our data do not traverse the production function AND we don't have multiple observations at a given time step


##################
# version 2 that assumes the states underlying the full observed time series are all random effects
##
parameters <- list(logr=0, logit_beta=0, log_sigma_proc=0, log_sigma_obs=0,
                   log_b=rep(0,nrow(lingcod))) #dimension of log_b is now equal to length of data
## Model function
lingcod_model_v2 <- function(parms) {
  #make things visible
  getAll(data, parms, warn=FALSE)
  ## 
  obs <- OBS(obs)  
  
  #transform parameters
  sigma_proc = exp(log_sigma_proc)
  sigma_obs = exp(log_sigma_obs)
  beta = -2.0 + 2.0*exp(logit_beta)/(1.0+exp(logit_beta))
  r = exp(logr)
  
  #initialize section (storage variables, etc.)
  nll <- 0
  #process model
  #set initial state 
  # loop over time
  # predict expected biomass 
  # do process error model - compare prediction to state
  m = -r/beta #initial state
  nll = nll - dnorm(log_b[1], m, sigma_proc, TRUE)
  for(i in 2:(length(obs))){
    m = r + (1+beta) * log_b[i-1] #// Gompertz
    nll <- nll - dnorm(log_b[i], m, sigma_proc, TRUE)
  }
  
  # # #storage vector
  # log_bpred <- rep(0, length(obs))
  # print(class(log_bpred))
  # log_bpred[1] = -r/beta
  # print(class(log_bpred))
  # for(i in 2:length(obs)){
  #   log_bpred[i] = log_b[i-1]
  # }
  # print(class(log_bpred))
  #// observation model
  nll <- nll - sum(dnorm(obs, log_b, sigma_obs, TRUE))
  
  #report_section
  #ADREPORT(log_bpred)
  ADREPORT(r)
  ADREPORT(beta)
  ADREPORT(sigma_proc)
  ADREPORT(sigma_obs)
  nll
}


## Make a function object
obj <- MakeADFun(lingcod_model_v2, parameters, random = "log_b") #need the random effects

## Call function minimizer
opt <-nlminb(obj$par, obj$fn, obj$gr)

opt$par #we notice different parameter estimates -including a reversal in the magnitude of the process andf observation error variance

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
summary(sdr)


biomass_pred <- tibble(mu = as.list(sdr, "Est")$log_b,
                       sd = as.list(sdr, "Std")$log_b) %>% 
  mutate(lo = qlnorm(0.025, mu, sdlog = sd),
         hi = qlnorm(0.975, mu, sdlog = sd),
         t = 1:33,
         y = exp(data$obs))
biomass_pred %>% 
  ggplot() +
  aes(x = t) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = gray(0.3), alpha = 0.2) +
  geom_point(aes(x = t, y=y), col = "blue") +
  geom_line(aes(y = exp(mu))) +
  theme_bw() + 
  labs(x = "time",
       y = "biomass",
       title = "lingcod state-space model") +
  ylim(0,1200)


# obj$simulate()
# 
set.seed(1)
chk <- checkConsistency(obj)
chk
# 
set.seed(1) 
chk <- checkConsistency(obj, estimate=TRUE)
summary(chk)
1-mean(summary(chk)$convergence) #this parameterization is more estimable, but the parameter estimates remain very undefined
GGally::ggpairs(slice(summary(chk)$estimate$par, which(summary(chk)$convergence==0)))

sim_param_ests <- summary(chk)$estimate$par |> 
  mutate(r = exp(logr),
         sigma_proc = exp(log_sigma_proc),
         sigma_obs = exp(log_sigma_obs),
         beta = -2.0 + 2.0*exp(logit_beta)/(1.0+exp(logit_beta))) |>
  select(-(1:4))
GGally::ggpairs(sim_param_ests)

