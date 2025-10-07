## surplus production model

## Load TMB TRUE
library(RTMB)
library(tidyverse)
library(readxl)


## Data and parameters

catch <- read_xlsx("data/flounder_index_data.xlsx","catch")
index <- read_xlsx("data/flounder_index_data.xlsx","survey")

index_years <- index$Year - min(catch$Year) + 1

data <- list(catches = catch$Catch,
             index = index$Index,
             index_years = index_years)
data

parameters <- list(log_k = 10, 
                   log_r = log(0.6), 
                   log_m_est = log(1),  #log(1) here is schaefer model
                   log_sigma_proc = log(0.02),
                   log_sigma_c = log(0.05), 
                   log_sigma_i = log(0.2),
                   #logq = -2, 
                   log_b = rep(5,length(data$catches)),
                   logit_u = rep(-1.6, length(data$catches)))
parameters


pella_t <- function(b1, r, K, m, u_i) {
  u_use = exp(u_i)/(1.0+exp(u_i))
  b2 = b1 + b1*r*(1/(m-1))*(1.0-(b1/K)^(m-1.0)) - u_use*b1
  log_b2 = log(b2)
}

## 
prod_model_fun <- function(parameters) {
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
    
  log_bpred <- pella_t(exp(log_k), exp(log_r), exp(log_k), m, logit_u[1])
  nll <- nll - dnorm(log_b[1], log_bpred, exp(log_sigma_proc), TRUE)
  # for (year in 3:(length(catches)+1)) {
  #   log_bpred <- pella_t(exp(log_b[year-2]), exp(log_r), exp(log_k), m, logit_u[year-1])
  #   nll <- nll - dnorm(log_b[year-1], log_bpred, exp(log_sigma_proc), TRUE)
  # }
  for (year in 2:(length(catches))) {
    log_bpred <- pella_t(exp(log_b[year-1]), exp(log_r), exp(log_k), m, logit_u[year])
    nll <- nll - dnorm(log_b[year], log_bpred, exp(log_sigma_proc), TRUE)
  }
  
  biomass <- rep(0,length(catches)+1)
  biomass[1] <- exp(log_k)
  biomass[2:(length(catches)+1)] <- exp(log_b)
  # for (year in 2:(length(catches)+1))
  #   biomass[year] <- exp(log_b[year-1])
  
  #observation model
  #predicted survey index, estimate survey q analytically
  #pred_bio <- 0.5*(biomass[index_years]+biomass[index_years+1])
  pred_bio <- biomass[index_years]
  log_q <- sum(log(index/pred_bio))/length(index)
  index_pred = exp(log_q)*pred_bio
  
  ##predicted catches
  #catch_pred <- AD(rep(0, length(catches)))
  catch_pred <- biomass[-(length(catches)+1)]*exp(logit_u)/(1.0+exp(logit_u))
  # for (year in 1:length(catches))
  #   catch_pred[year] <- biomass[year]*exp(logit_u[year])/(1.0+exp(logit_u[year]))
  
  
  # Likelihood
  ## CATCHES
  nll <- nll - sum(dnorm(log(catches),log(catch_pred),exp(log_sigma_c), TRUE))
  ##SURVEY
  nll <- nll - sum(dnorm(log(index),log(index_pred),exp(log_sigma_i), TRUE));
  
  q <- exp(log_q)

  ADREPORT(biomass)
  ADREPORT(q)
  ADREPORT(m)
  
  nll
}

prod_model_fun(parameters)


## Make a function object
obj <- MakeADFun(prod_model_fun, parameters, 
                 random = c("log_b"),
                 map = list(log_m_est=factor(NA), log_sigma_c = factor(NA),#)) #,
                 #map = list(log_sigma_c = factor(NA),
                            log_sigma_proc = factor(NA))) #,
                 #control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr, control= list(eval.max = 10000,
                                                     iter.max = 10000))

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
summary(sdr)

pl <- obj$env$parList(opt$par) 

results <- tibble(year = 1935:2017) |>
  mutate(biomass = as.list(sdr, "Est", report = TRUE)$biomass,
         sd_bio = as.list(sdr, "Std", report = TRUE)$biomass) |>
 left_join(index %>% rename(year = Year)) |>
  mutate(lo = qlnorm(0.025, log(biomass), sdlog = sd_bio/biomass),
         hi = qlnorm(0.975, log(biomass), sdlog = sd_bio/biomass))
         
results |>
  ggplot() +
  aes(x = year, y = biomass) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = gray(0.3), alpha = 0.2) +
  geom_point(aes(y = Index/as.list(sdr, "Est", report = TRUE)$q),
              col = "darkgreen") +
  theme_bw()
