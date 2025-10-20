
  library(RTMB)
  library(tidyverse)
  library(patchwork)
  
  ssl_data <- readRDS("../../2025-10-02/data/ssl-data.rds")
  
  forrester_pups <- ssl_data$forrester_pups
  forrester_nonpups <- ssl_data$forrester_nonpups
  
  ##
  fyr <- min(forrester_pups$Year, forrester_nonpups$Year)
  #fyr <- 1955
  pups <- forrester_pups %>% 
    mutate(Year = Year - fyr + 1)
  nonpups <- forrester_nonpups %>% 
    mutate(Year = Year - fyr + 1)
  
  data <- list(pups = log(pups$pupcount),
               pupyrs = as.integer(pups$Year),
               nonpups = log(nonpups$adultcount),
               nonpupyrs = as.integer(nonpups$Year))
  data$lyr <- max(c(data$pupyrs,data$nonpupyrs))
  data$lyr <- max(c(data$pupyrs,data$nonpupyrs)) + 12
  data
  lyr <- data$lyr
  
  
  parameters <- list(log_n0 = log(10000),
                     logit_beta = boot::logit(c(0.25,0.6,0.9,0.3)),
                     log_sigma = log(c(0.1,0.1)),
                     log_tau = log(rep(0.1,2)),
                     log_pups = rep(7,lyr),
                     log_nonpups = rep(7.5,lyr))
  parameters

  
  
  ssl_brs <- function(parameters) {
    "[<-" <- ADoverload("[<-")
    
    #make things visible
    getAll(data, parameters, warn=FALSE)
    ## 
    pups <- OBS(pups)
    nonpups <- OBS(nonpups)  
    
    nobs_pups = length(pups)
    nobs_nonpups = length(nonpups)
    
    beta_par = exp(logit_beta)/(1.+exp(logit_beta))
    sigma = exp(log_sigma)
    tau = exp(log_tau)
    #N0 = exp(log_n0)
    
    N <- matrix(0, nrow = 2, ncol = lyr)
    
    nll = 0
    pnp = 0
    
    #//population in first year
    nll <- nll - dnorm(log_pups[1], log(beta_par[1])+log_n0, tau[1], TRUE);
    #nll <- nll - dnorm(log_pups[1], log_n0[1], tau[1], TRUE);
    N[1,1] <- exp(log_pups[1])
    pnp <- log((beta_par[2]*beta_par[1] + beta_par[3])) + log_n0
    #pnp = log(beta_par[2]*exp(log_n0[1]) + beta_par[3]*exp(log_n0[2]))
    nll <- nll - dnorm(log_nonpups[1], pnp, tau[2], TRUE)
    N[2,1] <- exp(log_nonpups[1])

    for (t in 2:lyr) {
      #//pups
      nll <- nll - dnorm(log_pups[t], log(beta_par[1])+log(N[2,t-1]), tau[1], TRUE);
      N[1,t] = exp(log_pups[t])
      #//nonpups
      pnp = log(beta_par[2]*N[1,t-1] + beta_par[3]*N[2,t-1])
      nll <- nll - dnorm(log_nonpups[t], pnp, tau[2], TRUE)
      N[2, t] = exp(log_nonpups[t])
    }
    
    #//observation model
    pup_pred <- rep(0,nobs_pups)
    for (i in 1:nobs_pups)
      pup_pred[i] = log(N[1,pupyrs[i]])
    nonpup_pred <- rep(0,nobs_nonpups)
    for (i in 1:nobs_nonpups)
      nonpup_pred[i] = log(beta_par[4]*N[2,nonpupyrs[i]]);
    
    nll <- nll - sum(dnorm(pups, pup_pred, sigma[1], log = TRUE));
    nll <- nll - sum(dnorm(nonpups, nonpup_pred, sigma[2], log = TRUE));
    
    #//pup survival rate penalty
    nll <- nll - dbeta(beta_par[2], 57., 38., log = TRUE)
    #//haulout fraction penalty
    nll <- nll - dbeta(beta_par[4], 6., 14., log = TRUE)

    pup_nums <- N[1,]
    nonpup_nums <- N[2,]
    
    REPORT(pup_pred)
    REPORT(nonpup_pred)
    
    ADREPORT(pup_nums);
    ADREPORT(nonpup_nums);
    ADREPORT(beta_par);
    ADREPORT(sigma);
    ADREPORT(tau);
    
    nll
    
  }

  
  ssl_brs(parameters)
  

  # ## Make a function object
  # cmb <- function(f, d) function(p) f(p, d)
  
  parameters <- list(#log_n0 = log(c(2000,5000)),
                     log_n0 = log(10000),
                     logit_beta = boot::logit(c(0.25,0.75,0.9,0.2)),
                     log_sigma = log(c(0.01,0.05)),
                     log_tau = log(c(0.1,0.1)),
                     log_pups = rep(8,lyr),
                     log_nonpups = rep(7.5,lyr))
  parameters
  
  ## Make a function object
  obj <- MakeADFun(ssl_brs, parameters, 
                   random = c("log_pups","log_nonpups"),
                   map = list(logit_beta = factor(c(1,2,3,4)),
                              log_sigma = factor(c(1,2)),  #obs err estimated, unequal
                              #log_sigma = factor(c(1,1)), #obs err estimated, equal
                              #log_sigma = factor(c(NA,NA)), #obs err fixed
                              #log_sigma = factor(c(NA,1)), #obs err for pups fixed, nonpups estimated
                              #log_tau = factor(c(1,1))), #process err estimated, equal
                              log_tau = factor(c(1,2))), #process err estimated, unequal
                              #log_tau = factor(c(NA,1))), #process err for pups fixed, nonpups estimated
                              #log_tau = factor(c(NA,NA))), #process err fixed for both
  #)
                   control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  opt
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  sdr
  summary(sdr)

## look at fitted values
report <- obj$report()
plot(report$pup_pred, data$pups-report$pup_pred,
     xlab = "fitted values", ylab = "raw residuals", main = "pup counts")
abline(0,0, lty=2)
plot(report$nonpup_pred, data$nonpups-report$nonpup_pred,
     xlab = "fitted values", ylab = "raw residuals", main = "non-pup counts")
abline(0,0, lty=2)

  
  
#PLOT THE PREDICTED NUMBERS  
    
results <-  as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
std.results <- as.list(sdr, "Std", report=TRUE) ## ADREPORT uncertainties
  
results$year <- seq(from = fyr, by = 1, length.out = lyr)
results2 <- tibble(year = results$year,
                   pup_nums = results$pup_nums,
                   nonpup_nums = results$nonpup_nums,
                   lo_pups = qlnorm(0.025, log(pup_nums), sdlog = std.results$pup_nums/pup_nums),
                   hi_pups = qlnorm(0.975, log(pup_nums), sdlog = std.results$pup_nums/pup_nums),
                   lo_nonpups = qlnorm(0.025, log(nonpup_nums), 
                                       sdlog = std.results$nonpup_nums/nonpup_nums),
                   hi_nonpups = qlnorm(0.975, log(nonpup_nums),
                                       sdlog = std.results$nonpup_nums/nonpup_nums))
                  
p1 <- results2 |>
  ggplot() +
  aes(x = year, y = pup_nums) +
  geom_line() +
  geom_ribbon(aes(ymin = lo_pups, ymax = hi_pups), fill = gray(0.3), alpha = 0.2) +
  geom_point(data = forrester_pups,aes(x = Year, y = pupcount)) +
  ylim(0,NA)
p2 <- results2 |>
  ggplot() +
  aes(x = year, y = nonpup_nums) +
  geom_line() + 
  geom_ribbon(aes(ymin = lo_nonpups, ymax = hi_nonpups), fill = gray(0.3), alpha = 0.2) +
  geom_point(data = forrester_nonpups,aes(x = Year, y = adultcount/results$beta_par[4])) +
  ylim(0,NA)
p1 + p2


#goodness of fit using One step ahead residuals
osa <- oneStepPredict(obj, method="fullGaussian", discrete=FALSE)
qqnorm(osa$res); abline(0,1)

set.seed(42)
chk <- checkConsistency(obj, estimate=TRUE)
summary(chk)

summary(chk)$estimate$par |> 
  t() |>
  as.data.frame() |>
  rownames_to_column() |>
  pivot_longer(-1) |>
  ggplot() +
  aes(x = value) +
  geom_histogram(col="white") +
  facet_wrap(~rowname, scales = "free") +
  theme_bw()

summary(chk)$estimate$par |> 
  as_tibble(.name_repair = "unique") |>
  janitor::clean_names() |>
  GGally::ggpairs()


####### observation error only model
parameters <- list(log_n0 = log(10000),
                   logit_beta = boot::logit(c(0.25,0.6,0.9,0.3)),
                   log_sigma = log(c(0.2,0.2)),
                   log_tau = log(c(0.001,0.001)),
                   log_pups = rep(7,lyr),
                   log_nonpups = rep(7.5,lyr))
parameters

## Make a function object
obj <- MakeADFun(ssl_brs, parameters, 
                 random = c("log_pups","log_nonpups"),
                 map = list(logit_beta = factor(c(1,2,3,4)),
                            log_sigma = factor(c(1,2)),  #obs err estimated, unequal
                            log_tau = factor(c(NA,NA))), #process err fixed for both
#)
                control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
summary(sdr)


####### process error only model
parameters <- list(log_n0 = log(10000),
                   logit_beta = boot::logit(c(0.25,0.6,0.9,0.3)),
                   log_sigma = log(c(0.001,0.001)),
                   log_tau = log(c(0.1,0.1)),
                   log_pups = rep(8,lyr),
                   log_nonpups = rep(7.5,lyr))
parameters

## Make a function object
obj <- MakeADFun(ssl_brs, parameters, 
                 random = c("log_pups","log_nonpups"),
                 map = list(logit_beta = factor(c(1,2,3,4)),
                            log_sigma = factor(c(NA,NA)), #obs err fixed
                            log_tau = factor(c(1,2))), #process err estimated, unequal
                 #log_tau = factor(c(NA,1))), #process err for pups fixed, nonpups estimated
                 #log_tau = factor(c(NA,NA))), #process err fixed for both
)
#                control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
summary(sdr)



# #### getting time series of predictions from stochastic projections
# # the below is wrong because the historical time series are drawn from the statistical model, not from the posterior
# 
# # first refit the original model
# parameters <- list(
#   log_n0 = log(10000),
#   logit_beta = boot::logit(c(0.25,0.75,0.9,0.2)),
#   log_sigma = log(c(0.01,0.05)),
#   log_tau = log(c(0.1,0.1)),
#   log_pups = rep(8,lyr),
#   log_nonpups = rep(7.5,lyr))
# parameters
# ## Make a function object
# obj <- MakeADFun(ssl_brs, parameters, 
#                  random = c("log_pups","log_nonpups"),
#                  map = list(logit_beta = factor(c(1,2,3,4)),
#                             log_sigma = factor(c(1,2)),  #obs err estimated, unequal
#                             log_tau = factor(c(1,2))), #process err estimated, unequal
#                  control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))
# ## Call function minimizer
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# opt
# nsim = 100
# #simulate data sets, store the time series of random effects
# my_sims <- map(1:nsim, function(x) {set.seed(x);sim <- obj$simulate()})
# #grab the terminal year numbers of pups, compare to historical average count
# logpups_2030 <- map_dbl(my_sims, ~pluck(.x$log_pups[length(.x$log_pups)]))
# pups_2030 <- exp(logpups_2030)
# hist(pups_2030, breaks = 10)
# abline(v=mean(exp(data$pups)),lty=2, col = "blue", lwd = 3)
# 
# plot(1:71, exp(my_sims[[1]]$log_pups), type = 'l', col = gray(0.5), alpha = 0.05,
#      ylim = c(0,15000))
# for (i in 2:100)
#   lines(1:71, exp(my_sims[[i]]$log_pups), type = 'l', col = gray(0.5), alpha = 0.05)


#as.list(sdr, "Est")

