# age structured production model
library(RTMB)
library(tidyverse)
library(readxl)
library(compResidual)


## Data and parameters
catch <- read_xlsx("data/flounder_index_data.xlsx","catch")
index <- read_xlsx("data/flounder_index_data.xlsx","survey")

nages <- scan("data/floundah_biology.txt",n=1,skip=9)
maturity <- scan("data/floundah_biology.txt",n=nages,skip=11)
selex <- scan("data/floundah_biology.txt",n=nages,skip=13)
weight <- scan("data/floundah_biology.txt",n=nages,skip=15)

fish_acomp <- read_xlsx("data/floundah_agecomps.xlsx","fishery")[,-1]
surv_acomp <- read_xlsx("data/floundah_agecomps.xlsx","survey")[,-1]

#catch <- dplyr::filter(catch, Year >= 1963)
index_years <- index$Year - min(catch$Year) + 1

maturity[1] <- 0

data <- list(catches = catch$Catch, 
             index = index$Index,
             index_years = index_years,
             selex = selex,
             weight = weight,
             maturity = maturity,
             k=0,
             fishery_acomp = as.matrix(fish_acomp,nrow=nrow(fish_acomp)),
             survey_acomp = as.matrix(surv_acomp,nrow=nrow(surv_acomp)),
             effN_fish = 100,
             nage_fish = ncol(fish_acomp),
             effN_surv = 100,
             nage_surv = ncol(surv_acomp)) #,
#             nproj_yrs = 30,
#             fac2 =round(100*as.matrix(fish_acomp,nrow=nrow(fish_acomp))),
#             sac2 = round(100*as.matrix(surv_acomp,nrow=nrow(surv_acomp))))
#             projF = 0.3

parameters <- list(#dummy = 0,
  logR0 = log(10000),
  logfracR1 = log(1),
  logith = 0, # 0.5108256, # 0,
  logM = log(0.4),
  logSigmaC = log(0.05),
  logSigmaI = log(0.2),
  logSelpar = log(c(4,1,2,1)),
  logF = rep(log(0.2),nrow(catch)),
  recdev = rep(0.5*0.5*0.5,nrow(catch)),
  logSigmaR = log(0.5))


aspm_fun <- function(data, parameters) {
  "[<-" <- ADoverload("[<-")
  
  #make things visible
  getAll(data, parameters, warn=FALSE)
  ## 
  index <- OBS(index)
  catches <- OBS(catches)  

  nyrs <- length(catches)
  nages <- length(weight)
  
  #transform parameters
  Rzero = exp(logR0);
  R1 = Rzero*exp(logfracR1)
  h = 0.2 + 0.8*exp(logith)/(1.0+exp(logith))
  M = exp(logM)
  SigmaC = exp(logSigmaC)
  SigmaI = exp(logSigmaI)

  #Selectivity estimates
  Selpar = exp(logSelpar)
#  fselex <- selex
#  sselex <- selex
  fselex <- rep(0,nages)
  sselex <- rep(0,nages)
  for (age in 1:nages) {
    fselex[age] = 1./(1.+exp(-log(19)*(age+1.0-Selpar[1])/Selpar[2]));
    sselex[age] = 1./(1.+exp(-log(19)*(age+1.0-Selpar[3])/Selpar[4]));
  }
  
  nll <- 0
  
  
  #// F & Z calculations
  F <- matrix(0,nrow=nyrs,ncol=nages)
  Z <- matrix(0,nrow=nyrs,ncol=nages)
  for (iyr in 1:nyrs)
    F[iyr,] <- exp(logF[iyr])*fselex
  Z = M + F
  
  
  #/* POPULATION DYNAMICS */
  Bio <- rep(0, nyrs)
  SSB <- rep(0,nyrs+1)
  R <- rep(0,nyrs+1)
  N <- matrix(0,nrow=nyrs+1,ncol=nages)
  
  #//   Initial numbers at age & 
  #  //   unfished Spawning biomass calculation
  N[1,1] = R1
  for (age in 2:nages) {
    N[1,age] = N[1,age-1]*exp(-M)
  }
  N[1,nages] = N[1,nages]/(1-exp(-M))
  
  nb0 <- rep(0, nages)
  nb0[1] = Rzero
  for (age in 2:nages) 
    nb0[age] = nb0[age-1]*exp(-M)
  nb0[nages] = nb0[nages]/(1-exp(-M))
  SSB0 = sum(0.5*nb0*weight*maturity)
  
  #// 1st year stock quantities
  for (age in 1:nages) {
    Bio[1] = Bio[1] + sselex[age]*N[1,age]*weight[age]*exp(-0.5*Z[1,age]);
    if (age != 0) SSB[1] = SSB[1] + 0.5*N[1,age]*maturity[age]*weight[age];  
  }
  R[1] = N[1,1]
  
  #//   Begin Loop over time
  for (iyr in 2:(nyrs+1)) {
    
    #//     Numbers at age update 
    for (age in 2:nages)
      N[iyr,age] = N[iyr-1,age-1]*exp(-Z[iyr-1,age-1])
    N[iyr,nages] = N[iyr,nages] + N[iyr-1,nages]*exp(-Z[iyr-1,nages])
    
    #//Biomass predictions
    Bio[iyr] = 0.  #//bio is survey available biomass
    SSB[iyr] = 0.
    for (age in 1:nages) {
      if (iyr != (nyrs+1)) Bio[iyr] <- Bio[iyr] + sselex[age]*N[iyr,age]*weight[age]*exp(-0.5*Z[iyr,age])
      if (age != 1) SSB[iyr] <- SSB[iyr] + 0.5*N[iyr,age]*maturity[age]*weight[age]
    }
    
    #// Recruitment for this year given annual stock size
    R[iyr] <- (4*h*Rzero*SSB[iyr-k]/(SSB0*(1-h)+SSB[iyr-k]*(5*h-1)))
    R[iyr] <- R[iyr]*exp(recdev[iyr-1]-0.5*exp(logSigmaR)*exp(logSigmaR))
    N[iyr,1] <- R[iyr]
    
    #//       End Loop over time
  }

  
  #//       Contribution to the Likelihood function
  
  # estimate q analytically
  nobs <- length(index_years)
  index_pred <- rep(0,nobs)
  q = exp(sum(log(index/Bio[index_years]))/nobs)
  index_pred = q*Bio[index_years]
  #print(q)
  
  #  #estimate sigmaI analytically
  # SigmaI <- 0.;
  # // for (int iobs=0; iobs<nobs; iobs++) {
  #   //   SigmaI += pow(log(index(iobs))-log(index_pred(iobs)),2);
  #   // }
  # // SigmaI = sqrt(SigmaI/nobs);
  
  #catch predictions
  catch_pred <- rep(0, nyrs)
  for (iyr in 1:nyrs) {
    catch_pred[iyr] = sum(weight*F[iyr,]*N[iyr,]*(1.0-exp(-1.*Z[iyr,]))/Z[iyr,])
    # for (age in 1:nages)
    #   catch_pred[iyr] = catch_pred[iyr] + weight[age]*F[iyr,age]*N[iyr,age]*(1.0-exp(-1.*Z[iyr,age]))/Z[iyr,age]
  }
  
  REPORT(index_pred)
  REPORT(catch_pred)
  
  # //predictions of catch at age
  # matrix<Type> paa(nyrs,nage_fish);
  # paa.setZero();
  # // calculate predicted catch at ages in numbers
  p_caa <- matrix(0, nrow = nyrs, ncol = nage_fish)
  for (iyr in 1:nyrs) {
    temp_catage = F[iyr,]*N[iyr,]*(1.0-exp(-1.*Z[iyr,]))/Z[iyr,]
    p_caa[iyr,] = c(temp_catage[1:(nage_fish-1)], 
                    sum(temp_catage[nage_fish:nages]))/sum(temp_catage)
  }
  
  # //predictions of survey age comp
  p_saa <- matrix(0, nrow = nyrs, ncol = nage_surv)
  for (iyr in 1:nyrs) {
    temp_survage = N[iyr,]*exp(-0.5*Z[iyr,])*sselex
    p_saa[iyr,] = c(temp_survage[1:(nage_surv-1)], 
                    sum(temp_survage[nage_surv:nages]))/sum(temp_survage)
  }
  
  REPORT(p_caa)
  REPORT(p_saa)
  
  #// Likelihood
  #//CATCHES
  #print(catch_pred)
  #print(N[,1])
  nll <- nll - sum(dnorm(log(catches),log(catch_pred),SigmaC, TRUE))
  #print(nll)
  #//SURVEY
  nll <- nll - sum(dnorm(log(index),log(index_pred),SigmaI, TRUE))
  #print(nll)
  
  #e.g.
  # dmultinom(observed_n, effN, predicted_proportions, log = TRUE)
  #fishery age comp likelihood
  for (iyr in index_years[1]:nyrs) {
    temp_obs <- as.numeric(effN_fish*fish_acomp[iyr-index_years[1]+1,])
    #print(class(temp_obs))
    nll <- nll - dmultinom(temp_obs, #size = effN_fish, 
                           prob = p_caa[iyr,], log = TRUE)
  }
  #survey age comp likelihood
  for (iyr in index_years[1]:nyrs) {
    temp_obs <- as.numeric(effN_surv*surv_acomp[iyr-index_years[1]+1,])
    nll <- nll - dmultinom(temp_obs, #size = effN_surv,
                             prob = p_saa[iyr,], log = TRUE)
  }
  
  # 
  # // age-composition data contribution!
  #   // first calculate predictions
  # // then do likelihood
  # // LL = -effN*obs_p*log(pred_p)
  # 
  # //recruitment penalty
  # neglogL -= sum(dnorm(rdev,0,exp(logSigmaR),true));
  nll <- nll - sum(dnorm(recdev, 0, exp(logSigmaR), log = TRUE))
  # 
  ADREPORT(SSB)
  nll
}  


aspm_fun(data, parameters)


cmb <- function(f, d) function(p) f(p, d)

## Make a function object
obj <- MakeADFun(cmb(aspm_fun, data), parameters, 
                 #random = c("rdev"),
                 map = list(logM = factor(NA),
                            #logfracR1 = factor(NA),
                            #logith = factor(NA),
                            logSigmaC = factor(NA), #)) #,
                            logSigmaI = factor(NA),
                            logSigmaR = factor(NA)#,
                            #recdev = factor(rep(NA,nrow(catch)))
                            ))
#control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr, control= list(eval.max = 10000,
                                                     iter.max = 10000))

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
summary(sdr)



hist(sdr$gradient.fixed)

library(reshape2)
correlation_matrix <- cov2cor(sdr$cov)

cor_df <- melt(correlation_matrix)

ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + # Add white borders between tiles
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(min(cor_df$value), max(cor_df$value)),
                       space = "Lab", name="Correlation") +
  theme_minimal() +
  labs(title = "Correlation Matrix Heatmap",
       x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_fixed() # Ensures square tiles


# set.seed(1)
# chk <- checkConsistency(obj)
# chk

# look at the residuals


# #can also extract results as list objects
# as.list(sdr, "Est") ## parameter estimates
# as.list(sdr, "Std") ## parameter uncertainties
# 
# as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
# as.list(sdr, "Std", report=TRUE) ## ADREPORT uncertainties



fish_age_osa <- resMulti(data$fishery_acomp+1e-07, obj$report()$p_caa)
plot(fish_age_osa)
surv_age_osa <- resMulti(data$survey_acomp+1e-07, obj$report()$p_saa)
plot(surv_age_osa)

catch_at_age <- obj$report()$p_caa |>
  as.data.frame() |>
  mutate(year = 1:nrow(catch)) |>
  pivot_longer(cols = 1:5, names_prefix = "V", 
               names_to = "age",
               values_to = "pred")
temp <- data$fishery_acomp |>
  as.data.frame() 
names(temp) <- 1:5
temp <- temp |>
  mutate(year = index_years) |>
  pivot_longer(cols = 1:5, 
               names_to = "age",
               values_to = "obs")
catch_at_age <- catch_at_age |>
  right_join(temp)

catch_at_age |>
  #filter(year == 29) |>
  ggplot() +
  aes(x = age, y = obs, group = year) +
  #geom_line() +
  geom_area(aes(y = pred), fill = gray(0.8)) +
  geom_point() +
  geom_line() +
  facet_wrap(~year, dir = "v") +
  ylim(0,NA) +
  theme_bw()


catch_at_age <- obj$report()$p_saa |>
  as.data.frame() |>
  mutate(year = 1:nrow(catch)) |>
  pivot_longer(cols = 1:5, names_prefix = "V", 
               names_to = "age",
               values_to = "pred")
temp <- data$survey_acomp |>
  as.data.frame() 
names(temp) <- 1:5
temp <- temp |>
  mutate(year = index_years) |>
  pivot_longer(cols = 1:5, 
               names_to = "age",
               values_to = "obs")
catch_at_age <- catch_at_age |>
  right_join(temp)

catch_at_age |>
  #filter(year == 29) |>
  ggplot() +
  aes(x = age, y = obs, group = year) +
  #geom_line() +
  geom_area(aes(y = pred), fill = gray(0.8)) +
  geom_point() +
  geom_line() +
  facet_wrap(~year, dir = "v") +
  ylim(0,NA) +
  theme_bw()





