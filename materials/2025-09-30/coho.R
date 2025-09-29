## univariate coho salmon example, section 3.2.1

library(RTMB)

#generate some data
set.seed(1138)
num.times <- 40
Ricker.a <- 1.5
Ricker.b <- 0.0003
equilibrium <-  log(Ricker.a)/Ricker.b
cat("equilibrium=", round(equilibrium),"\n")

initial.n <- round(0.1*equilibrium)
juveniles <- estimates <- numeric(num.times)
juveniles[1] <- rpois(1,lambda = Ricker.a*initial.n*exp(-Ricker.b*initial.n))
CV.obs <- .30
obs.sd <- sqrt(log(CV.obs^2+1))
estimates[1] <- round(rlnorm(1,meanlog=log(juveniles[1]),sdlog=obs.sd))
for(i in 2:num.times) {
  juveniles[i] <- rpois(1,lambda=Ricker.a*juveniles[i-1]*exp(-Ricker.b*juveniles[i-1]))
  estimates[i] <-  round(rlnorm(1,meanlog=log(juveniles[i])-obs.sd^2/2,sdlog=obs.sd))
}
my.ylim <- range(c(juveniles,estimates))
plot(1:num.times,juveniles,xlab="Year",ylab="",ylim=my.ylim,type="l")
lines(1:num.times,estimates,lty=2,col=2)
legend(2,0.95*max(my.ylim),legend=c("State","Observation"),lty=1:2,col=1:2)


#create the data list
data <- list(y = estimates)

#parameters
parameters <- list(log_n0 = 0,
                   logit_alpha = 0,
                   log_beta = 0,
                   log_sigma = -1,
                   log_n = rep(0,length(data$y)-1))

# model
coho_model <- function(parms) {
  #make objects visible
  getAll(data, parms, warn=FALSE)
  ## 
  y <- OBS(y)  
  
  nll <- 0
  
  #transform parameters
  alpha <- 2.0*exp(logit_alpha)/(1.0+exp(logit_alpha)) #keeps alpha between 0 & 2
  beta <- exp(log_beta)
  sigma <- exp(log_sigma)
  
  #storage vectors
  n <- rep(0,length(estimates))
  lamda <- n
  
  #initialize the population
  n[1] = exp(log_n0);
  
  #process model
  for (t in 2:length(estimates)) {
    n[t] <- exp(log_n[t-1])  #population size at time t
    lamda[t] <- alpha*n[t-1]*exp(-beta*n[t-1]) #prediction of population size at time t
    nll <- nll - dpois(n[t],lamda[t], log = TRUE) #process error contribution
  }
  
  #observation model
  nll <- nll - sum(dnorm(log(y), log(n)-0.5*sigma*sigma, sigma, log = TRUE))
  
  # #report some stuff
  ADREPORT(alpha)
  ADREPORT(beta)
  ADREPORT(sigma)
  ADREPORT(n)
  #ADREPORT(lamda)
  
  #return
  nll
  
}


## Make a function object
obj <- MakeADFun(coho_model,
                 parameters, 
                 random = "log_n")

## Call function minimizer
opt <-nlminb(obj$par, obj$fn, obj$gr)
opt

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
summary(sdr)

#extract the predictions, compute approximate confidence intervals and plot
n_pred <- summary(sdr, "report") %>% as_tibble() %>% slice(-c(1:3)) %>% 
  mutate(cv = `Std. Error`/Estimate,
           sdlog = sqrt(log(1+cv^2)),
         obs = data$y,
         true = juveniles) %>% 
  rowid_to_column(var = "t")
n_pred$lo <- qlnorm(0.025, log(n_pred$Estimate), sdlog = n_pred$sdlog)
n_pred$hi <- qlnorm(0.975, log(n_pred$Estimate), sdlog = n_pred$sdlog)
n_pred %>% 
  ggplot() +
  aes(x = t) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = gray(0.3), alpha = 0.2) +
  geom_line(aes(y = Estimate)) +
  geom_point(aes(x = t, y=obs), col = "blue") +
  geom_line(aes(x = t, y = true), col = "orange") +
  theme_bw() + 
  labs(x = "time",
       y = "number of coho",
       caption = "Model-estimated numbers (black line) & 95% CI (gray shading)\n compared to data (blue points), and true values used to generate data (orange)")

