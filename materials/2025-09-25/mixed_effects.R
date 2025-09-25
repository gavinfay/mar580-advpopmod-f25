
library(RTMB)
library(tidyverse)
library(lattice)
library(nlme)
library(lme4)
library(broom.mixed)

### STREAMS EXAMPLE

Streams <- scan('data/streams.dat',what=list(Stream=0,Density=0),n=3*18,skip=5)
streams.lm1 <- lm(Density ~ 1, data = Streams)
streams.lm2 <- lm(Density ~ factor(Stream) - 1, data = Streams)
streams.lme1 <- lmer(Density ~ 1 + (1 | Stream), data = Streams, REML=FALSE)


#summary of results
summary(streams.lme1)

#compare residuals among models
par(las=0)
boxplot(split(residuals(streams.lm1)/summary(streams.lm1)$sigma,Streams$Stream))
mtext(side=3,"lm, intercept only",line=2)
mtext(side=2,"Residuals",line=2)
abline(h=0)
boxplot(split(residuals(streams.lm2)/summary(streams.lm2)$sigma,Streams$Stream))
mtext(side=3,"lm, stream effect",line=2)
mtext(side=1,"Stream",line=2)
abline(h=0)
boxplot(split(residuals(streams.lme1)/summary(streams.lme1)$sigma,Streams$Stream))
mtext(side=3,"lmer, random stream effect",line=2)
abline(h=0)

# plot the fitted values
# library(lattice)
# plot(streams.lme1,Density~fitted(.)|Stream, abline = c(0,1))
lme1_results <- streams.lme1 |>
  augment() |>
  janitor::clean_names()
lme1_results |>
  ggplot() +
  aes(x = fitted, y = density) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~stream)

# alternative using lme
streams.lme1 <- lme(Density~1,random=~1|Stream,data=Streams,meth="ML")
# diagnostic plots of residuals, random effects
par(mfrow=c(2,3))
# QQ plot
qqnorm(streams.lme1$residuals/streams.lme1$sigma,ylab="Quantiles of Residuals")
qqline(streams.lme1$residuals/streams.lme1$sigma)

plot(residuals(streams.lme1)/streams.lme1$sigma,ylab="Standardized Residuals")
hist(residuals(streams.lme1)/streams.lme1$sigma,xlab="Standardized Residuals",main="")

#homogeneity of within group variance
boxplot(split(residuals(streams.lme1)/streams.lme1$sigma,Streams$Stream),ylab="Standardized Residual",xlab="Stream",csi=0.2)
abline(0,0,lwd=3)
#normality of the between-group residuals
re.sigma <- as.numeric(VarCorr(streams.lme1)[1,2])
print(streams.lme1$coefficients$random$Stream)
re<-streams.lme1$coefficients$random$Stream/re.sigma  ##/0.6184
qqnorm(re,ylab="Quantiles of random effects")
qqline(re)
hist(streams.lme1$coefficients$random$Stream/re.sigma,xlab="random effects",main="")



##### to the TMB Robin! ########

## parameters
parameters <- list(beta=0, 
                   log_sigma_o=0, 
                   log_sigma_b=0, 
                   b=rep(0,length(unique(Streams$Stream))))

## model
streams_model <- function(parms) {
  getAll(Streams, parms, warn=FALSE)
  ## Optional (enables extra RTMB features)
  Density <- OBS(Density)
  ## contribution of the random (stream) effects
  b %~% dnorm(rep(0, length(b)), sd = exp(log_sigma_b))
  ## Data
  pred_density <- beta + b[Stream]
  Density %~% dnorm(pred_density, sd = exp(log_sigma_o))
  # ## Get predicted  uncertainties
  sigma_o <- exp(log_sigma_o)
  sigma_b <- exp(log_sigma_b)
  ADREPORT(sigma_o)
  ADREPORT(sigma_b)
  ADREPORT(pred_density)
  # ## Return
  # nll
}

#process the objective function
obj <- MakeADFun(streams_model, parameters, random = c("b"))

#fit the model using nlminb
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit
#TMBAIC(fit)


# get variance estimates for model parameters
sdr <- sdreport(obj)
sdr
summary(sdr) #gives the ADREPORT variables too.

summary(sdr,select="random")
summary(sdr,select="fixed")


#consistency check
set.seed(1)
chk <- checkConsistency(obj)
chk
#chk <- checkConsistency(obj, estimate=TRUE)
#summary(chk)

# look at the residuals
#goodness of fit using One step ahead residuals
osa <- oneStepPredict(obj, 
                      method="fullGaussian",
                      discrete=FALSE)
qqnorm(osa$res); abline(0,1)





