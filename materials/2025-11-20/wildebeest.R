#wildebeest
set.seed(8675309)

yr=c(61,63,65,67,71,72,77,78,82,84,86,89)
Nhat=c(263,357,439,483,693,773,1444,1249,1209,1338,1146,1686)/1000
se=c(NA,NA,NA,NA,28.8,76.7,200,355,272,138,133,176)/1000
cvhat=se/Nhat
names(Nhat)=names(se)=names(cvhat)=yr
# randomly select from exisitng CVs to create CVs fsro first four years
nas=which(is.na(cvhat))
cvhat[nas]=sample(na.omit(cvhat),length(nas),replace=TRUE)
sehat=cvhat*Nhat
allys=60:89
ymin=min(allys)
y=yr-ymin
ys=allys-ymin

# Rainfall:
rains=c(100,38,100,104,167,167,165,79,91,77,134,192,235,159,211,257,204,300,187,84,99,163,97,228,208,83,
        44,112,191,202)/100
names(rains)=allys
rain=c(38,104,167,79,192,235,300,187,97,208,44,202)/100 # just in years with N estimates
names(rain)=yr

# Catch:
Ca=c(rep(0,length(allys)))
Ca[allys>=77]=40/1000 # Hilborn & Mangel's best estimate of poaching: 40,000 poached all years from 1977
names(Ca)=allys

# Calf survival estimates:
phi1hat=c(.5,.25,.3,.26,.36,.32,.35,.36) # calf survival estimates
yrc=c(64,65,66,67,68,70,71,72) # calf survival years
yc=yrc-ymin # shift origin to same as for other years

# Adult survival estimates:
phi2hat=c(1-0.017,1-0.014,1-0.008,1-0.005,1-0.027,1-0.021)^4 # adult survival estimates
yra=c(68,69,71,72,82,83) # adult survival years
ya=yra-ymin # shift origin to same as for other years

lower=Nhat-1.96*sehat
upper=Nhat+1.96*sehat