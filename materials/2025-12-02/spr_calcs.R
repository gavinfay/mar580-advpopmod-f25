

## read in biol data
# nages <- scan("lectures/floundah_biology.txt",n=1,skip=9)
# maturity <- scan("lectures/floundah_biology.txt",n=nages,skip=11)
# selex <- scan("lectures/floundah_biology.txt",n=nages,skip=13)
# weight <- scan("lectures/floundah_biology.txt",n=nages,skip=15)
# M <- 0.4

## Numbers-per-recruit

NPR <- function(F,M,selex) {
  nages <- length(weight)
  N <- rep(0,nages)
  N[1] <- 1
  for (age in 2:nages)
    N[age] <- N[age-1]*exp(-1*(selex[age-1]*F+M))
  N[nages] <- N[nages]/(1-exp(-1*(selex[age]*F+M)))
  return(N)
}
npr <- NPR(F=0,M,selex)

## Yield-per-recruit

YPR <- function(F,M,weight,selex) {
  N <- NPR(F,M,selex)
  Z <- selex*F+M
  yield <- sum(weight*selex*F*N*(1-exp(-Z))/Z)
  return(yield)
}
#YPR(F=0.1,npr,M,weight,selex)
YPR(F=0.1,M,weight,selex)

## Spawning Biomass per recruit

SBPR <- function(F,M,weight,selex,maturity) {
  N <- NPR(F,M,selex)
  sbpr <- sum(N*weight*maturity)
  return(sbpr)
}
SBPR(0,M,weight,selex,maturity)
SBPR(0.4,M,weight,selex,maturity)

