nages <- scan("data/floundah_biology.txt",n=1,skip=9)
maturity <- scan("data/floundah_biology.txt",n=nages,skip=11)
selex <- scan("data/floundah_biology.txt",n=nages,skip=13)
weight <- scan("data/floundah_biology.txt",n=nages,skip=15)

M = 0.4

YPR <- function(F = 0, M = 0.4, selex, weight) {
  nages <- length(weight)
  N <- rep(1, nages)
  Z <- F*selex+M
  for (i in 2:nages)
    N[i] <- N[i-1]*exp(-Z[i-1])
  N[nages] <- N[nages]/(1-exp(-Z[i]))

  Y = sum(weight * N * selex * F * (1 - exp(-Z)) / Z)  
  return(Y) 
}
#example call
YPR(F = 1, M = 0.4, selex = selex, weight = weight)

SBPR <- function(F = 0, M = 0.4, selex, weight, maturity) {
  nages <- length(weight)
  N <- rep(1, nages)
  Z <- F*selex+M
  for (i in 2:nages)
    N[i] <- N[i-1]*exp(-Z[i-1])
  N[nages] <- N[nages]/(1-exp(-Z[i]))
  
  SBPR <- sum(weight * maturity * N)
  return(SBPR)  
}
SBPR(F = 0, M = 0.4, selex = selex, weight = weight, maturity = maturity)
SBPR(F = 0.5, M = 0.4, selex = selex, weight = weight, maturity = maturity)






