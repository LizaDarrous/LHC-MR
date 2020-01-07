#!/usr/bin/env Rscript

### Function to create SNP summary statistics of two non overlapping samples X,Y, without the need for genomic data.
## Given M (total number of SNPs), nX (sample size 1), nY (sample size 2), h2U (heritbaility of confounder U),
# h2X (heritability of trait X), h2Y (heritability of trait Y), qX (effect of U on trait X), qY (effect of U on trait Y),
# a (causal effect of X on Y), b (causal effect of XY on X), p(X/Y/U) (proportion of effective SNPs for X/Y/U).

# In this data simulation, SNPs can be sampled from a normal or a student-t distribution.

data_SIM <- function(pX = 0.1, 
                     pU = 0.1, 
                     pY = 0.1, 
                     h2X = 0.25,
                     h2Y = 0.2,
                     h2U = 0.3,
                     qX = 0.3,
                     qY = 0.2,
                     a = 0.3,
                     b = 0.1,
                     M = 1e4,
                     nX = 1e6,
                     nY = 1e6,
                     NormDist = TRUE,
                     df = 5) {
  
  tX = qX * sqrt(h2U)
  tY = qY * sqrt(h2U)
  #tot_h2X = h2X + tX^2 + (b^2 * h2Y)
  #tot_h2Y = h2Y + tY^2 + (a^2 * h2X)
  
  if (NormDist) {
    ## normal
    piX = as.integer(runif(M, 0, 1) < pX) * (sqrt(h2X / (M * pX)) * rnorm(M))
    piU = as.integer(runif(M, 0, 1) < pU) * (sqrt(h2U / (M * pU)) * rnorm(M))
    piY = as.integer(runif(M, 0, 1) < pY) * (sqrt(h2Y / (M * pY)) * rnorm(M))
    piNx = sqrt(1 / nX) * rnorm(M) #Error term
    piNy = sqrt(1 / nY) * rnorm(M) #Error term
  } else{
    ## Student
    piX = as.integer(runif(M, 0, 1) < pX) * (sqrt(h2X / (M * pX)) * rt(M, df=df) /
                                               sqrt(df / (df - 2)))
    piU = as.integer(runif(M, 0, 1) < pU) * (sqrt(h2U / (M * pU)) * rt(M, df=df) /
                                               sqrt(df / (df - 2)))
    piY = as.integer(runif(M, 0, 1) < pY) * (sqrt(h2Y / (M * pY)) * rt(M, df=df) /
                                               sqrt(df / (df - 2)))
    piNx = sqrt(1 / nX) * rnorm(M) #Error term
    piNy = sqrt(1 / nY) * rnorm(M) #Error term
  }
  
  ######
  
  bX = piX + (qX + (b * qY)) * piU + (b * piY) + piNx
  bY = piY + (qY + (a * qX)) * piU + (a * piX) + piNy
  
  zX = sqrt(nX) * (bX) / sqrt(1 - bX ^ 2)
  zY = sqrt(nY) * (bY) / sqrt(1 - bY ^ 2)
  
  summary_list <-
    list(
      "bX" = bX,
      "bY" = bY,
      "zX" = zX,
      "zY" = zY
    )
  
  return(summary_list)
}
