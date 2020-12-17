library(ggplot2)
library(NMOF)
library(mvtnorm)
library(stats)
library(TwoSampleMR)

args <- commandArgs(trailingOnly = TRUE)
Zval <- as.numeric(as.character(args[1]))
jobID <- as.numeric(as.character(args[2]))

### Structural Equation Model, calculates log likelihood from parameters
LHC_SEM = function(x) {
  pX = x[1]
  pU = x[2]
  pY = x[3]
  h2X = x[4]
  h2Y = x[5]
  tX = x[6]
  tY = x[7]
  a = x[8]
  b = x[9]
  
  L = NA
  
  S0 = matrix(c(1/nX, 0, 0, 1/nY), ncol = 2, nrow = 2)
  S1 = (h2X / (M*pX)) * matrix(c(1, a, a, a^2), ncol = 2, nrow = 2)
  S2 = (1 / (M*pU)) * matrix(c((tX + (b*tY))^2, (tX + (b*tY)) *
       (tY + (a*tX)), (tX + (b*tY)) * (tY + (a*tX)), (tY + (a*tX))^2), ncol = 2, nrow = 2)
  S3 = (h2Y / (M*pY)) * matrix(c(b^2, b, b, 1), ncol = 2, nrow = 2)
  
  mu = c(0,0)
  
  L = (pX*(1-pU)*(1-pY)*dmvnorm(D, mu, S1+S0))+
      (pU*(1-pX)*(1-pY)*dmvnorm(D, mu, S2+S0))+
      (pY*(1-pX)*(1-pU)*dmvnorm(D, mu, S3+S0))+
      (pX*pU*(1-pY)*dmvnorm(D, mu, S1+S2+S0))+
      (pX*pY*(1-pU)*dmvnorm(D, mu, S1+S3+S0))+
      (pY*pU*(1-pX)*dmvnorm(D, mu, S3+S2+S0))+
      (pY*pU*pX*dmvnorm(D, mu, S1+S2+S3+S0))+
      (1-pX)*(1-pU)*(1-pY)*dmvnorm(D, mu, S0)
  L[which(L==0)]=10e-100
  L2 = -sum(log(L))
  return(L2)
}

### Pruning function to select instruments to be used in MR
prune_X = function(zX,p_limit=1e-5){
  zX=zX
  z_limit=abs(qnorm(0.5*p_limit))
  ind_keep=which(abs(zX)>z_limit)
  ind_keep=unique(ind_keep)
  ind_keep=list(ind_keep)
  return(ind_keep)
}


source("DataSIM.R") # data genetration functions (takes following parameters)
pX = 0.1; pU = 0.05; pY = 0.15
h2X = 0.25; h2U = 0.3; h2Y = 0.2
qX = 0.3; qY = 0.2; a = 0.3; b = 0
M = 5e4; nX = 5e5; nY = 5e5
NormDist = TRUE; kurt = 5

D_list = data_SIM(
  pX = pX, pU = pU, pY = pY,
  h2X = h2X, h2Y = h2Y, h2U = h2U,
  qX = qX, qY = qY, a = a, b = b,
  M = M, nX = nX, nY = nY,
  NormDist = TRUE
)

bX = D_list$bX  #Effects of trait X
bY = D_list$bY  #Effects of trait X
zX=D_list$zX; zY=D_list$zY  #Standardised effects of traits X and Y
D = cbind(bX, bY)  # used directly in LHC_SEM as a global parameter

## Repeat optimisation function 100 times with different starting points to explore likelihood space
value.mat = matrix(data = NA, nrow = 100, ncol = 10)
for (i in 1:100) {
  para=c(runif(5,0,0.5),runif(4,-1,1)) # Generate random starting parameters
  test = optim(para, LHC_SEM, method = "L-BFGS-B", lower = c(0,0,0,0,0,0,-1,-1,-1), upper = c(1,1,1,1,1,1,1,1,1),
               control = list(maxit = 5e3))
  test_min_value=test$value
  value.mat[i,1]=test_min_value
  test_min_array=test$par
  value.mat[i,2:10]=test_min_array
}

## Run standard MR analysis and store both results in this array
result_array = array(data = NA, dim = 20)

## Prune SNPs for standard MR methods
pval=2*pnorm(-abs(Zval))
ind_prunX=unlist(prune_X(zX,pval))

## Get Standard /error for the twi traits to be used in standard MR
SE.x=sqrt(1-bX^2)/sqrt(nX)
SE.y=sqrt(1-bY^2)/sqrt(nY)

### Create M fake SNP names needed for MRbase
labe=c(paste("rs",seq(M)))
ale=rep(c("A","C","T","G"),M)
ale=sample(ale)

random_exp_df = data.frame(SNP = labe,
                           beta = bX,
                           se = SE.x,
                           effect_allele = ale)
### Select only SNPs with a Z value higher than cutoff specified in bash script
random_exp_df = random_exp_df[ind_prunX,]

random_exp_dat <- format_data(random_exp_df, type="exposure")
#random_exp_dat2 <- clump_data(random_exp_dat)
random_exp_dat2 <- random_exp_dat

random_out_df = data.frame(SNP = labe,
                           beta = bY,
                           se = SE.y,
                           effect_allele = ale)

random_out_dat <- format_data(random_out_df, type="outcome")
random_out_dat2=random_out_dat[random_out_dat$SNP %in% random_exp_dat2$SNP,]

### Harmonise the exposure and outcome data
dat <- harmonise_data(
  exposure_dat = random_exp_dat2, 
  outcome_dat = random_out_dat2, action = 1
)

### Perform MR
res <- mr(dat)

result_array[1:10]=value.mat[which.min(value.mat[,1]),]
result_array[11:12]=res[1,7:8] # Egger
result_array[13:14]=res[2,7:8] # weighted median
result_array[15:16]=res[3,7:8] # IVW
result_array[17:18]=res[4,7:8] # simple mode
result_array[19:20]=res[5,7:8] # weighted mode

## Write output to file
result_array=matrix(result_array,nrow = 1)
outFile2=paste(jobID,"-temp.csv", sep="")
write.table(result_array, file = outFile2, sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE)
