library(ggplot2)
library(NMOF)
library(mvtnorm)
library(stats)
library(data.table)
.libPaths(c(.libPaths(), "/data/sgg2/ninon/bin/R-3.4.3_Packages/"))
library(TwoSampleMR)

args <- commandArgs(trailingOnly = TRUE)
nX <- as.numeric(as.character(args[1]))
nY <- as.numeric(as.character(args[2]))
Zval <- as.numeric(as.character(args[3]))
RAPS <- as.character(args[4])
jobID <- as.numeric(as.character(args[5]))

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

### Read in data from previously generated data
gen_data = fread(paste0("./Data/Data_",jobID,".csv"))
M = dim(gen_data)[1] #Number of SNPs
bX = gen_data$bX  #Effects of trait X
bY = gen_data$bY  #Effects of trait X
zX = gen_data$zX  #Standardised effects of traits X and Y
zY = gen_data$zY  #Standardised effects of traits X and Y
D = cbind(bX, bY)  # used directly in LHC_SEM as a global parameter

### Repeat optimisation function 100 times with different starting points to explore likelihood space
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

### Run standard MR analysis and store both results in this array
if(RAPS){
  result_array = array(data = NA, dim = 24)
}else{
  result_array = array(data = NA, dim = 20)
}

## Prune SNPs for standard MR methods
pval=2*pnorm(-abs(Zval))
ind_prunX=unlist(prune_X(zX,pval))

## Get Standard /error for the two traits to be used in standard MR
seX = gen_data$seX
seY = gen_data$seY

## Create M fake SNP names needed for MRbase
labe = c(paste("rs", seq(M)))
ale = sample(c("A","C","T","G"), M, replace = T)

random_exp_df = data.frame(SNP = labe,
                           beta = bX,
                           se = seX,
                           effect_allele = ale)
## Select only SNPs with a Z value higher than cutoff specified in bash script
random_exp_df1 = random_exp_df[ind_prunX,]
random_exp_dat <- format_data(random_exp_df1, type="exposure")
#random_exp_dat1 <- clump_data(random_exp_dat)  ## Don't clump because the SNPs are independent. 
random_exp_dat1 <- random_exp_dat

## Get outcome data from SNPs that are also found in the exposure
random_out_df = data.frame(SNP = labe,
                           beta = bY,
                           se = seY,
                           effect_allele = ale)
random_out_dat <- format_data(random_out_df, type="outcome")
random_out_dat1=random_out_dat[random_out_dat$SNP %in% random_exp_dat1$SNP,]

## Harmonise the exposure and outcome data
dat <- harmonise_data(
  exposure_dat = random_exp_dat1, 
  outcome_dat = random_out_dat1, action = 1
)

## Perform MR
res <- mr(dat)

if(RAPS){
  ## Clean workspace
  rm(random_exp_df1,random_exp_dat,random_exp_dat1,random_out_dat,random_out_dat1,dat)
  
  ### MR RAPS filtered
  pval2 = 1e-4
  ind2_prunX=unlist(prune_X(zX,pval2))
  
  random_exp_df1 = random_exp_df[ind2_prunX,]
  random_exp_dat <- format_data(random_exp_df1, type="exposure")
  #random_exp_dat1 <- clump_data(random_exp_dat)  ## Don't clump because the SNPs are independent. 
  random_exp_dat1 <- random_exp_dat
  
  random_out_dat <- format_data(random_out_df, type="outcome")
  random_out_dat1=random_out_dat[random_out_dat$SNP %in% random_exp_dat1$SNP,]
  
  # Harmonise the exposure and outcome data
  dat <- harmonise_data(
    exposure_dat = random_exp_dat1, 
    outcome_dat = random_out_dat1, action = 1
  )
  
  # Perform MR
  res1 <- mr(dat, method_list = "mr_raps")
  
  ### MR-RAPS - unfiltered
  ## CLean workspace
  rm(random_exp_df1,random_exp_dat,random_exp_dat1,random_out_dat,random_out_dat1,dat)
  
  random_exp_dat1 <- format_data(random_exp_df, type="exposure")
  
  random_out_dat <- format_data(random_out_df, type="outcome")
  random_out_dat1=random_out_dat[random_out_dat$SNP %in% random_exp_dat1$SNP,]
  
  # Harmonise the exposure and outcome data
  dat <- harmonise_data(
    exposure_dat = random_exp_dat1, 
    outcome_dat = random_out_dat1, action = 1
  )
  
  # Perform MR
  res3 <- mr(dat, method_list = "mr_raps")
}

### Save results
result_array[1:10]=value.mat[which.min(value.mat[,1]),]
result_array[11:12]=res[1,7:8] # Egger
result_array[13:14]=res[2,7:8] # weighted median
result_array[15:16]=res[3,7:8] # IVW
result_array[17:18]=res[4,7:8] # simple mode
result_array[19:20]=res[5,7:8] # weighted mode
if(RAPS){
  result_array[21:22]=res2[1,7:8] # MR-RAPS
  result_array[23:24]=res3[1,7:8] # MR-RAPS unfiltered
}

### Write output to file
result_array=matrix(result_array,nrow = 1)
outFile2=paste(jobID,"-temp.csv", sep="")
write.table(result_array, file = outFile2, sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE)