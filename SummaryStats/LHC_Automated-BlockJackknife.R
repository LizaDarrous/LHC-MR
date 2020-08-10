### Script that runs the LHC analysis, takes advantage of rslurm to run ~30 different optimisations with distinct starting points.
# Obtains SE of estimated parameters using block jackknife. Out of ~200 blocks, each block is left out, and the optimisation
# is run for the same set of 30 random starting points
library(data.table)
library(rslurm)
library(tidyverse)
library(stringr)
library(extRemes)
library(MASS)


args <- commandArgs(trailingOnly = TRUE)
Gdir <- as.character(args[1])
EXP <- as.character(args[2])
OUT <- as.character(args[3])

### updated function from rslurm that checks if slurm jobs are still running on cluster
get_job_status <- function (slr_job) 
{
  if (!(class(slr_job) == "slurm_job")) 
    stop("input must be a slurm_job")
  stat <- suppressWarnings(system(paste("squeue -n", slr_job$jobname), 
                                  intern = TRUE))
  if (length(stat) > 1) {
    res = "Job running or in queue."
  }
  else {
    res = "Job completed or stopped."
    # tmpdir <- paste0("rslurm", slr_job$jobname)
    #  out_files <- file.path(tmpdir, paste0("slurm_", 0:(slr_job$nodes - 
    #                                                      1), ".out"))
    #for (outf in out_files) {
    #  cat(paste("\n----", outf, "----\n\n"))
    #  cat(paste(readLines(outf), collapse = "\n"))
  }
  return(res)
}

setwd(Gdir)
grand_dir=getwd()

fX = fread(paste0(EXP,"_uniq.tsv")) ### Exposure preprocessed summary stats from UKBB
fY = fread(paste0(OUT,"_uniq.tsv")) ### Outcome preprocessed summary stats from UKBB

system(paste0("mkdir filtered_JK_30sp"))
setwd(paste0(getwd(),"/filtered_JK_30sp/"))
parent_dir = getwd()

# for jura + UKBB, get the two traits, and mkdir with those two traits + log file
sink(paste0(EXP,"-",OUT,"_log.txt"), append=FALSE, split=TRUE)

#LD scores and regression weights found in this file were computed based on 3781 individuals from the UK10K study.
#For the LD score calculation, squared correlations were computed between the target SNP and all sequence variants within its 2 Mb neighbourhood 
#with MAF>=0.5 (in the UK10K). For the regression weight calculation, we restricted the neighbouring variants to those part of the 4'650'107 well-imputed SNPs.
LDfile = fread("~/data/LDscores.txt", sep="\t", header=TRUE)
## Extract SNPs high-quality SNPs, defined as being present in both UK10K and UK Biobank, having MAF>1 in both data sets, 
#info>0.99 in the UK Biobank, non-significant (P_{diff}>0.05) allele frequency difference between UK Biobank and UK10K and 
#residing outside the HLA region (chr6:28.5-33.5Mb)
attach(LDfile)
mafT = .005
selF = which(info>.99 & abs(mafUK10K-mafUKBB)<(2*sqrt(1/4e3+1/360e3)) & mafUKBB>mafT & mafUK10K>mafT & 
               !(chr==6 & pos>=28.5e6 & pos<=33.5e6))
LDfile = LDfile[selF,]

### Since the summary statistics of the two traits are already filtered out for overlapping SNPs (in get_tstatXY_*.R), we can use either
#as reference for LD data.
snp_tokeep = intersect(LDfile$rs, fX$rsid)
ld_wd = LDfile[LDfile$rs %in% snp_tokeep,]
fX = fX[fX$rsid %in% snp_tokeep,]
fY = fY[fY$rsid %in% snp_tokeep,]

ld_wd=ld_wd[order(ld_wd$rs),] 
fX=fX[order(fX$rsid),] 
fY=fY[order(fY$rsid),] #All 3 ought to be the same size
print("traitX rsid match ld_wd rsid")
all(fX$rsid==ld_wd$rs)
print("traitX rsid match traitY rsid")
all(fX$rsid==fY$rsid)

fX = separate(data = fX, col = variant, into = c("chr", "pos","A1","A2"), sep = ":")
fY = separate(data = fY, col = variant, into = c("chr", "pos","A1","A2"), sep = ":")
fX = transform(fX, chr = as.numeric(chr), pos = as.numeric(pos))
fX$uniqPos = (fX$chr*5e8)+(fX$pos) #Get unique position
fY = transform(fY, chr = as.numeric(chr), pos = as.numeric(pos))
fY$uniqPos = (fY$chr*5e8)+(fY$pos) #Get unique position

ld_wd$uniqPos = (ld_wd$chr*5e8)+(ld_wd$pos) #Get unique position

nX = mean(fX$n_complete_samples)  #Get sample size for trait X
nY = mean(fY$n_complete_samples)  #Get sample size for trait Y

bX = fX$tstat/sqrt(nX)   #Get standardised beta for trait X
bY = fY$tstat/sqrt(nY)   #Get standardised beta for trait Y

#Take every 10th element for faster computation
bX = bX[seq(1, length(bX), 10)]
bY = bY[seq(1, length(bY), 10)]
chr = fX1[seq(1, nrow(fX1), 10),'chr']
pos = fX1[seq(1, nrow(fX1), 10),'pos']
uniqPos = fX1[seq(1, nrow(fX1), 10),'uniqPos']
ld = as.numeric(unlist(ld_wd[seq(1, nrow(ld_wd), 10),'LDSC']))
wd = as.numeric(unlist(ld_wd[seq(1, nrow(ld_wd), 10),'weight']))

big_df= cbind(bX, bY, ld, wd, chr, pos, uniqPos)
### Divide the snps into 202 blocks with 11 snps leftover added to the last chunk
nBlock = 200
nSNP = nrow(big_df) %/% nBlock  #2343
limit = nBlock * nSNP
leftover = as.numeric(nrow(big_df) %% nSNP )
start_ind = which(big_df$uniqPos %in% unlist(big_df[seq(1, limit, nSNP),'uniqPos']))
end_ind = start_ind+(nSNP-1)
end_ind[length(end_ind)]=end_ind[length(end_ind)]+leftover

JK_index=cbind(start_ind,end_ind) # 200 rows/bins

M = 10106833  #Number of SNPs

### First run of the complete model. To be compared against the nested models in the following section

#phen_corr = 0  #Standard value is 0, could be replaced by phenotypic correlation between X and Y. 0 value converges well enough
#nXY = (nX/380e3)*(nY/380e3)*380e3 #Estimated fraction of sample overlap
#rho = phen_corr*nXY/(as.numeric(nX)*as.numeric(nY)) #Calculated phenotypic correlation due to sample overlap

# Create a theta matrix with the following columns to use as input for the rlsurm parameters:
# iX - iY - pX - pU - pY - h2X - h2Y - tX - tY - a - b - rho
theta.df = matrix( data = NA, nrow=30, ncol=12) #nrows indicates how many different set of starting points to run in rslurm
theta.df1 = t(apply(theta.df,1,FUN = function(x){c(1,1,0.5*runif(3,0,1),0.5*runif(3,0,1),runif(3,0,1)-0.5,0)} )) #Random starting points for the 12 parameters
par.df = data.frame(par=I(apply(theta.df1,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm
#Set of arguments to be used in optim function, upper and lower bounds, parscale factor. 
args.df = rbind(lower = c(0,0,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,-1,-1,-1,(-1/max(nX,nY))*10),
                upper = c(2,2,1,1,1,1,1,1,1,1,1,(1/max(nX,nY))*10),
                parscale=c(1,1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-7))
colnames(args.df) = c("iX", "iY", "pX","pU","pY","h2X","h2Y","tX","tY","a","b","rho") #naming them for later reference
args.df=as.data.frame(args.df)
par.df2 = merge(par.df,JK_index)  #202 rows, for each bin same sp 
#Source the complete optimisation script
source("/data/sgg2/liza/SEM_Real/auto/Scripts/optim_comp.R")

start.time <- Sys.time()
sjob = slurm_apply(f = run_optim, params = par.df2, jobname = "JKcomp", nodes = 400, cpus_per_node = 1,
                   add_objects = c("bX","bY","M","nX","nY","ld","wd","args.df"),
                   slurm_options = list(partition = "X", time="2-00:00:00", `cpus-per-task`=1, share=TRUE, mem="3G"), #note: partition X is representative of a single parition in cluster.
                   submit = TRUE)

rm(run_optim)
#Keep a loop open till the job is completely done.
wait_counter = 0
while (wait_counter < 1) {
  wait_counter = 0
  if (tryCatch(str_detect(get_job_status(sjob),"completed"), warning = function(w){FALSE},
               error = function(e){FALSE})) {
    wait_counter = wait_counter + 1
  } else{
    wait_counter = wait_counter
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste0("Time taken:",time.taken))

#Get output of minus log likelihood (mLL) and estimated parameters from rslurm in the form of a table with nrows equal to nrows(par.df)
#setwd(parent_dir)
res_temp = get_slurm_out(sjob, outtype = 'table')
write.csv(res_temp, "AllRes_parscale.csv", row.names = FALSE)

res_min = res_temp %>%
  group_by(start_ind) %>%
  slice(which.min(mLL)) %>%
  ungroup()

write.csv(res_min, "MinRes_parscale.csv", row.names = FALSE) 

### get bimodal values, SE and mean
res_minFil = dplyr::select(as.data.frame(res_min), -c(mLL,start_ind,end_ind))
JK_res = as.data.frame(matrix(data=NA, nrow=ncol(res_minFil),ncol=17))
colnames(JK_res)=c("Parameter","mean","median","se","se_JK","var","loglik_1comp","AIC_1comp","convergance_2comp","mu1","mu2","sigma1","sigma2","lambda1","lambda2","loglik_2comp","AIC_2comp")
JK_res$Parameter=colnames(res_minFil)
JK_res$mean = colMeans(res_minFil)
JK_res$median = apply(res_minFil, 2, median)
JK_res$se = apply(res_minFil, 2, sd)
JK_res$se_JK = (apply(res_minFil, 2, sd))*sqrt(nBlock-1)
JK_res$var = apply(res_minFil, 2, var)

for (x in 1:ncol(res_minFil)){
  param = res_minFil[,x]
  Xf1 = MASS::fitdistr(param, "normal")
  Xf2 = capture.output(mixtools::normalmixEM(param, k=2, maxit=100000))
  AIC1 = 2*2 - 2*(Xf1$loglik)
  AIC2 = 2*4 - 2*(as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "loglik")[1]+1]), " " )[[1]][2]))
  JK_res$loglik_1comp[x] = Xf1$loglik
  JK_res$loglik_2comp[x] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "loglik")[1]+1]), " " )[[1]][2])
  JK_res$AIC_1comp[x] = AIC1
  JK_res$AIC_2comp[x] = AIC2
 
  print(colnames(res_minFil)[x])
  JK_res$convergance_2comp[x] = !any(str_detect(Xf2, "WARNING! NOT CONVERGENT!"))
  JK_res[x,10:11] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "mu")+1]), " " )[[1]][2:3])
  JK_res[x,12:13] = (as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "sigma")+1]), " " )[[1]][2:3]))
  JK_res[x,14:15] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "lambda")+1]), " " )[[1]][2:3])
}

tstat = (JK_res$mu1 - JK_res$mu2) / sqrt(JK_res$sigma1^2 + JK_res$sigma2^2)
t.pval = 2*pnorm(-abs(tstat))

JK_res$tstat = tstat
JK_res$tstat_pval = t.pval
JK_res$sigma1_corr = JK_res$sigma1*sqrt(nBlock-1)*sqrt(10)
JK_res$sigma2_corr = JK_res$sigma2*sqrt(nBlock-1)*sqrt(10)

JK_res$ci_lower = JK_res$mean - (1.96*JK_res$se_JK)
JK_res$ci_upper = JK_res$mean + (1.96*JK_res$se_JK)

bimo = which(JK_res$tstat_pval == 0)

JK_res$bimod = "FALSE"
JK_res$bimod[bimo] = "TRUE"

for(x in bimo){
  if(JK_res$lambda1[x]>JK_res$lambda2[x]){
    JK_res$ci_lower[x] = JK_res$mu1[x] - (1.96*JK_res$sigma1_corr[x])
    JK_res$ci_upper[x] = JK_res$mu1[x] + (1.96*JK_res$sigma1_corr[x])
  }else{
    JK_res$ci_lower[x] = JK_res$mu2[x] - (1.96*JK_res$sigma2_corr[x])
    JK_res$ci_upper[x] = JK_res$mu2[x] + (1.96*JK_res$sigma2_corr[x])
  }
}

write.csv(JK_res,paste0("JKres_",EXP,"-",OUT,".csv"))

cov_matrix = cov(res_minFil)
write.csv(cov_matrix,paste0("VarCovMatrix_",EXP,"-",OUT,".csv"))
print("Done!")
sink()
