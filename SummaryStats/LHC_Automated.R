### Script that runs the LHC analysis, takes advantage of rslurm to run ~100 different optimisations with distinct starting points.
# Runs a nested version that eliminates one of each parameters (a, b, tX, tY, U) and compares the subsequent model
# with the parent model using LRT. 
library(data.table)
library(rslurm)
library(tidyverse)
library(stringr)
library(extRemes)

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

grand_dir="~/LDL-BMI/" ### Directory specific to the trait pair to be analysed.
#Should contain the summary stat files if non-UKBB SNPs were made to overlap in this directory.
setwd(grand_dir)
comp_run=TRUE   ### Change to false if complete model was already ran.

fX = fread(paste0(grand_dir,"/LDL_uniq.tsv")) ### Exposure file if non-UKBB (overlapping SNPs between EXP-OUT)
fY = fread(paste0(grand_dir,"/BMI_uniq.tsv")) ### Outcome file if non-UKBB (overlapping SNPs between EXP-OUT)
#fX = fread(paste0("~/data/", EXP ,"_uniq.tsv")) ### Exposure preprocessed summary stats from UKBB
#fY = fread(paste0("~/data/", OUT ,"_uniq.tsv")) ### Outcome preprocessed summary stats from UKBB

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
all(fX$rsid==ld_wd$rs)	
all(fX$rsid==fY$rsid)

nX = mean(fX$n_complete_samples)  #Get sample size for trait X
nY = mean(fY$n_complete_samples)  #Get sample size for trait Y

bX = fX$tstat/sqrt(nX)   #Get standardised beta for trait X
bY = fY$tstat/sqrt(nY)   #Get standardised beta for trait Y

#Take every 10th element for faster computation
bX = bX[seq(1, length(bX), 10)]
bY = bY[seq(1, length(bY), 10)]
ld = as.numeric(unlist(ld_wd[seq(1, nrow(ld_wd), 10),'LDSC']))
wd = as.numeric(unlist(ld_wd[seq(1, nrow(ld_wd), 10),'weight']))

M = 10106833  #Number of SNPs

### First run of the complete model. To be compared against the nested models in the following section
if(comp_run){
  system(paste0("mkdir Comp"))
  setwd(paste0(getwd(),"/Comp/"))
  
  phen_corr = 0  #Standard value is 0, could be replaced by phenotypic correlation between X and Y. 0 value converges well enough
  nXY = (nX/380e3)*(nY/380e3)*380e3 #Estimated fraction of sample overlap
  rho = phen_corr*nXY/(as.numeric(nX)*as.numeric(nY)) #Calculated phenotypic correlation due to sample overlap
  
  # Create a theta matrix with the following columns to use as input for the rlsurm parameters:
  # iX - iY - pX - pU - pY - h2X - h2Y - tX - tY - a - b - rho
  theta.df = matrix( data = NA, nrow=120, ncol=12) #nrows indicates how many different set of starting points to run in rslurm
  theta.df1 = t(apply(theta.df,1,FUN = function(x){c(1,1,0.5*runif(3,0,1),0.5*runif(3,0,1),runif(3,0,1)-0.5,rho)} )) #Random starting points for the 12 parameters
  par.df = data.frame(par=I(apply(theta.df1,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm
  #Set or arguments to be used in optim function, upper and lower bounds, parscale factor. 
  args.df = rbind(lower = c(0,0,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,-1,-1,-1,(-1/max(nX,nY))*10),
                  upper = c(2,2,1,1,1,1,1,1,1,1,1,(1/max(nX,nY))*10),
                  parscale=c(1,1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-7))
  #Source the complete optimisation script
  source("/data/Scripts/optim_comp.R")
  start.time <- Sys.time()
  sjob = slurm_apply(f = run_optim, params = par.df, jobname = "Comp", nodes = 100, cpus_per_node = 2,
                     add_objects = c("bX","bY","M","nX","nY","ld","wd","args.df"),
                     slurm_options = list(partition = "X"), #note: partition X is representative of a single parition in cluster. 
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
  #Get output of minus log likelihood (mLL) and estimated parameters from rslurm in the form of a table with nrows equal to nrows(par.df)
  res_temp = get_slurm_out(sjob, outtype = 'table')
  write.csv(res_temp, "AllRes_parscale.csv", row.names = FALSE) ## Will be used to read new sp_mat (starting points in nested models)
  
} else {
  if(dir.exists("./Comp")){
    setwd("./Comp")
  }
} # comp_run = true/false

### Nested models
# Assuming we're in a folder that ran a complete model minimisation (/Comp) and has its results saved in Allres_parscale.csv
# These are the parameters to nest across
parameters = c("a","b","tX","tY","U")
parameters = sort(parameters) #Already sorted but just to keep order of parameter importance (luckily it's alphabetical). 
to_remove=c() #Vector that stores parameters that are removed after each nesting iteration.

parent_dir = getwd() #If we're starting in /Comp, then that's our parent model. 
setwd(parent_dir)

Nesting_boolean = TRUE
while(length(parameters)>0 & Nesting_boolean){
  #While(parameters !empty & Nesting=TRUE), if parameters is empty then you have tested all parameters and 
  #ended up with a non-model, if Nesting = False, all the LRT for that parent were non-significant.
  sp_mat=fread("AllRes_parscale.csv") #Read results from parent model, and use estimated parameters as starting points for faster computation
  
  args.df = rbind(lower = c(0,0,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,-1,-1,-1,(-1/max(nX,nY))*10),
                  upper = c(2,2,1,1,1,1,1,1,1,1,1,(1/max(nX,nY))*10),
                  parscale=c(1,1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-7)
                  )
  colnames(args.df) = c("iX", "iY", "pX","pU","pY","h2X","h2Y","tX","tY","a","b","rho") #naming them for later reference
  args.df=as.data.frame(args.df)
  part=c("X","Y","Z") #Partition used in cluster if there exists any. 
  
  res_list=list() #List to hold results of various paramters that are being nested
  counter=0 #Counter to track parameters
  
  for (param in parameters) {
    #If either tX or tY are already removed in an earlier nesting level, and the next parameter to remove is tY or tX, then skip and go staright to U
    if((param=="tX" & "tY" %in% to_remove) || (param=="tY" & "tX" %in% to_remove)){
      parameters = parameters[!parameters %in% c(param)]
      next
    }
    counter = counter + 1 #Current parameter position in list
    path = paste0(param,"/") #Create a subfolder inside parent model folder for each parameter to nest
    system(paste0("mkdir ", path))
    setwd(path)
    if( param =="U" & length(intersect(c("tX","tY"),to_remove))>0){ #Select script name with U and no tX/tY in it (doesn't exist)
      script = paste0("_",paste0(sort(c(param,to_remove[!to_remove %in% c('tX','tY')])),collapse=""))
    }else{ #Generate script name from previoulsy removed/nested parameters and current parameter to nest
      script = paste0("_",paste0(sort(c(param,to_remove)),collapse=""))
    }
    if(param =='U'){param=c('pU','tX','tY')} # When removing/nesting confounder U, its proportion of effective SNPs (pU)
    #tX, and tY must also be removed from arguments
    if("U" %in% to_remove){ 
      # If U is already removed, then column selection must take pU,tX,tY into account and disregard the non-existant U column
      param=c(param,'pU','tX','tY')
      args.df1 = select (args.df, -c(param,to_remove[-match("U",to_remove)]))
    }else{
      args.df1 = select (args.df, -c(param,to_remove))
    }
    
    sp_mat1 = select (sp_mat,-c(mLL,param[param %in% colnames(sp_mat)])) #Parameter in to_remove already gone, nested starting points
    par.df = data.frame(par=I(apply(sp_mat1,1,as.list)))
    
    source(paste0("/data/Scripts/optim",script,".R"))
    
    sjob = slurm_apply(f = run_optim, params = par.df, jobname = script, nodes = 60,cpus_per_node = 2,
                       add_objects = c("bX","bY","M","nX","nY","ld","wd","args.df1"),
                       slurm_options = list(partition = part[(counter%%length(part))+1]),
                       submit = TRUE)
    rm(run_optim) #To remove the old function and make sure we're sourcing a new one each time
    res_list[[counter]]=sjob
    setwd(parent_dir)
  }
  
  #Keep a loop open till jobs are completely done.
  wait_counter = 0
  while (wait_counter < length(parameters)) {
    wait_counter = 0
    for (x in 1:length(parameters)) {
      setwd(paste0("./",parameters[x])) ## must enter directory to access sjob
      if (tryCatch(str_detect(get_job_status(res_list[[x]]),"completed"), warning = function(w){FALSE},
                   error = function(e){FALSE})) {
        wait_counter = wait_counter + 1
      } else{
        wait_counter = wait_counter
      }
      setwd(parent_dir)
    }
  }
  
  print("Jobs are completed.")
  
  ### Get the current/parent min mLL, run LRT against the parent mLL and compare
  parent_min = sp_mat[which(sp_mat$mLL==min(sp_mat$mLL)),]
  parent_min1=cbind(model = "Parent", parent_min)
  f1="Pair_LRT.csv"
  write.table(parent_min1, file=f1, sep = ",", append = TRUE, row.names = FALSE)
  write.table("Comparison against Parent Model", file=f1, sep = ",", append = TRUE, row.names = FALSE)
  res_pval=list()
  for (x in 1:length(parameters)) {
    setwd(paste0("./", parameters[x]))
    res_temp = get_slurm_out(res_list[[x]], outtype = 'table')
    write.csv(res_temp, "AllRes_parscale.csv", row.names = FALSE) #Will be used to read new sp_mat if nesting continues
    res_min = res_temp[which(res_temp$mLL == min(res_temp$mLL)), ]
    df1 = length(parent_min) - length(res_min) #get degrees of freedom
    LRT = lr.test(x = res_min$mLL, y = parent_min$mLL, df = df1)
    res_pval[[x]] = LRT$p.value
    LRT1 = cbind(model = paste0("No ", parameters[x]), res_min, pval = LRT$p.value)
    setwd(parent_dir) #Write LRT comparison results in parent directory
    write.table(LRT1, file = f1, sep = ",", append = TRUE, row.names = FALSE)
  }
  
  print("LRT performed.")
  
  ### Save space, get rid of temp folders created by rslurm
  for (x in 1:length(parameters)) {
    setwd(paste0("./", parameters[x]))
    cleanup_files(res_list[[x]])
    setwd(parent_dir)
  }
  print("File cleanup done.")
  
  ## The order of the parameters (a, b, tX, tY, U) helps in giving priority to which element should be 
  ## removed should two of them have the same pvalue. The only issue might be when several have relatively close
  ## p-values to the threshold (0.05)
  ind_maxPval=which(unlist(res_pval)==max(unlist(res_pval)))[1] #Which parameter when removed/nested has the largest p-value
  
  if (res_pval[ind_maxPval] > 0.05) { #Bonferroni Correction: 0.05/total number of pairs tested x2 for bidirectionality
    #If nested model with largest p-value is non-significant in comparison to parent (LRT), it must become parent model
    setwd(paste0("./", parameters[ind_maxPval])) #Enter that nested model and turn it to a parent model
    parent_dir = getwd() #Nested model is now parent model
    to_remove = c(to_remove, parameters[ind_maxPval]) #Add parameter from the nested model to the to_remove list
    parameters = parameters[-ind_maxPval] #and remove it from the parameters list
    ## Two if sttements below are unlikely to happen, have been addressed earlier.
    if ("tX" %in% to_remove & "tY" %in% to_remove) { #If either tX or tY was previously removed, and now either tY or tX got removed, skip U
      parameters = parameters[-match("U", parameters)]
    }
    if ("U" %in% to_remove) { #If U was removed, remove tX and tY
      parameters = parameters[-match(c("tX", "tY"), parameters)]
    }
  } else{ #If nested model with largest p-value is significant in comparison to parent (LRT), no more nesting needed
    if (length(to_remove) == 0) {
      pt = "/Comp"
    } else{
      pt = paste0(sep = "/", to_remove, collapse = "")
    }
    write.table(
      paste0(grand_dir, pt),
      file = paste0(grand_dir, "/Nesting.txt"),
      append = TRUE,
      row.names = FALSE
    )
    Nesting_boolean = FALSE
    #system("Rscript /data/Scripts/MakeFigures.R") #Possible script to geenrate figures of estimated parameters of the winning model
    print("Made Figures")
    print(time.taken) #Time to run a acomplete model
  }
} ## while loop depending on Nesting_boolean and length(parameters)>0
