#!/bin/bash

#SBATCH --job-name= 								# Job name
#SBATCH --workdir= 									# The Working Directory of the job
#SBATCH --mail-type=ALL								# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user= 								# Where to send mail, email address
#SBATCH --ntasks=2									# Run 32 tasks (one per core)
#SBATCH --time=24:00:00								# Time limit hrs:min:sec
#SBATCH --output=slurm-array-%A-%a-%j.out 			# Standard output and error log
### %j: job allocation number, %A: SLURM_ARRAY_JOB_ID, %a: SLURM_ARRAY_TASK_ID (as defined above)
#SBATCH --account=									# runs on the X nodes. 
#SBATCH --partition=								# runs on the X node partition. 
#SBATCH --array=1-100%30							# Runs arrays jobs from [a] to [b], with a limit of X at a time - [a]-[a]%X

Zval=4												# Z-value used to prune SNPs, above which are used as instrumental values in MR 

### Data generation for simulation using R, parameters changed inside this script.
Rscript optim_SIM9param.R ${Zval} ${SLURM_ARRAY_TASK_ID}


## run with terminal command: sbatch LHC_SIMarray.sh
## This is an sbatch script that will run the Rscript (optim_SIM9param.R) 100 different times (100 different data generations), outputing the results into files whose names are
## specified in the Rscript with a number attached that corresponds to the SLURM_ARRAY_TASK_ID (runs from 1-100) to separate the different resulting files.

## The 100 different files can then be joined into one csv when this scipt is done using cat *-temp.csv > final_simulation.csv and analysed
## Header for the resulting files is: mLL - pX - pU - pY - h2X - h2Y - tX - tY - a - b - EggerEst - EggerErr - WMedianEst - WMedianErr - IVWEst - IVWErr - ModeEst - ModeErr - WModeEst - WModeErr
## mLL : maximum log likelihood, this is the smallest value (maximum likelihood) out of the 100 different starting points used for the optimistaion of a single data geenration
## pX/pU/pY : estimated parameters corresponding to mLL, represent the proportion (pi) of effective SNPs on trait X, confounder U, and trait Y
## h2X - h2Y : estimated parameters corresponding to mLL, represent the direct heritability on trait X and trait Y respectively 
## tX - tY : estimated parameters corresponding to mLL, represent the confounder effect on trait X and trait Y respectively 
## a - b : estimated parameters corresponding to mLL, represtn the causal effect from trait X onto trait Y and the reverse causal effect from trait Y on to trait X
## EggerEst - EggerErr, etc... : TwoSampleMR causal estimation of trait X onto trait Y and the standard error reported by the function. Methods selected are MR Egger, Weighted Median, 
## Inverse Variance Weighted, Simple Mode, Weighted Mode