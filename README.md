# LHC-MR

:warning: **LHC-MR is undergoing MAJOR changes and developments. A new implementation can be found [HERE](https://github.com/LizaDarrous/lhcMR). This version is no longer maintained. Thank you for understanding!** :warning:

Latent Heritable Confounder MR (LHC-MR) is a method applicable to association summary statistics, which estimates bi-directional causal effects, direct heritability, confounder effects and sample overlap. 

LHC-MR extends the standard Mendelian Randomisation model to incorporate the presence of a latent (unmeasured) heritable confounder and estimates its contribution to the exposure and outcome traits, while simultaneously estimating the bi-directional causal effect between the two traits.


## Usage

The source code for LHC-MR was written in R version 3.6.0. No other language is used in the computation, and thus only bash and R are needed. 
Several R-packages are needed to run the analysis, which will be detailed below.


### Simulations
#### Prerequisites
R Packages needed to run the simulations include:
``` 
install.packages(ggplot2);  library(ggplot2)
install.packages(NMOF); library(NMOF)
install.packages(mvtnorm);  library(mvtnorm)
install.packages(TwoSampleMR);  library(TwoSampleMR)

install.packages(reshape2); library(reshape2)
install.packages(dplyr);  library(dplyr)
install.packages(stringr);  library(stringr)
install.packages(ggpubr); library(ggpubr)
install.packages(extrafont);  library(extrafont)
install.packages(SimDesign);  library(SimDesign)
```

To run simulations for LHC-MR, the scripts needed are found in `LHC-MR/Simulation/`:

   * **`DataSIM.R`**
     This script generates SNP summary statistics for two non-overlapping samples X and Y with the parameters specified for the simulation analysis. It is sourced and run in `optim_SIM9param.R`. 
     
   * **`optim_SIM9param.R`**
   This script first sets up the parameters for the simulation, generates the SNP summary statistics with the specified parameters and then optimises the likelihood function for 100 different starting points. For the same SNP summary statistics, the causal estimate from trait X onto Y is calculated using 5 standard MR methods (MR-Egger, weighted median, IVW, simple mode, weighted mode). It returns the estimated parameters of the largest likelihood from the 100 different optimisations, as well the causal estimate from the standard MR method in a csv file.
   
   * **`LHC_SIMarray.sh`**
   This bash script runs `optim_SIM9param.R` as an array job 100 different times. It takes a single parameter input, the Z-value used as a threshold to prune X-significantly associated SNPs for standard MR analysis. This array job script thus results in 100 different output csv files, numbered by the jobID and ending in "-temp.csv".

#### Results
The resulting 100 csv files can be joined into a single one with the following command: `cat *-temp.csv > results.csv`.
The result file will have the following columns (mLL - pX - pU - pY - h2X - h2Y - tX - tY - a - b - EggerEst - EggerErr - WMedianEst - WMedianErr - IVWEst - IVWErr - ModeEst - ModeErr - WModeEst - WModeErr) comprised of maximum likelihood for that specific run, the corresponding parameter estimates from LHC-MR, the standard MR estimates and their corresponding standard error.

   * **`boxplots.R`**
   The resulting encompassing single file is used as input in this script to generate boxplots of the estimated parameters and calculate data summary of the estimations (bias, variance, RMSE). 


### Association Summary Statistics
#### Prerequisites
R Packages needed to run the summary stat analysis include:
``` 
install.packages(data.table);  library(data.table)
install.packages(data.table);  library(data.table)
install.packages(DescTools); library(DescTools)
install.packages(rslurm);  library(rslurm)
install.packages(tidyverse); library(tidyverse)
install.packages(stringr); library(stringr)
install.packages(extRemes);  library(extRemes)
 
install.packages(GGally);  library(GGally)
install.packages(mixtools);  library(mixtools)
```
A very important package used in the analysis is [rslurm](https://cran.r-project.org/web/packages/rslurm/rslurm.pdf), which allows us to submit array jobs from within R without having to create a bash script to do so. Moreover, this parallelisation step takes advantage of the presence of a cluster with several partitions, onto which it can simultaneously distribute and run array jobs. Once rslurm is installed, it's important to edit the template files that come with the package to reflect the info of the cluster being used.

The LHC-MR association summary statistics analysis for multiple traits was written with a certain folder hierarchy assummed, such as the following (note that this is just a suggestion):
```
~ --|---- data
    |       |---- Scripts (optim.R and variations)
    |       |--- LDscores.txt
    |       |--- Unprocessed summary stat files
    |
    |---- DataPreprocess
    |       |--- get_tstatXY_[GWAS].R scripts
    |
    |---- pair1 (trait1-trait2)
    |---- pair2 (trait1-trait3)
    |---- pair3 (trait2-trait3)
    |---- so on...
```
Raw summary stat files downloaded in `data` require certain columns to continue with analysis: rsid, EffectAllele/A1/alt, OtherAllele/A2/ref, number of samples (effective sample size), beta/effect size of effective allele,
se/standard error of effect size.

These files ought to be processed using the `get_tstatXY_[GWAS].R` scripts in order to calculate the t-statistics (beta/se) of the SNPs, which will be used in the main LHC-MR analysis.

The script `get_tstatXY_UKBB.R` is slightly different than the rest. It first selects the set of SNPs for which imputation quality info>0.99 and MAF>0.005 using the `LDscore.txt` file. 
For LD scores you can either use the classical ones (https://data.broadinstitute.org/alkesgroup/LDSCORE/) or we have also created LD scores based on UK10K sequencing data (https://drive.google.com/file/d/1uua6zIournPvcJAs-QVWs8pTUVLtZ5Ym/view?usp=sharing), where the list of SNPs are optimised for UKB-based summary stats (e.g. by Neale).
Then for these set of SNPs, it uses the `variants.tsv.bgz` obtained from Neale's [UKBB GWAS Imputed v3](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=178908679) to merge variant information to variant details of the summary statistics files. This leads to the addition of the A1 and A2 columns as well as the rsid, all of which are needed for further analysis. Some duplication occurs when `variants.tsv.bgz` is merged with UKBB summary stat files, which is handled in the scripts.
This script is used once per UKBB trait, allowing all the traits to then have the same set of overlapping SNPs when running trait pair analysis between UKBB. 

The other scripts in this folder contain the processing procedure for trait pairs that are not both UKBB. In this case, an overlapping set of SNPs ought to be found for each trait pair to be analysed, and this is where the folder structure presented above comes in handy. Each trait pair's summary stat data is read from the original data folder, then a set of overlapping SNPs is found, columns are renamed to include A1, A2 and rsid. T-stat and effective sample size is calculated for the non-UKBB trait, and the strands of the two traits are harmonised. 
This smaller set of overlapping and harmonised SNPs is written into the pair folder (trait1-trait2) to be used in the LHC-MR analysis. 

For the LHC-MR analysis, the script `LHC_Automated.R` first reads in the summary stat files for the exposure (X) and outcome (Y). It then obtains LD scores and weights for the set of overlapping SNPs, calculates the standardised effects of the SNPs and runs the complete model in a new directory called `Comp` in order to estimate the 12 unknown parameters:-
- iX, iY: population structure.
- pX, pU, pY: proportion of SNPs with an effect on X, confounder U, and Y.
- h2X, h2Y: direct heritabilities of X and Y.
- tX, tY: confounder effect of X and Y.
- a, b: causal effect of X->Y and Y->X.
- rho: phenotypic correlation due to sample overlap.

Once the complete model is ran, nested models are subsequently run, each time removing one parameter from the following {a,b,tX,tY,U(tX,tY,pU)} and running the analysis in created directories named after the removed parameter. The parent model (in this case the complete model) is then tested using LRT against each of the nested models. If the estimated parameters of the nested model explained the fit better than the parent model, then that nested model becomes the parent model and the remaining parameters are each nested in turn and tested against the new parent model using LRT. This is repeated until a model that is better than all its nested models is found or all 5 parameters are exhausted.

The folder hierarchy within a pair-trait folder is as follows:
```
---- pair1 (trait1-trait2)
      |--- Nesting.txt (reveals the nesting levels) 
      |--- trait1_uniq.tsv
      |--- trait2_uniq.tsv
      |---- Comp
             |--- Pair_LRT.csv (results of the LRT between the parent/complete model and all the nested models)
             |--- AllRes_Parscale.csv (results of the parameter estimations from the 100 random starting points)
             |---- b
             |---- U
             |---- tX
             |---- tY
             |---- a
                   |--- Pair_LRT.csv (results of the LRT between the parent model and all the new nested models)
                   |--- AllRes_Parscale.csv (results of the parameter estimations from the 100 random starting points)
                   |--- Pairs.pdf (plots from MakeFigures.R found in /data/Scripts)
                   |--- Hist.pdf (plots from MakeFigures.R found in /data/Scripts)
                   |---- b
                   |---- U
                   |---- tX
                   |---- tY
```
#### Results
The results of the LRT performed between each parent model and its nested models is found in the directory of the parent model under the name `Pair_LRT.csv`. In that file, the p-value of the LRT between the parent model and the nested model that removes a single parameter is reported and is used as a measure for the significance of that parameter.

The results of the optimisation run in rslurm with 100 different starting points is found in `AllRes_parscale.csv`. These parameter estimates are used as future starting points for nested models. Furthermore, the smallest negative likelihood (largest likelihood) is reported back for that (nested) model as well as its estimated parameters.

Finally, when a final model is reached, the loop of nesting ends, and the script `MakeFigures.R` is run on the final `AllRes_parscale.csv` to obtain histogram plots for all the remaining estimated parameters as well as pairwise plots between them that can reveal the presence of patterns (i.e. bimodal causal effect).

## Citation

Preprint: https://www.medrxiv.org/content/10.1101/2020.01.27.20018929v2.article-metrics 

## Authors

Liza Darrous - liza.darrous@unil.ch

Ninon Mounier

Zoltán Kutalik

## Acknowledgments

SGG Lab members

Package creators for those used
