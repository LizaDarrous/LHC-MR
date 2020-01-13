# LHC-MR

Latent Heritable Confounder MR (LHC-MR) is a method applicable to association summary statistics, which estimates bi-directional causal effects, direct heritability, confounder effects and sample overlap. 

LHC-MR extends the standard Mendelian Randomisation model to incorporate the presence of a latent (unmeasured) heritable confounder and estimates its contribution to the exposure and outcome traits, while simultaneously estimating the bi-directional causal effect between the two traits.

## Overview

! To be filled out !

## Usage

The source code for LHC-MR was written in R version 3.6.0. No othr language is used in computation, and thus only bash and R are needed. 
Several R-packages are needed to run the analysis, which they will be detailed below.
```

```
### Simulations
#### Prerequisites
R Packages needed to run the simualtions include:
``` 
install.packages(ggplot2);  library(ggplot2)
install.packages(NMOF); library(NMOF)
install.packages(mvtnorm);  library(mvtnorm)
install.packages(stats);  library(stats)
install.packages(TwoSampleMR);  library(TwoSampleMR)

install.packages(reshape2); library(reshape2)
install.packages(dplyr);  library(dplyr
install.packages(stringr);  library(stringr)
install.packages(ggpubr); library(ggpubr)
install.packages(extrafont);  library(extrafont)
install.packages(SimDesign);  library(SimDesign)
```

To run simulations for LHC-MR, the scripts needed are found in `LHC-MR/Simulation/`:

   * **`DataSIM.R`**
     This script generates SNP summary statistics for two non overlapping samples X and Y with the parameters specified for the simulation analysis. It is sourced and run in `optim_SIM9param.R`. 
     
   * **`optim_SIM9param.R`**
   This script first sets up the parameters for the simulation, generates the SNP summary statistics with the specified parameters and then optimises the likelihood function for 100 different starting points. For the same SNP summary statistics, the causal estimate from trait X onto Y is caluclated using 5 standard MR methods (MR-Egger, weighted median, IVW, simple mode, weighted mode). It returns the estimated parameters of the largest likelihood from the 100 different optimisations, as well the causal estimate from the standard MR method in a csv file.
   
   * **`LHC_SIMarray.sh`**
   This bash script runs `optim_SIM9param.R` as an array job 100 different times. It takes a single paramter input, the Z-value used as a threshold to prune X-significantly associated SNPs for standard MR analysis. This array job script thus results in 100 different output csv files, numbered by the jobID and ending in "-temp.csv".

#### Results
The resulting 100 csv files can be joined into a single one with the following command: `cat *-temp.csv > results.csv`.
The result file will have the follwoing columns (mLL - pX - pU - pY - h2X - h2Y - tX - tY - a - b - EggerEst - EggerErr - WMedianEst - WMedianErr - IVWEst - IVWErr - ModeEst - ModeErr - WModeEst - WModeErr) comprised of maximum likelihood for that specific run, the corresponding parameter estimates from LHC-MR, the standard MR estimates and their corresponding standard error.

   * **`boxplots.R`**
   The resulting encompassing single file is used as input in this script to generate boxplots of the estimated parameters and calculate data summary of the estimations (bias, variance, RMSE). 


### Association Summary Statistics
#### Prerequisites
R Packages needed to run the summary stat analysis include:
``` 
library(data.table);
library(dplyr);
library(DescTools);
library(rslurm);
library(tidyverse);
library(stringr);
library(extRemes);
 
library(GGally);
library(mixtools);

```


#### Results


## Citation

Manuscript in preparation 

## Authors

Liza Darrous - liza.darrous@unil.ch

Ninon Mounier

Zoltán Kutalik

## Acknowledgments

SGG Lab members

Package creators for those used
