### Generating N different data sets to be used in LHC-MR and standard MR methods
# This is to avoid having to repeat the data generation process each time a new method is to be compared
# A seed is set to generate the different data, they are saved in their respective folders (named after the seeting) 
# under a subfloder called Data
# Need to source DataSim, set up the parameters for that specific setting and generate data N times (100 usually)
library(NMOF)
library(mvtnorm)
library(stats)

## Set wd as the folder of the setting to be tested
setwd("/data/sgg3/liza/SEM_Sim/Unified_Sim/Standard")
## Set number of data to be generated
nD = 100
## Make a directory to store the generated data in 
system(paste0("mkdir Data"))
setwd(paste0(getwd(),"/Data/"))
## Settings used to geenrate data
pX = 0.1
pU = 0.05
pY = 0.15
h2X = 0.25
h2U = 0.3
h2Y = 0.2
qX = 0.3
qY = 0.2
a = 0.3
b = 0
M = 5e4
nX = 5e5
nY = 5e5
NormDist = TRUE
df = 5

params = data.frame(pX = pX, pU = pU, pY = pY,
                    h2X = h2X, h2U = h2U, h2Y = h2Y, 
                    qX= qX, qY = qY, a = a, b = b,
                    M = M, nX = nX, nY = nY, 
                    NormDist = NormDist, df = df)
## Save them to file
write.csv(params, "params_settings.csv", row.names=F)
## Source data generation file
source("/data/sgg3/liza/SEM_Sim/Unified_Sim/Scripts/DataSIM.R") # data genetration functions (takes following parameters)
print(paste0("Simulating data ",nD," times..."))

## Set seed for each data generation and save data in folder
for(ID in 1:nD){
  
  set.seed(ID)
  if(ID%%10 == 0) print(paste0("... ", ID))
  
  Data = data_SIM(
    pX = pX, pU = pU, pY = pY,
    h2X = h2X, h2Y = h2Y, h2U = h2U,
    qX = qX, qY = qY, a = a, b = b,
    M = M, nX = nX, nY = nY,
    NormDist = NormDist, df = df)
  
  bX = Data$bX  #Effects of trait X
  bY = Data$bY  #Effects of trait X
  seX = sqrt(1-bX^2)/sqrt(nX)  #Standard error of trait X
  seY = sqrt(1-bY^2)/sqrt(nY)  #Standard error of trait Y
  zX = Data$zX  #Standardised xffects of trait X
  zY = Data$zY  #Standardised xffects of trait X
  D = cbind(bX, bY, seX, seY, zX, zY)
  write.table(D, paste0("Data_", ID, ".csv"), sep=",", quote=F, row.names=F)
}

