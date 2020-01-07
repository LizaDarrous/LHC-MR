### Script to obtain the SNPs for which we had LD and regression weight calculation from
# the overlap of summary statistics SNPs between UKBiobank and GLGC
library(data.table)
library(dplyr)
library(DescTools)

## Working directory should be a unique directory to the trait pair to be compared, not the same one as where
# the original data files are held.
fileX = "~/data/LDL_uniq.tsv" # GLGC trait
fileY = "~/data/SHeight_uniq.tsv" # UKBB trait

fX = fread(fileX) #fread(cmd = " zcat < ~/data/jointGwasMc_LDL.txt.gz", header=TRUE)
fY = fread(fileY)

#VARinfo = fread(cmd = " zcat < ~/data/variants.tsv.bgz", sep="\t", header=TRUE)
#fX = inner_join(fX, VARinfo[,c(4:6)])  #Joining, by = "rsid" , 4:6 = ref, alt, rsid
#fY = inner_join(fY, VARinfo[,c(4:6)])  #Joining, by = "rsid" , 4:6 = ref, alt, rsid - done in get_tstatXY_UKBB.R

#Get overlapping SNPS between the two traits 
comm=intersect(fY$rsid,fX$rsid)
fX_sub=fX[fX$rsid %in% comm,]
fY_sub=fY[fY$rsid %in% comm,]

fX1=fX[order(fX$rsid),] 
fY1=fY_sub[order(fY_sub$rsid),]

## Add info to GLGC file to match those of UKBB
#Get tsat (standerdised effects) and average sample size to be used as Ny
fX1$tstat=fX1$beta/fX1$se
fX1$n_complete_samples=mean(fX1$N)
fX1$A1=toupper(fX1$A1)
fX1$A2=toupper(fX1$A2)

#Check alignment between alleles
aligned = which(fY1$ref==fX1$A2 &
                  fY1$alt==fX1$A1)
swapped = which(fY1$ref==fX1$A1 &
                  fY1$alt==fX1$A2)
#Correct the effect of swapped alleles
fX1[swapped,'tstat']=fX1[swapped,'tstat']*-1

fX2=fX1[c(aligned,swapped),]
fY2=fY1[c(aligned,swapped),]

#These files will be used as input in LHC Analysis
fwrite(fX2, "LDL_uniq.tsv", sep="\t")
fwrite(fY2, "SHeight_uniq.tsv", sep="\t")

