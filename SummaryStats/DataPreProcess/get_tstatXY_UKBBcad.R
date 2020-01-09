### Script to obtain the SNPs for which we had LD and regression weight calculation from
# the overlap of summary statistics SNPs between UKBiobank and CARDIoGRAMplusC4D
library(data.table)
library(dplyr)
library(DescTools)

## Working directory should be a unique directory to the trait pair to be compared, not the same one as where
# the original data files are held.
fileX = "~/data/CAD_uniq.tsv" # CARDIo trait
fileY = "~/data/BMI_uniq.tsv" # UKBB trait

fX = fread(fileX) #fread(cmd = " zcat < ~/data/CAD_META.gz.txt.gz", header=TRUE)
fY = fread(fileY)

#VARinfo = fread(cmd = " zcat < ~/data/variants.tsv.bgz", sep="\t", header=TRUE)
#fX = inner_join(fX, VARinfo[,c(4:6)])  #Joining, by = "rsid"/c("rsid"="oldID"), 4:6 = ref, alt, rsid ## CAD has oldID as rsid
#fY = inner_join(fY, VARinfo[,c(4:6)])  #Joining, by = "rsid" , 4:6 = ref, alt, rsid - done in get_tstatXY_UKBB.R

#Get overlapping SNPS between the two traits 
comm=intersect(fY$rsid,fX$oldID)
fX_sub=fX[fX$oldID %in% comm,]
fY_sub=fY[fY$rsid %in% comm,]

fX1=fX_sub[order(fX_sub$oldID),] 
fY1=fY_sub[order(fY_sub$rsid),] 

## Add info to CAD file to match those of UKBB
#Get tsat (standerdised effects) and average sample size to be used as Ny
fX1$tstat=fX1$Effect/fX1$StdErr
fX1$n_complete_samples=380831 #Effective sample size calculated from (4*N_controls*N_cases)/N_control+N_cases
fX1$A1=toupper(fX1$Allele1)
fX1$A2=toupper(fX1$Allele2)
fX1$rsid=fX1$oldID #In analysis file, rsid is needed to match SNPs to LD file


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
fwrite(fX2, "CAD_uniq.tsv", sep="\t")
fwrite(fY2, "BMI_uniq.tsv", sep="\t")

