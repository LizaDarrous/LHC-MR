### Script to obtain the SNPs for which we had LD and regression weight calculation from
# the overlap of summary statistics SNPs between UKBiobank and CARDIoGRAMplusC4D
library(data.table)
library(dplyr)
library(DescTools)

## Working directory should be a unique directory to the trait pair to be compared, not the same one as where
# the original data files are held.
fX = fread(cmd = " zcat < ~/data/CAD_META.gz.txt.gz", header=TRUE)      # CARDIo trait
fY = fread("~/data/BMI_uniq.tsv") # UKBB trait

#VARinfo = fread(cmd = " zcat < ~/data/variants.tsv.bgz", sep="\t", header=TRUE)
#fY = inner_join(fY, VARinfo[,c(1,4,5)])  #Joining, by = "variant" , 1,4,5,6 = variant, ref, alt, rsid - done in get_tstatXY_UKBB.R, if not then uncomment

## Update column to have consistent comparison
fX$rsid=fX$oldID

#Get overlapping SNPS between the two traits 
comm=intersect(fY$rsid,fX$rsid)
fX_sub=fX[fX$rsid %in% comm,]
fY_sub=fY[fY$rsid %in% comm,]

## Check for duplication and delete if found
length(which(duplicated(fX_sub$rsid)==TRUE)) # 2 dups
length(which(duplicated(fY_sub$rsid)==TRUE))

ALL_dups = which((duplicated(fX_sub$rsid) | duplicated(fX_sub$rsid, fromLast=TRUE))==TRUE) # get ALL indices that are dups
fX_sub2 = fX_sub[-ALL_dups,]

## Second set of overlapping SNPs between the two traits after duplicate deletion
comm2=intersect(fY_sub$rsid,fX_sub2$rsid)
fY_sub2=fY_sub[fY_sub$rsid %in% comm2,]

## Order SNPs
fX1=fX_sub2[order(fX_sub2$rsid),] 
fY1=fY_sub2[order(fY_sub2$rsid),] 
#all(fX1$rsid==fY1$rsid) #check all to be true


## Add info to CAD file to match those of UKBB
#Get tsat (standerdised effects) and average sample size to be used as Ny
fX1$beta=fX1$Effect
fX1$se=fX1$StdErr
fX1$tstat=fX1$beta/fX1$se
fX1$n_complete_samples=380831 #Effective sample size calculated from (4*N_controls*N_cases)/N_control+N_cases
fX1$A1=toupper(fX1$Allele1)
fX1$A2=toupper(fX1$Allele2)

fY1$A1=fY1$alt
fY1$A2=fY1$ref

#Check alignment between alleles
aligned = which(fY1$A1==fX1$A1 &
                  fY1$A2==fX1$A2)
swapped = which(fY1$A1==fX1$A2 &
                  fY1$A2==fX1$A1)
#Correct the effect of swapped alleles for one of the two traits
fX1[swapped,'tstat']=fX1[swapped,'tstat']*-1
fX1[swapped,'beta']=fX1[swapped,'beta']*-1
temp_A1=fX1$A1
fX1[swapped,'A1']=fX1[swapped,'A2']
fX1[swapped,'A2']=temp_A1[swapped]

fX2=fX1[c(aligned,swapped),]
fY2=fY1[c(aligned,swapped),]

## test swapping
length(which(fY2$A1==fX2$A2 & fY2$A2==fX2$A1))

#These files will be used as input in LHC Analysis
fwrite(fX2, "CAD_uniq.tsv", sep="\t")
fwrite(fY2, "BMI_uniq.tsv", sep="\t")

