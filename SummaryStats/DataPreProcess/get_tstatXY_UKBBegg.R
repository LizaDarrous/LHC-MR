### Script to obtain the SNPs for which we had LD and regression weight calculation from
# the overlap of summary statistics SNPs between UKBiobank and EGG
library(data.table)
library(dplyr)
library(DescTools)

## Working directory should be a unique directory to the trait pair to be compared, not the same one as where
# the original data files are held.
fileX = "~/data/BW-corrected_uniq.tsv"  # EGG trait
fileY = "~/data/MyocInf_uniq.tsv"  # UKBB trait

fX = fread(fileX)
fY = fread(fileY)

#VARinfo = fread(cmd = " zcat < ~/data/variants.tsv.bgz", sep="\t", header=TRUE)
#fY = inner_join(fY, VARinfo[,c(1,4,5)])  #Joining, by = "variant" , 1,4,6 = variant, ref, alt - done in get_tstatXY_UKBB.R, if not then uncomment

#Get overlapping SNPS between the two traits
comm=intersect(fY$rsid,fX$RSID)
fX_sub=fX[fX$RSID %in% comm,]
fY_sub=fY[fY$rsid %in% comm,]

fX1=fX_sub[order(fX_sub$RSID),] 
fY1=fY_sub[order(fY_sub$rsid),] 

## Add info to EGG file to match those of UKBB
#Get tsat (standerdised effects) and average sample size to be used as Ny
fX1$tstat=fX1$beta/fX1$se
fX1$n_complete_samples=mean(fX1$n_ownBW)
fX1$A1=toupper(fX1$ea)
fX1$A2=toupper(fX1$nea)
fX1$rsid=fX1$RSID #In analysis file, rsid is needed to match SNPs to LD file


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
fwrite(fX2, "BW-corrected_uniq.tsv", sep="\t")
fwrite(fY2, "MyocInf_uniq.tsv", sep="\t")

