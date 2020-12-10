### Script to obtain the SNPs for which we had LD and regression weight calculation from
# the overlap of summary statistics SNPs between UKBiobank and CARDIoGRAMplusC4D
library(data.table)
library(dplyr)
library(DescTools)

## Working directory should be a unique directory to the trait pair to be compared, not the same one as where
# the original data files are held.

fX = fread("~/data/GWAS_CP_all.txt") # Different cohort trait
fY = fread(cmd = " zcat < ~/data/1160.gwas.imputed_v3.both_sexes.tsv.bgz", header=TRUE) # raw file of UKBB trait (from Neale's download)

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
LDrs = rs[selF] ## Selected SNPs above

## Complete columns for fX and fY 
VARinfo = fread(cmd = " zcat < ~/data/variants.tsv.bgz", sep="\t", header=TRUE) ## Variant info file from UKBB GWAS Imputed v3 - File Manifest Release 20180731
fY1 = inner_join(fY, VARinfo[,c(1:6)]) #Joining, by = "variant"
fY1 <- fY1%>% rename(A1 = alt, A2 = ref )

fX1 <- fX %>% rename(rsid = MarkerName, chr = CHR, beta = Beta, se = SE, pval = Pval )
fX1$tstat=fX1$beta/fX1$se
fX1$n_complete_samples = 257828

## Get common SNPs 
commonSNPs = intersect(fX1$rsid,fY1$rsid)
SNPs_to_keep = intersect(commonSNPs, LDrs)
fX2 = fX1[fX1$rsid %in% SNPs_to_keep,]
fY2 = fY1[fY1$rsid %in% SNPs_to_keep,]

dupsX = which(duplicated(fX2$rsid)==TRUE)
dupsY = which(duplicated(fY2$rsid)==TRUE)

if(length(dupsX)!=0){
  dups_index = NA
  ## Remove all SNP duplicates by checking which of each set is closest to the MAF in the LD file
  for(i in dupsX){
    tab=fX2[fX2$rsid==fX2$rsid[i],]
    ld_maf=LDfile[which(LDfile$rs==fX2$rsid[i]),5] #UKBB_AF
    dups_index=c(dups_index,as.numeric(setdiff(rownames(tab),rownames(tab[Closest(tab$minor_AF,as.numeric(ld_maf), which = TRUE),]))))
  }
  ## Filter out unwanted duplicate SNP rows
  fX3=fX2[!rownames(fX2) %in% dups_index, ]
  # Confirm unique SNPs
  length(which(duplicated(fX3$rsid)==TRUE))
}else{fX3 = fX2}

if(length(dupsY)!=0){
  dups_index = NA
  ## Remove all SNP duplicates by checking which of each set is closest to the MAF in the LD file
  for(i in dupsY){
    tab=fY2[fY2$rsid==fY2$rsid[i],]
    ld_maf=LDfile[which(LDfile$rs==fY2$rsid[i]),5] #UKBB_AF
    dups_index=c(dups_index,as.numeric(setdiff(rownames(tab),rownames(tab[Closest(tab$minor_AF,as.numeric(ld_maf), which = TRUE),]))))
  }
  ## Filter out unwanted duplicate SNP rows
  fY3=fY2[!rownames(fY2) %in% dups_index, ]
  # Confirm unique SNPs
  length(which(duplicated(fY3$rsid)==TRUE))
}else{fY3 = fY2}

## Order based on rsid
fX3=fX3[order(fX3$rsid),] 
fY3=fY3[order(fY3$rsid),] 
all(fX3$rsid==fY3$rsid)

#Check alignment between alleles
aligned = which(fY3$A1==fX3$A1 &
                  fY3$A2==fX3$A2)
swapped = which(fY3$A1==fX3$A2 &
                  fY3$A2==fX3$A1)
#Correct the effect of swapped alleles as well as the t-stat for one of the two traits
fX3[swapped,'tstat']=fX3[swapped,'tstat']*-1
fX3[swapped,'beta']=fX3[swapped,'beta']*-1
temp_A1=fX3$A1
fX3[swapped,'A1']=fX3[swapped,'A2']
fX3[swapped,'A2']=temp_A1[swapped]

fX4=fX3[c(aligned,swapped),]
fY4=fY3[c(aligned,swapped),]

## test swapping
length(which(fY4$A1==fX4$A2 & fY4$A2==fX4$A1))

system("mkdir aligned")
#These files will be used as input in LHC Analysis
fwrite(fX4, "./aligned/CP_uniq.tsv", sep="\t")
fwrite(fY4, "./aligned/sleepD_uniq.tsv", sep="\t")

