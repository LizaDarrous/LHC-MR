### Script to obtain the SNPs for which we had LD and regression weight calculation from UKBiobank summary statistics only!
library(data.table)
library(dplyr)
library(DescTools)

#LD scores and regression weights found in this file were computed based on 3781 individuals from the UK10K study.
#For the LD score calculation, squared correlations were computed between the target SNP and all sequence variants within its 2 Mb neighbourhood 
#with MAF>=0.5 (in the UK10K). For the regression weight calculation, we restricted the neighbouring variants to those part of the 4'650'107 well-imputed SNPs.
LDfile=fread("~/data/LDscores.txt", sep="\t", header=TRUE)

## Extract SNPs high-quality SNPs, defined as being present in both UK10K and UK Biobank, having MAF>1 in both data sets, 
#info>0.99 in the UK Biobank, non-significant (P_{diff}>0.05) allele frequency difference between UK Biobank and UK10K and 
#residing outside the HLA region (chr6:28.5-33.5Mb)
attach(LDfile)
mafT = .005
selF = which(info>.99 & abs(mafUK10K-mafUKBB)<(2*sqrt(1/4e3+1/360e3)) & mafUKBB>mafT & mafUK10K>mafT & 
               !(chr==6 & pos>=28.5e6 & pos<=33.5e6))

Xfile = fread(cmd = " zcat < 21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", header=TRUE) ## Summary statistics file, shoudl be repeated for each trait.
VARinfo = fread(cmd = " zcat < variants.tsv.bgz", sep="\t", header=TRUE) ## Variant info file from UKBB GWAS Imputed v3 - File Manifest Release 20180731

Xfile1 = inner_join(Xfile, VARinfo[,c(1:6)]) #Joining, by = "variant"
LDrs = rs[selF] ## Selected SNPs above

Xfile2=Xfile1[Xfile1$rsid %in% LDrs,]
dups=which(duplicated(Xfile2$rsid)==TRUE) #Some duplication arrises due to joining with VARinfo file that has multiple entery for the same SNP
dups_index = NA

## Remove all SNP duplicates by checking which of each set is closest to the MAF in the LD file
for(i in dups){
  tab=Xfile2[Xfile2$rsid==Xfile2$rsid[i],]
  ld_maf=LDfile[which(LDfile$rs==Xfile2$rsid[i]),4]
  dups_index=c(dups_index,as.numeric(setdiff(rownames(tab),rownames(tab[Closest(tab$minor_AF,as.numeric(ld_maf), which = TRUE),]))))
}

## Filter out unwanted duplicate SNP rows
Xfile3=Xfile2[!rownames(Xfile2) %in% dups_index, ]
# Confirm unique SNPs
length(which(duplicated(Xfile3$rsid)==TRUE))
## Write trait file to be used on LHC-MR analysis
fwrite(Xfile3, "BMI_uniq.tsv", sep="\t")


