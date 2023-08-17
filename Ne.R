####################################################################
############# VALE INSTITUTE OF TECHNOLOGY - ITV ###################
####################################################################

#------------------------------------------------------------------#
#               Effective population size (Ne)                     #
#------------------------------------------------------------------#

#  Script by: Jeronymo Dalapicolla, Jamille Viega & Rodolfo Jaffé

## 1. LOAD THE PACKAGES:
#remotes::install_github("nspope/r2vcftools")
library(r2vcftools)
library(adegenet)
library(vcfR)
library(dartR)
library(poppr)
library(strataG)
library(vcfR)
library(here)

#to install strataG:
library(devtools)
has_devel()
devtools::install_github('ericarcher/strataG', build_vignettes = FALSE, force = T)


## 2. AUXILIARY FUNCTIONS:
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#------------------------------------------------------------------#
#                        Loading Files                             #
#------------------------------------------------------------------#

## 3. Load VCF file
snpsR = read.vcfR("mod-selection_K2_phi100_m100_seed1_H50_r60_it-0_t-100_spp-spp_0.vcf", verbose = T)
snpsR = snpsR[,1:1001]

## 4. Defining populations using the metafile:
#snps =  vcfLink("ipomoea_filtered_ld_hw_neutral.vcf", overwriteID=T)
#VCFsummary(snps) ## 115 individuals and 13167 SNPs
population = rep(as.factor(1), ncol(snpsR@gt) - 1)


#------------------------------------------------------------------#
#                      Converting Files                            #
#------------------------------------------------------------------#

## 5. Converting VCF to Genind
snps_genind = vcfR2genind(snpsR)
class(snps_genind)


## 6. Adding strata (pops) into Genind
snps_genind@pop = population


## 7. Removing Missing Data
#snps_genind = missingno(snps_genind, type = "loci", cutoff = 0)


## 8. Converting Genind to Gtypes
snps_gtypes = genind2gtypes(snps_genind)
class(snps_gtypes)


## 9. Converting Gtypes to GENEPOP to run in NeEstimator:
#genepopWrite(snps_gtypes, "ipomea_genepop.txt")


#------------------------------------------------------------------#
#                      Running Analysis                            #
#------------------------------------------------------------------#

## 10.  Estimating Ne using ldNe
Ne = ldNe(snps_gtypes, maf.threshold = 0, by.strata = TRUE, ci = 0.95, drop.missing = TRUE, num.cores = 4)
Ne

## 11. Save the results:
#write.csv(Ne, "Ne_results_missing.csv")

