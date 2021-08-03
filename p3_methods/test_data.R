source("general_functions.R")

library("here")
library("vcfR")
library("lfmm")

gsd_df <- get_gsd("/Users/Anusha/Documents/GitHub/LandGenSamp/p3_methods/data/mod-K2_phi50_m25_seed1_H50_r60_it--1_t-1000_spp-spp_0.csv")
gen <- get_gen("/Users/Anusha/Documents/GitHub/LandGenSamp/p3_methods/data/mod-K2_phi50_m25_seed1_H50_r60_it--1_t-1000_spp-spp_0.vcf")
loci_df <- read.csv("/Users/Anusha/Documents/GitHub/LandGenSamp/p3_methods/data/nnloci_K2_phi50_m25_seed1_H50_r60.csv")


subIDs <- as.character(samples)
gen <- gen[subIDs,]
gsd_df <- gsd_df[subIDs,]