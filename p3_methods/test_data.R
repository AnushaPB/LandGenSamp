source("general_functions.R")

library("here")
library("vcfR")
library("lfmm")

gsd_df <- get_gsd("/Users/Anusha/Documents/GitHub/LandGenSamp/p3_methods/data/GNX_mod-K3_phi50_m25_seed1_H5_r30/it--1/spp-spp_0/mod-K3_phi50_m25_seed1_H5_r30_it-1_t-1000_spp-spp_0.csv")
gen <- get_gen("/Users/Anusha/Documents/GitHub/LandGenSamp/p3_methods/data/GNX_mod-K3_phi50_m25_seed1_H5_r30/it--1/spp-spp_0/mod-K3_phi50_m25_seed1_H5_r30_it-1_t-1000_spp-spp_0.vcf")

