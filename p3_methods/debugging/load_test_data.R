#script to just load some test files

library("vcfR")
library("tidyverse")
library("here")

source(here("p3_methods", "general_functions.R"))

loci_df <- read.csv(here("p3_methods", "test_data/old/nnloci_K2_phi50_m100_seed1_H50_r60.csv"))
gen <- get_gen(here("p3_methods", "test_data/old/mod-K2_phi50_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.vcf"))
gsd_df <- get_gsd(here("p3_methods", "test_data/old/mod-K2_phi50_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.csv"))

#sample to make smaller/easier to work with
s <- sample(1:nrow(gen), 2000)
gen <- gen[s,]
gsd_df <- gsd_df[s,]

filepath <- ("test_data/mod-K2_phi50_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.vcf")

ldim <- 100
edge_buffer <- 5
nsite <- 9
inc <- (ldim - 2*edge_buffer)/(sqrt(nsite)-1)
xgrid <- ygrid <- seq(edge_buffer, (ldim - edge_buffer), inc)
cgrid <- expand.grid(xgrid, ygrid)
colnames(cgrid) <- c("x","y")
coordinates(cgrid) <- ~x+y

