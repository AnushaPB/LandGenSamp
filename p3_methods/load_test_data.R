#script to just load some test files

library("vcfR")

#Get gen data
get_gen <- function(filepath){
  #read vcf
  vcf <- read.vcfR(filepath)
  #convert to genlight from vcf
  genlight <- vcfR2genlight(vcf) #CHECK THIS
  #convert to matrix
  genmat <- as.matrix(genlight)
  #assign IDs from genlight to matrix rownames
  rownames(genmat) <- genlight@ind.names
  return(genmat)
}

#Get geospatial data
get_gsd <- function(filepath){
  gsd_df <- read.csv(filepath)
  gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
  gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
  rownames(gsd_df) <- gsd_df$idx
  return(gsd_df)
}

loci_df <- read.csv("test_data/nnloci_K4_phi50_m25_seed1_H50_r30.csv")
gen <- get_gen("test_data/mod-K4_phi50_m25_seed1_H50_r30_it--1_t-1000_spp-spp_0.vcf")
gsd_df <- get_gsd("test_data/mod-K4_phi50_m25_seed1_H50_r30_it--1_t-1000_spp-spp_0.csv")

#sample to make smaller/easier to work with
s <- sample(1:nrow(gen), 1000)
gen <- gen[s,]
gsd_df <- gsd_df[s,]
