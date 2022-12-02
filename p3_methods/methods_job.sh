
# Calculate population size and phenotypic mismatch for each simulation
Rscript popsize_count.R
Rscript phenotype_mismatch.R

# Run landscape genomic methods

## LFMM
lfmm_indsampling.R
lfmm_sitesampling.R
lfmm_indsampling_fullK.R
lfmm_sitesampling_fullK.R

## RDA
rda_indsampling.R
rda_sitesampling.R

## MMRR
mmrr_indsampling.R
mmrr_sitesampling.R

## GDM
gdm_indsampling.R
gdm_sitesampling.R
