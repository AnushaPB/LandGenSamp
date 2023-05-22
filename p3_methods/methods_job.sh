
# Calculate population size and phenotypic mismatch for each simulation
Rscript popsize_count.R
Rscript phenotype_mismatch.R

# Run landscape genomic methods

## LFMM
Rscript lfmm_indsampling.R
Rscript lfmm_sitesampling.R
Rscript lfmm_indsampling_fullK.R
Rscript lfmm_sitesampling_fullK.R

## RDA
Rscript rda_indsampling.R
Rscript rda_sitesampling.R

## MMRR
Rscript mmrr_indsampling.R
Rscript mmrr_sitesampling.R

## GDM
Rscript gdm_indsampling.R
Rscript gdm_sitesampling.R
