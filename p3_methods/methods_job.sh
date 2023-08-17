# You must be in p3_methods to run this file
# Create folder for outputs
mkdir -p outputs

# Calculate population size and phenotypic mismatch for each simulation
Rscript popsize_count.R
Rscript phenotype_mismatch.R

# Run landscape genomic methods
Rscript run_methods.R