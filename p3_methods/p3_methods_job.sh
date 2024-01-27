cd p3_methods

# Create folder for outputs
mkdir -p outputs

# Calculate population size and phenotypic mismatch for each simulation
Rscript popsize_count.R
Rscript phenotype_mismatch.R

# Run landscape genomic methods
Rscript run_methods.R

# Compress outputs for git
zip -r outputs.zip outputs
# to decompress: unzip outputs

cd ..