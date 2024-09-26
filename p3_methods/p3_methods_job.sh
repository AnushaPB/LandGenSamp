cd p3_methods

# Create folder for outputs
mkdir -p outputs

# Calculate genotype-environment/phenotype-environment and other stats for each simulation
Rscript phenotype_environment.R
Rscript genotype_environment.R
Rscript simulation_stats.R

# Run landscape genomic methods
Rscript run_methods.R

# Compress outputs for git
zip -r outputs.zip outputs
# to decompress: unzip outputs

cd ..