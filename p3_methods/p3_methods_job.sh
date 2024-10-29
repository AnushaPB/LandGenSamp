cd p3_methods

# Create folder for outputs
mkdir -p outputs

# Calculate genotype-environment/phenotype-environment and other stats for each simulation
Rscript phenotype_environment.R
Rscript genotype_environment.R
Rscript simulation_stats.R

# Run landscape genomic methods
Rscript run_methods.R

# Move outputs from full GDM runs to a seperate folder
mkdir -p outputs/fullGDM
mv outputs/*fullGDM.csv outputs/fullGDM

# Compress outputs for git
zip -r outputs.zip outputs
# split the zip file into 95MB parts for GitHub
split -b 95m outputs.zip outputs.zip.part_
rm outputs.zip
# to decompress: unzip outputs using the following code
# cat outputs.zip.part_* > outputs.zip
# unzip outputs.zip
cd ..