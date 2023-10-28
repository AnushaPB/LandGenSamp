cd p2_gnxsims

# You must be in p2_sampling to run this file
# Create folder for outputs
mkdir -p outputs

# Create individual based datasets
Rscript envgeo_indsampling.R
Rscript grid_indsampling.R
Rscript random_indsampling.R
Rscript transect_indsampling.R

# Create site based datasets
Rscript envgeo_sitesampling.R
Rscript equi_sitesampling.R
Rscript random_sitesampling.R

cd ..

