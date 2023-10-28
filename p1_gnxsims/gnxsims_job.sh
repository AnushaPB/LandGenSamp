# Create MNLMs
Rscript MNLM/create_MNLM.R
# The generated MNLMs can be visualized with the p1_gnxsims/MNLM/view_MNLM.Rmd notebook

# Run gnx simulations
cd gnx
## Create genomic architecture for gnx simulations
Rscript create_genomic_architecture.R

## Create and activate conda env
## This conda env was created with the p1_gnxsims/gnx/conda_create.sh script
conda env create -f gnx.yml
source activate gnx

## Run (note: this takes several weeks and is parallelized)
## Contact anusha.bishop@berkeley.edu if you would like the simulation results
python3 run_gnx.py