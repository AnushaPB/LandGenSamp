Landscape genomic sampling project
================

## 1. Main Pipeline

### 1.1 [p1_simulations](https://github.com/AnushaPB/LandGenSamp/tree/main/p1_gnxsims)

Run [geonomics](https://geonomics.readthedocs.io/en/latest/) simulations

### 1.2 [p2_sampling](https://github.com/AnushaPB/LandGenSamp/tree/main/p2_sampling)

Generate datasets for different sampling strategies

### 1.3 [p3_methods](https://github.com/AnushaPB/LandGenSamp/tree/main/p3_methods)

Test different landscape genomic analyses

### 1.4 [p4_analysis](https://github.com/AnushaPB/LandGenSamp/tree/main/p4_analysis)

Analyze and visualize landscape genomic results

## To run everything:

``` bash
# Create MNLMs
Rscript p1_gnxims/MNLM/create_MNLM.R
# The generated MNLMs can be visualized with the p1_gnxsims/MNLM/view_MNLM.Rmd notebook

# Run gnx simulations
cd p1_gnxsims/gnx
## Create genomic architecture for gnx simulations
Rscript create_genomic_architecture.R
## Run simulations (note: this takes several weeks and is parallelized)
## contact anusha.bishop@berkeley.edu if you would like the simulation results
bash parallel_LGS_job.sh

# Create subsampled datasets
cd ../p2_gnxsims
bash sampling_job.sh

# Run landscape genomic methods
cd ../p3_methods
bash methods_job.sh

# Results can be visualized with the Rmd notebooks in the p4_analysis folder
```
