Optimizing sampling design for landscape genomics
================

## 1. Files

I have tried to create a complete list of files in all of the
directories and what they are for, but shoot me an email if I have
missed anything or you have questions (<anusha.bishop@berkeley.edu>).

Software versions:

- Ubuntu 20.04.6

- R 4.3.0

- Python 3.9.7

[renv](https://rstudio.github.io/renv/articles/renv.html) was used to
manage R package dependencies. If you clone this repository and open the
project (i.e., the .Rproj file), renv will automatically start up and
ask you if you want to install all the required packages by running
`renv::restore()`.

*Note: most of the scripts are parallelized and written to run on 10-25
processors with 126 GB RAM*

### 1.1 [p1_simulations](https://github.com/AnushaPB/LandGenSamp/tree/main/p1_gnxsims)

Run [geonomics](https://geonomics.readthedocs.io/en/latest/) simulations

    [p1_gnxsims]
    |
    └───[MNLM]
    │   │   create_MNLM.R - create NLMs for gnx simulations
    │   │   view_MNLM.Rmd - view NLMs created by create_MNLM.R
    |   |   
    │   └───[layers] - directory with NLM csvs (used to create landscapes for gnx simulations)
    |   
    └───[gnx]
        │   create_genomic_architecture.R - create genomic architecture file for gnx simulations
        │   genomic_architecture.csv - genomic architecture file created by create_genomic_architecture.R
        │   run_gnx.py - run geonomics simulations
        │   data_cleanup.sh - run after parallel_LGS.py to clean up and reorganize the files
        │   run_make_dosage.R - run after data_cleanup.sh to create dosage matrices from simulation vcf files
        │   
        └───[LGS data] - directory created by data_cleanup.sh and filled with simulation results from parallel_LGS.py

### 1.2 [p2_sampling](https://github.com/AnushaPB/LandGenSamp/tree/main/p2_sampling)

Generate datasets for different sampling strategies

    [p2_sampling]
    |   sampling_job.sh - bash script to create sampling datasets
    |   sampling_functions.R - functions used to create different sampling schemes
    |   {sampling scheme}_{ind/sitesampling}.R - all files following this pattern are used to create sampling datasets 
    |   
    └───[outputs] - directory for storing sampling outputs

### 1.3 [p3_methods](https://github.com/AnushaPB/LandGenSamp/tree/main/p3_methods)

Test different landscape genomic analyses

    [p3_methods]
    |   methods_job.sh - bash script to run landscape genomic methods and calculate summary stats for simulation results
    |   general_run_functions.R - functions used for running landscape genomic methods
    |   GEA_functions.R - functions used to run GEA analyses
    |   IBDIBE_functions.R - functions used to run IBD/IBE analyses
    |   run_methods.R - script to run all analyses (run by methods_job.sh)
    |   popsize_counts.R - script to calculate population size for all simulations (run by methods_job.sh)
    |   phenotype_mismatch.R - script to calculate phenotypic mismatch (run by methods_job.sh)
    |  
    └───[outputs] - directory for storing methods outputs

### 1.4 [p4_analysis](https://github.com/AnushaPB/LandGenSamp/tree/main/p4_analysis)

Analyze and visualize landscape genomic results

    [p4_analysis]
    |   analysis_functions.R - functions used for running summary analyses
    |   analysis_GEA.Rmd - notebook for running and visualizing summary analyses for GEA methods (File S2)
    |   analysis_IBDIBE.Rmd - notebook for running and visualizing summary analyses for IBD/IBE methods (File S3)
    |   analysis_simulations.Rmd - notebook for visualizing summary statistics for simulations
    |   analysis_megatable.Rmd - notebook for creating giant tables of mixed model results
    |   example_simulations.Rmd - notebook used to visualize example simulations for Figure 1
    |  
    └───[outputs] - directory for storing analysis outputs
    |  
    └───[example_data] - directory for storing example data used by example_simulations.Rmd

# 2. Running everything

``` bash
# clone repo
git clone https://github.com/AnushaPB/LandGenSamp.git
cd LandGenSamp

# Run gnx simulations
bash p1_gnxsims_job.sh

# Create subsampled datasets
bash p2_sampling_job.sh

# Run landscape genomic methods
bash p3_methods_job.sh

# Results can be visualized with the .Rmd notebooks in the p4_analysis folder
```
