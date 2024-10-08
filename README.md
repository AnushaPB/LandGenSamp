Optimizing sampling design for landscape genomics
================
Anusha P. Bishop, Drew E. Terasaki Hart, Ian J. Wang

- [1. Dependencies](#1-dependencies)
- [2. Files](#2-files)
  - [2.1 p1_simulations](#21-p1_simulations)
  - [2.2 p2_sampling](#22-p2_sampling)
  - [2.3 p3_methods](#23-p3_methods)
  - [2.4 p4_analysis](#24-p4_analysis)
- [3. Running everything](#3-running-everything)

## 1. Dependencies

Software versions:

- R 4.3.0

- Python 3.9.7

[renv](https://rstudio.github.io/renv/articles/renv.html) was used to
manage R package dependencies. If you clone this repository and open the
project (i.e., the .Rproj file), renv will automatically start up and
ask you if you want to install all the required packages. The packages
can also be installed by running `renv::restore()`.

[conda](https://docs.conda.io/en/latest/) was used to manage the python
package dependencies. The conda environment can be recreated using the
p1_gnxsims/gnx/gnx.yml file (`conda env create -f gnx.yml -n gnx`).

## 2. Files

I have tried to create a complete list of files in all of the
directories and what they are for, but shoot me an email if I have
missed anything or you have questions (<anusha.bishop@berkeley.edu>).

### 2.1 [p1_simulations](https://github.com/AnushaPB/LandGenSamp/tree/main/p1_gnxsims)

Run [geonomics](https://geonomics.readthedocs.io/en/latest/) simulations

    [p1_gnxsims]
    |   p1_gnxsims_job.sh - bash script to create MNLMs and run geonomics simulations
    |
    └───[MNLM]
    │   │   create_MNLM.R - create NLMs for gnx simulations (run by p1_gnxsims_job.sh)
    │   │   view_MNLM.Rmd - view NLMs created by create_MNLM.R
    │   │   File_S2.Rmd - generate environmental correlation distributions, rendered as .html and converted to PDF (rendered by p1_gnxsims_job.sh)
    |   |   
    │   └───[layers] - directory with NLM csvs (used to create landscapes for gnx simulations)
    |   
    └───[gnx]
        │   conda_create.sh - create conda environment for gnx simulations (used to create the gnx.yml file)
        │   create_genomic_architecture.R - create genomic architecture file for gnx simulations (run by p1_gnxsims_job.sh)
        │   FileS3_gnx_parameters.py - supplementary file of gnx parameters
        │   genomic_architecture.csv - genomic architecture file created by create_genomic_architecture.R
        │   gnx.yml - conda yml used to create conda environment (created by conda_create.sh)
        │   run_gnx.py - run geonomics simulations (run by p1_gnxsims_job.sh)
        │   data_cleanup.sh - run after parallel_LGS.py to clean up and reorganize the files (run by p1_gnxsims_job.sh)
        │   run_make_dosage.R - run after data_cleanup.sh to create dosage matrices from simulation vcf files (run by p1_gnxsims_job.sh)
        │   
        └───[LGS data] - directory created by data_cleanup.sh and filled with simulation results from parallel_LGS.py.

### 2.2 [p2_sampling](https://github.com/AnushaPB/LandGenSamp/tree/main/p2_sampling)

Generate datasets for different sampling strategies

    [p2_sampling]
    |   sampling_job.sh - bash script to create sampling datasets
    |   sampling_functions.R - functions used to create different sampling schemes
    |   {sampling scheme}_{ind/sitesampling}.R - all files following this pattern are used to create sampling datasets 
    |   
    └───outputs.zip - compressed folder of p2_sampling outputs

### 2.3 [p3_methods](https://github.com/AnushaPB/LandGenSamp/tree/main/p3_methods)

Test different landscape genomic analyses

    [p3_methods]
    |   methods_job.sh - bash script to run landscape genomic methods and calculate summary stats for simulation results
    |   general_run_functions.R - functions used for running landscape genomic methods
    |   GEA_functions.R - functions used to run GEA analyses
    |   IBDIBE_functions.R - functions used to run IBD/IBE analyses
    |   run_methods.R - script to run all analyses (run by methods_job.sh)
    |   popsize_counts.R - script to calculate population size for all simulations (run by methods_job.sh)
    |   phenotype_environment.R - script to calculate phenotypic mismatch and phenotype-environment associations (run by methods_job.sh)
    |  
    └───outputs.zip - compressed folder of p3_methods outputs

### 2.4 [p4_analysis](https://github.com/AnushaPB/LandGenSamp/tree/main/p4_analysis)

Analyze and visualize landscape genomic results

    [p4_analysis]
    |   analysis_functions.R - functions used for running summary analyses
    |   analysis_GEA.Rmd - notebook for running and visualizing summary analyses for GEA methods (File S4; rendered by p4_analysis.job)
    |   analysis_IBDIBE.Rmd - notebook for running and visualizing summary analyses for IBD/IBE methods (File S5; rendered by p4_analysis.job)
    |   analysis_simulations.Rmd - notebook for visualizing summary statistics for simulations
    |   analysis_megatable.Rmd - notebook for creating giant tables of mixed model results
    |   example_simulations.Rmd - notebook used to visualize example simulations for Figure 1
    |  
    └───[outputs] - directory for storing analysis outputs
    |  
    └───[example_data] - directory for storing example data used by example_simulations.Rmd

# 3. Running everything

**Note: most of the scripts are parallelized and were written to run on
10-25 processors with 126 GB RAM**

``` bash
# Clone repo
git clone https://github.com/AnushaPB/LandGenSamp.git
cd LandGenSamp

# Run gnx simulations
bash p1_gnxsims_job.sh

# Create subsampled datasets
bash p2_sampling_job.sh

# Run landscape genomic methods
bash p3_methods_job.sh

# Run final analyses and generate figures/supplemental files
# Note: see job file for additional information about code notebooks
bash p4_analysis_job.sh
```
