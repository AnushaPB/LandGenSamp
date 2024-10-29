Optimizing sampling design for landscape genomics
================
Anusha P. Bishop, Drew E. Terasaki Hart, Ian J. Wang

## 1\. Dependencies

Software versions:

  - R 4.3.0

  - Python 3.9.7

[conda](https://docs.conda.io/en/latest/) was used to manage the python
package dependencies. The conda environment can be recreated using the
p1\_gnxsims/gnx/gnx.yml file (`conda env create -f gnx.yml -n gnx`).

To install all of the necessary R packages run the following:

``` r
source("general_functions.R")
get_packages(install = TRUE)
```

    ##  [1] "here"       "vcfR"       "adegenet"   "stringr"    "dplyr"     
    ##  [6] "tidyr"      "purrr"      "lfmm"       "AssocTests" "gdm"       
    ## [11] "vegan"      "robust"     "qvalue"     "raster"     "hierfstat" 
    ## [16] "tess3r"     "devtools"   "PRROC"      "lme4"       "lmerTest"  
    ## [21] "gt"         "ggplot2"

## 2\. Files

I have tried to create a complete list of files in all of the
directories and what they are for, but shoot me an email if I have
missed anything or you have questions (<anusha.bishop@berkeley.edu>).

### 2.1 [p1\_simulations](https://github.com/AnushaPB/LandGenSamp/tree/main/p1_gnxsims)

Run [geonomics](https://geonomics.readthedocs.io/en/latest/) simulations

    [p1_gnxsims]
    |   p1_gnxsims_job.sh - bash script to create MNLMs and run geonomics simulations
    |
    └───[MNLM]
    │   │   run_MNLM.R - create NLMs for gnx simulations (run by p1_gnxsims_job.sh)
    |   |   MNLM_functions.R - functions for creating and viewing MNLMs
    │   │   view_MNLM.Rmd - view NLMs created by create_MNLM.R
    │   │   FileS2_landscape_parameterization.Rmd - generate environmental correlation distributions, rendered as .html and converted to PDF (rendered by p1_gnxsims_job.sh)
    |   |   FileS2
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
        │   run_gnx.py - run geonomics test simulations (run by p1_gnxsims_job.sh)
        |   FileS4_simulation_equilibrium_tests.Rmd - visualize simulation equilibrium test results (run by p1_gnxsims_job.sh)
        │   data_cleanup.sh - run after parallel_LGS.py to clean up and reorganize the files (run by p1_gnxsims_job.sh)
        │   run_make_dosage.R - run after data_cleanup.sh to create dosage matrices from simulation vcf files (run by p1_gnxsims_job.sh)
        │   
        └───[LGS data] - directory created by data_cleanup.sh and filled with simulation results from run_gnx.py

### 2.2 [p2\_sampling](https://github.com/AnushaPB/LandGenSamp/tree/main/p2_sampling)

Generate datasets for different sampling strategies

    [p2_sampling]
    |   sampling_job.sh - bash script to create sampling datasets
    |   sampling_functions.R - functions used to create different sampling schemes
    |   {sampling scheme}_{ind/sitesampling}.R - all files following this pattern are used to create sampling datasets 
    |   
    └───outputs.zip - compressed folder of p2_sampling outputs

### 2.3 [p3\_methods](https://github.com/AnushaPB/LandGenSamp/tree/main/p3_methods)

Test different landscape genomic analyses

    [p3_methods]
    |   p3_methods_job.sh - bash script to run landscape genomic methods and calculate summary stats for simulation results
    |   general_run_functions.R - functions used for running landscape genomic methods
    |   GEA_functions.R - functions used to run GEA analyses
    |   IBDIBE_functions.R - functions used to run IBD/IBE analyses
    |   run_methods.R - script to run all analyses (run by methods_job.sh)
    |   simulation_stats.R - script to summary statistics for all simulations (run by methods_job.sh)
    |   phenotype_environment.R - script to calculate phenotypic mismatch and phenotype-environment associations (run by methods_job.sh)
    |   genotype_environment.R - script to calculate genotype-environment associations (run by methods_job.sh)
    |  
    └───outputs_*.zip - compressed folder of p3_methods outputs. The outputs are split into multiple zip files due to GitHub maximum file size requirements. To combine them run the following code: cat outputs.zip.part_* > outputs.zip

### 2.4 [p4\_analysis](https://github.com/AnushaPB/LandGenSamp/tree/main/p4_analysis)

Analyze and visualize landscape genomic results

    [p4_analysis]
    |   analysis_functions.R - functions used for running summary analyses
    |   FileS6_GEA.Rmd - notebook for running and visualizing summary analyses for GEA methods (File S4; rendered by p4_analysis.job)
    |   FileS7_IBDIBE.Rmd - notebook for running and visualizing summary analyses for IBD/IBE methods (File S5; rendered by p4_analysis.job)
    |   FileS5_simulation_summary_statistics.Rmd - notebook for visualizing summary statistics for simulations
    |   analysis_megatable.Rmd - notebook for creating giant tables of mixed model results
    |   example_simulations.Rmd - notebook used to visualize example simulations
    |  
    └───[outputs] - directory for storing analysis outputs (empty; filled by running .Rmd files)
    |  
    └───[example_data] - directory for storing example data used by example_simulations.Rmd (empty; example data is moved here by p1_gnxsims_job.sh)

# 3\. Running everything

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
