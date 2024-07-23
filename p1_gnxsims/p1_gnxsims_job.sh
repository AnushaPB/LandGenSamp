# SIMULATED LANDSCAPES -------------------------------------------------------------------
## Create MNLMs
Rscript p1_gnxsims/MNLM/run_MNLM.R
## The generated MNLMs can be visualized with the p1_gnxsims/MNLM/view_MNLM.Rmd notebook

## Render File S2 (MNLM parameterization)
Rscript -e "rmarkdown::render(here::here('p1_gnxsims', 'MNLM', 'FileS2.Rmd')"

# GEONOMICS SIMULATIONS ------------------------------------------------------------------
cd p1_gnxsims/gnx
## Create genomic architecture for gnx simulations
Rscript create_genomic_architecture.R

## Create and activate conda env
## This conda env was created with the p1_gnxsims/gnx/conda_create.sh script
conda env create -f gnx.yml -n gnx
source activate gnx

## Run geonomics simulation tests (see File S# for more information)
python3 run_gnx_test.py > run_gnx_test.pyout
mkdir -p test
mv GNX_mod-test* test

### Render File S# (results of simulation tests)
Rscript -e "rmarkdown::render(here::here('p1_gnxsims', 'gnx', 'FileS#_gnx_test.Rmd')"

## Run full simulations (note: this takes several weeks and is parallelized)
python3 run_gnx.py > run_gnx.pyout
mkdir -p LGS_data
mv GNX_mod* LGS_data

# Make everything in the LGS_data folder read-only
sudo chmod -R 744 LGS_data

# Create dosage files
mkdir -p dosage
Rscript run_make_dosage.R

# Create folder with simulation results for archive
mkdir -p LGS_simulation_archive
## This only contains the raw csv and vcf mod outputs  from geonomics from the final time step
find LGS_data -name 'mod*t-6000*' -type f -exec cp --parents \{\} LGS_simulation_archive \;
## Compress the archive
tar -czvf LGS_simulation_archive.tar.gz LGS_simulation_archive/
## Create directory for dryad data
mkdir -p dryad
## Move the archive to the dryad directory
mv LGS_simulation_archive.tar.gz dryad/

cd ..
