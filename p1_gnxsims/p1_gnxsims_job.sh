# SIMULATED LANDSCAPES -------------------------------------------------------------------
## Render File S2 (environmental correlations)
Rscript -e "rmarkdown::render(here::here('p1_gnxsims', 'MNLM', 'FileS2.Rmd')"

## Create MNLMs
Rscript p1_gnxsims/MNLM/run_MNLM.R
## The generated MNLMs can be visualized with the p1_gnxsims/MNLM/view_MNLM.Rmd notebook

# GEONOMICS SIMULATIONS ------------------------------------------------------------------
cd p1_gnxsims/gnx
## Create genomic architecture for gnx simulations
Rscript create_genomic_architecture.R

## Create and activate conda env
## This conda env was created with the p1_gnxsims/gnx/conda_create.sh script
conda env create -f gnx.yml -n gnx
source activate gnx

## Run geonomics simulation tests (see File S# for more information)
python3 run_gnx_test1.py > run_gnx_test1.pyout
mkdir -p test1
mv GNX_mod-test1* test1

python3 run_gnx_test2.py > run_gnx_test2.pyout
mkdir -p test2
mv GNX_mod-test2* test2

### Render File S# (results of simulation tests)
Rscript -e "rmarkdown::render(here::here('p1_gnxsims', 'gnx', 'FileS#_gnx_test.Rmd')"

## Run full simulations (note: this takes several weeks and is parallelized)
## Contact anusha.bishop@berkeley.edu if you would like the simulation results
python3 run_gnx.py > run_gnx.pyout

# Create folder with simulation results for archive
mkdir -p LGS_simulation_archive
## This only contains the raw csv and vcf mod outputs from geonomics
cp LGS_data/mod* LGS_simulation_archive
## Compress the archive
tar -czvf LGS_simulation_archive.tar.gz LGS_simulation_archive/
## Create directory for dryad data
mkdir -p dryad

cd ..
