# SIMULATED LANDSCAPES -------------------------------------------------------------------
## Create MNLMs
Rscript p1_gnxsims/MNLM/run_MNLM.R
## The generated MNLMs can be visualized with the p1_gnxsims/MNLM/view_MNLM.Rmd notebook

## Render File S2 (MNLM parameterization)
Rscript -e "rmarkdown::render(here::here('p1_gnxsims', 'MNLM', 'FileS2_landscape_parameterization.Rmd')"

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
Rscript -e "rmarkdown::render(here::here('p1_gnxsims', 'gnx', 'FileS4_simulation_equilibrium_tests.Rmd')"

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
# Copy the files directly into LGS_simulation_archive without creating LGS_simulation_archive/LGS_data
find LGS_data -name 'mod*t-6000*' -type f -exec cp --parents \{\} LGS_simulation_archive/ \;
# Move the files from LGS_simulation_archive/LGS_data to LGS_simulation_archive
find LGS_simulation_archive/LGS_data -type f -exec mv \{\} LGS_simulation_archive/ \;
# Remove the now-empty LGS_simulation_archive/LGS_data directory
rm -r LGS_simulation_archive/LGS_data
# Compress the archive
tar -czvf LGS_simulation_archive.tar.gz LGS_simulation_archive/
# Create directory for dryad data
mkdir -p dryad
# Move the archive to the dryad directory
mv LGS_simulation_archive.tar.gz dryad/

# Copy example data into p4_analysis folder
mkdir -p ../../p4_analysis/example_data
## List of files to copy with their respective folder structures
files=(
  "GNX_mod-K2_phi100_m100_seed1_H50_r60/it-0/spp-spp_0/mod-K2_phi100_m100_seed1_H50_r60_it-0_t-6000_spp-spp_0.csv"
  "GNX_mod-K1_phi50_m100_seed1_H50_r60/it-0/spp-spp_0/mod-K1_phi50_m100_seed1_H50_r60_it-0_t-6000_spp-spp_0.csv"
  "GNX_mod-K1_phi100_m25_seed1_H50_r60/it-0/spp-spp_0/mod-K1_phi100_m25_seed1_H50_r60_it-0_t-6000_spp-spp_0.csv"
  "GNX_mod-K1_phi100_m100_seed1_H5_r60/it-0/spp-spp_0/mod-K1_phi100_m100_seed1_H5_r60_it-0_t-6000_spp-spp_0.csv"
  "GNX_mod-K1_phi100_m100_seed1_H50_r30/it-0/spp-spp_0/mod-K1_phi100_m100_seed1_H50_r30_it-0_t-6000_spp-spp_0.csv"
  "GNX_mod-K1_phi100_m100_seed1_H50_r60/it-0/spp-spp_0/mod-K1_phi100_m100_seed1_H50_r60_it-0_t-6000_spp-spp_0.csv"
)

# Loop through each file and copy it
for file in "${files[@]}"; do
  # Copy the file to the analysis directory
  cp LGS_data/"$file" "../../p4_analysis/example_data"
done

cd ..
