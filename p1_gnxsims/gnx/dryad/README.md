# Data from: Optimizing sampling design for landscape genomics

https://doi.org/10.5061/dryad.63xsj3v8s

This dataset contains a compressed tarball (.tar.gz) with the simulation data used in "Optimizing sampling design for landscape genomics"

## Description of the data and file structure

The tarball must first be unpacked. For example, this can be done using this bash code:
`tar -xzvf LGS_simulation_archive.tar.gz`

The tarball contains 960 pairs of CSV files and Variant Call Format (VCF) files with genomic data for each of the 960 simulations.

Each file is titled as such:
mod-K[1 or 2]_phi[50 or 100]_m[25 or 100]_H[5 or 50]_r[30 or 60]_it-[1-10]_t-6000_spp-spp_0
The values within brackets represent the different low/high parameter levels (e.g., K1 = small population and K2 = large population) or the iteration (1 through 10)

File name abbreviations are:
K = population size
phi = selection strength
m = migration rate
H = spatial autocorrelation
r = environmental correlatio
it = iteration
See the original paper for more information on these parameters.

The CSV files contain geospatial data for the simulated individuals. The columns are:
idx = individual ID (matches with VCF)
z = individual phenotypes ([trait1, trait2])
e = individual environments ([environment 0, environment 1, environment 2])
age = individual age
sex = individual sex
x = individual x coordinate
y = individual y coordinate

The VCF files follow standard VCF formatting and have the same IDs as the CSV files.

## Code/Software

The code and additional files used to generate this data are archived on Zenodo: [DOI HERE]

The most recent version of this code can be found on GitHub: https://github.com/AnushaPB/LandGenSamp