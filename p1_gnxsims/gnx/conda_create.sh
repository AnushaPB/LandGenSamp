# create conda env for geonomics
conda create --name gnx python=3

# for linux
source activate gnx

# install packages
conda install msprime
pip3 install geonomics geopandas rasterio matplotlib scipy bitarray tskit scikit-learn statsmodels psutil nlmpy

# export env
conda env export > gnx.yml

# to recreate conda env based on yml:
#conda env create -f gnx.yml
#source activate gnx

# create conda env for geonomics old version
conda create --name gnx_v1.3.9 python=3

# for linux
source activate gnx

# install packages
conda install msprime
pip3 install geopandas rasterio matplotlib scipy bitarray tskit scikit-learn statsmodels psutil nlmpy
pip3 install geonomics==1.3.9
