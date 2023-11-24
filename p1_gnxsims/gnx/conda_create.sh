# create conda env for geonomics
conda create --name gnx python=3

# for linux
source activate gnx

# install packages
conda install msprime
pip3 install geopandas rasterio matplotlib scipy bitarray tskit scikit-learn statsmodels psutil nlmpy
pip3 install geonomics==1.3.9

# export env
conda env export > gnx.yml

# to create based on yml
#conda env create -f gnx.yml -n gnx
