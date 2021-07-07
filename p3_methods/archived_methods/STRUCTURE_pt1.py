import pandas as pd
import ipyrad.analysis as ipa
import toyplot
import matplotlib.pyplot as plt

#following this tutorial: https://ipyrad.readthedocs.io/en/0.9.68/API-analysis/cookbook-structure.html

#Set working directory
import os
os.chdir('/Users/Anusha/Documents/GeonomicsAB/LandscapeGenomicSampling/data')

#read in data (CHANGE ONCE YOU HAVE REAL DATA)
gea_df = pd.read_csv("GNX_mod-islands/it--1/spp-spp_0/mod-islands_it--1_t-2000_spp-spp_0.csv")

#convert format of genomic data
#WRITE CODE TO REMOVE ADAPTIVE LOCI FROM VCF FILE

#init a conversion tool
converter = ipa.vcf_to_hdf5(
    name="modtest",
    data="GNX_mod-islands/it--1/spp-spp_0/mod-islands_it--1_t-2000_spp-spp_0.vcf",
    ld_block_size = 1
)
# run the converter
converter.run()

#STRUCTURE
struct = ipa.structure(
    name="test",
    data="analysis-vcf2hdf5/modtest.snps.hdf5"
)

#set params
struct.mainparams.burnin = 500 #MODIFY
struct.mainparams.numreps = 1000 #MODIFY

#TRUE K (SET BASED ON REAL ANALYSIS)
K = 3

#run stucture
struct.run(nreps=3, kpop=k, auto=True) #MODIFY nreps

#get admixture proportions
table = struct.get_clumpp_table(k)

#export admix proportions
#MODIFY OUTPUT FILE NAME
table.to_csv("outfile.csv", index=None)
