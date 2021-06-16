#!/bin/bash                                                                                                             
# Job name:                                                                                                             
#SBATCH --job-name=smallK_test
#
#Out/err directories
#SBATCH -o /global/home/users/anushabishop/scratch/stdout/smallK_test.sh.%J.out
#SBATCH -e /global/home/users/anushabishop/scratch/stderr/smallK_test.sh.%J.err
#
# Account:                                                                                                              
#SBATCH --account=fc_landgen                                                                                            
#                                                                                                                       
# Partition:                                                                                                            
#SBATCH --partition=savio                                                                                      
#                                                                                                                       
# Wall clock limit:                                                                                                     
#SBATCH --time=3-00:00:00                                                                                               
#                                                                                                                       
#SBATCH --mail-type=ALL                                                                             
#                                                                                                                                         
#SBATCH --mail-user=anusha.bishop@berkeley.edu                                                                                                                                                              

module load python
module load gsl
module load gcc
module load imagemagick

ipython /global/home/users/anushabishop/GeonomicsAB/LandscapeGenomicSampling/p0_tests/smallK_test.py 
