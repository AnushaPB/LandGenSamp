#!/bin/bash
# Job name:
#SBATCH --job-name=LGS
#
# Set out/err directories
#SBATCH -o /global/home/users/anushabishop/scratch/stdout/LGS.sh.%J.out
#SBATCH -e /global/home/users/anushabishop/scratch/stderr/LGS.sh.%J.err
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
#SBATCH --mail-type=ALL                                                                             #                                  
#SBATCH --mail-user=anusha.bishop@berkeley.edu  

module load python 
module load gsl 
module load gcc 
module load imagemagick

ipython /global/home/users/anushabishop/GeonomicsAB/LGS/p1_gnxsims/LGS.py > LGS.pyout


