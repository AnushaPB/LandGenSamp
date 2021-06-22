#!/bin/bash
# Job name:
#SBATCH --job-name=LGS_loop
#
# Set out/err directories
#SBATCH -o /global/home/users/anushabishop/scratch/stdout/LGS_loop.sh.%J.out
#SBATCH -e /global/home/users/anushabishop/scratch/stderr/LGS_loop.sh.%J.err
#
# Account:
#SBATCH --account=fc_landgen
#
# Partition:
#SBATCH --partition=savio3_bigmem
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
#SBATCH --mail-type=END,FAIL                                                                         #                                  
#SBATCH --mail-user=anusha.bishop@berkeley.edu  

module load python 
module load gsl 
module load gcc 
module load imagemagick

ipython /global/scratch/anushabishop/LandGenSamp/p1_gnxsims/LGS10k.py > LGS10k.pyout


