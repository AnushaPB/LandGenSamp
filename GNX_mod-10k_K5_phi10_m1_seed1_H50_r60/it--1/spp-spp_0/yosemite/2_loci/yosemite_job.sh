#!/bin/bash
# Job name:
#SBATCH --job-name=yosemite
#
# Set out/err directories
#SBATCH -o /global/home/users/anushabishop/scratch/stdout/yos.sh.%J.out
#SBATCH -e /global/home/users/anushabishop/scratch/stderr/yos.sh.%J.err
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
#SBATCH --mail-type=ALL                                                                             #                                  
#SBATCH --mail-user=anusha.bishop@berkeley.edu  
## Command(s) to run (example):

stdbuf -i0 -o0 -e0 command

module load python 
module load gsl 
module load gcc 
module load imagemagick

ipython /global/home/users/anushabishop/yosemite/run_yosemite_demo.py > yosemite_job.pyout
