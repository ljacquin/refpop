#!/bin/sh
#=====================================================#
#  script for launching genomic prediction for traits #
#=====================================================#
### Requirements
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64
#SBATCH --cpus-per-task=20

### Email
### --mail-user=laval.jacquin@inrae.fr
### --mail-type=ALL

R -q --vanilla < refpop_wiser_genomic_prediction_scab_gaussian_kernel.R
