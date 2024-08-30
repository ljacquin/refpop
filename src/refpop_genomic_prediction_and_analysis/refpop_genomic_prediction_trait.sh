#!/bin/sh
#======================================================#
#  script for launching genomic prediction for a trait #
#======================================================#
### Requirements
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64
#SBATCH --cpus-per-task=12
trait_num_par=$1
Rscript refpop_genomic_prediction_trait.R $trait_num_par
