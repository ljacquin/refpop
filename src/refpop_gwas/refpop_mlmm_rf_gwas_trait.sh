#!/bin/sh
#========================================#
#  script for launching gwas for a trait #
#========================================#
### Requirements
#SBATCH --partition=p01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8
#SBATCH --cpus-per-task=8
trait_num_par=$1
Rscript refpop_mlmm_rf_gwas_trait.R $trait_num_par
