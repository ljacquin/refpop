#!/bin/sh
#=================================================================#
#  script for launching genomic prediction for a trait and kernel #
#=================================================================#
### Requirements
#SBATCH --partition=p01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64
#SBATCH --cpus-per-task=12
env_num_par=$1
kernel_num_par=$2
trait_num_par=$3 
Rscript refpop_wiser_genomic_prediction_trait_env.R $env_num_par $kernel_num_par $trait_num_par  
