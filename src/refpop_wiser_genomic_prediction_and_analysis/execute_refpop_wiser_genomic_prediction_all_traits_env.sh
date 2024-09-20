#!/bin/sh
#===============================================================#
#  script for launching wiser genomic prediction for all traits #
#===============================================================#
n_env=1
n_kernel=1
n_trait=1
for env_num in $( seq 1 1 $n_env )
 do
 for kernel_num in $( seq 1 1 $n_kernel )
  do
  for trait_num in $( seq 1 1 $n_trait )
   do
     sbatch refpop_wiser_genomic_prediction_trait_env.sh $env_num $kernel_num $trait_num
  done
 done
done
