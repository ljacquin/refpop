#!/bin/sh
#===============================================================#
#  script for launching wiser genomic prediction for all traits #
#===============================================================#
n_kernel=1
n_trait=6
for kernel_num in $( seq 1 1 $n_kernel )
 do
 for trait_num in $( seq 1 1 $n_trait )
  do
    sbatch refpop_wiser_genomic_prediction_trait.sh $kernel_num $trait_num
 done
done
