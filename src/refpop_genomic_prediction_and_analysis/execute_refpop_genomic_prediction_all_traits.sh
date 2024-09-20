#!/bin/sh
#=========================================================#
#  script for launching genomic prediction for all traits #
#=========================================================#
n_trait=15
for trait_num in $( seq 1 1 $n_trait )
 do
  sbatch refpop_genomic_prediction_trait.sh $trait_num
done
