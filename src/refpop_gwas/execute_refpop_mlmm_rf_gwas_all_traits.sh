#!/bin/sh
#===========================================#
#  script for launching gwas for all traits #
#===========================================#
n_trait=14
for trait_num in $( seq 1 1 $n_trait )
 do
   sbatch refpop_mlmm_rf_gwas_trait.sh $trait_num
done
