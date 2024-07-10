#!/bin/sh
#=========================================================#
#  script for launching genomic prediction for all traits #
#=========================================================#
sbatch refpop_genomic_prediction_color_over.sh
sbatch refpop_genomic_prediction_flowering_begin.sh
sbatch refpop_genomic_prediction_flowering_end.sh
sbatch refpop_genomic_prediction_flowering_full.sh
sbatch refpop_genomic_prediction_flowering_intensity.sh
sbatch refpop_genomic_prediction_fruit_number.sh
sbatch refpop_genomic_prediction_fruit_weight.sh
sbatch refpop_genomic_prediction_fruit_weight_single.sh
sbatch refpop_genomic_prediction_harvest_date.sh
sbatch refpop_genomic_prediction_powdery_mildew.sh
sbatch refpop_genomic_prediction_russet_freq_all.sh
sbatch refpop_genomic_prediction_sample_size.sh
sbatch refpop_genomic_prediction_scab.sh
sbatch refpop_genomic_prediction_trunk_diameter.sh
sbatch refpop_genomic_prediction_trunk_increment.sh
sbatch refpop_genomic_prediction_weight_sample.sh