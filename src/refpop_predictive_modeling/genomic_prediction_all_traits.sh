#!/bin/sh
#=========================================================#
#  script for launching genomic prediction for all traits #
#=========================================================#
sbatch genomic_prediction_traits_color_over.sh
sbatch genomic_prediction_traits_flowering_begin.sh
sbatch genomic_prediction_traits_flowering_end.sh
sbatch genomic_prediction_traits_flowering_full.sh
sbatch genomic_prediction_traits_flowering_intensity.sh
sbatch genomic_prediction_traits_fruit_number.sh
sbatch genomic_prediction_traits_fruit_weight.sh
sbatch genomic_prediction_traits_fruit_weight_single.sh
sbatch genomic_prediction_traits_harvest_date.sh
sbatch genomic_prediction_traits_powdery_mildew.sh
sbatch genomic_prediction_traits_russet_freq_all.sh
sbatch genomic_prediction_traits_sample_size.sh
sbatch genomic_prediction_traits_scab.sh
sbatch genomic_prediction_traits_trunk_diameter.sh
sbatch genomic_prediction_traits_trunk_increment.sh
sbatch genomic_prediction_traits_weight_sample.sh