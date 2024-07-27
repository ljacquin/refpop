#!/bin/sh
#===============================================================#
#  script for launching wiser genomic prediction for all traits #
#===============================================================#
sbatch refpop_wiser_genomic_prediction_flowering_begin_gaussian_kernel.sh
sbatch refpop_wiser_genomic_prediction_flowering_begin_identity_kernel.sh
sbatch refpop_wiser_genomic_prediction_flowering_begin_linear_kernel.sh
sbatch refpop_wiser_genomic_prediction_harvest_date_gaussian_kernel.sh
sbatch refpop_wiser_genomic_prediction_harvest_date_identity_kernel.sh
sbatch refpop_wiser_genomic_prediction_harvest_date_linear_kernel.sh
sbatch refpop_wiser_genomic_prediction_scab_gaussian_kernel.sh
sbatch refpop_wiser_genomic_prediction_scab_identity_kernel.sh
sbatch refpop_wiser_genomic_prediction_scab_linear_kernel.sh