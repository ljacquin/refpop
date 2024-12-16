#!/bin/sh
#==================================================#
#  script which executes refpop tasks and analyses #
#==================================================#

# refpop_data_treatment_and_analysis tasks and analyses
cd src/refpop_data_treatment_and_analysis/
R -q --vanilla < refpop_0_raw_phenotype_data_manage_type_correction.R
R -q --vanilla < refpop_1_raw_phenotype_data_outlier_detection_per_env.R
R -q --vanilla < refpop_2_raw_phenotype_data_spat_hetero_correct_and_h2_estim.R
R -q --vanilla < refpop_3_adjusted_blups_lsmeans_phenotypes_computation.R

# refpop_data_structure_analysis tasks and analyses
cd ../refpop_data_structure_analysis/
R -q --vanilla < refpop_genomic_data_structure_analysis.R
R -q --vanilla < refpop_pedigree_and_phenotype_data_structure_analysis.R
R -q --vanilla < refpop_gem_data_structure_analysis.R

# refpop_wiser_genomic_prediction_and_analysis tasks and analyses
cd ../refpop_wiser_genomic_prediction_and_analysis/
sbatch execute_refpop_wiser_genomic_prediction_all_traits.sh