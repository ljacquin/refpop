#!/bin/sh
#==================================================#
#  script which executes refpop tasks and analyses #
#==================================================#

# refpop_data_treatment_and_analysis tasks and analyses
R -q --vanilla < src/refpop_data_treatment_and_analysis/refpop_0_phenotype_outlier_detection_per_env.R
R -q --vanilla < src/refpop_data_treatment_and_analysis/refpop_1_spat_hetero_correct_per_env_trait_and_h2_estim.R
R -q --vanilla < src/refpop_data_treatment_and_analysis/refpop_2_adjusted_lsmeans_phenotype.R

# refpop_data_structure_analysis tasks and analyses
R -q --vanilla < src/refpop_data_structure_analysis/refpop_genomic_data_structure_analysis.R
R -q --vanilla < src/refpop_data_structure_analysis/refpop_pedigree_and_phenotype_data_structure_analysis.R

# refpop_predictive_modeling_and_analysis tasks and analyses
sbatch src/refpop_predictive_modeling_and_analysis/genomic_prediction_all_traits.sh