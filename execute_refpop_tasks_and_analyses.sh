#!/bin/sh
#==================================================#
#  script which executes refpop tasks and analyses #
#==================================================#

# refpop_data_treatment_and_analysis tasks and analyses
cd src/refpop_data_treatment_and_analysis/
R -q --vanilla < refpop_0_phenotype_outlier_detection_per_env.R
R -q --vanilla < refpop_1_spat_hetero_correct_per_env_trait_and_h2_estim.R
R -q --vanilla < refpop_2_adjusted_lsmeans_phenotype.R

# refpop_data_structure_analysis tasks and analyses
cd ../refpop_data_structure_analysis/
R -q --vanilla < refpop_genomic_data_structure_analysis.R
R -q --vanilla < refpop_pedigree_and_phenotype_data_structure_analysis.R

# refpop_predictive_modeling_and_analysis tasks and analyses
cd ../refpop_genomic_prediction_and_analysis/
sbatch refpop_genomic_prediction_all_traits.sh