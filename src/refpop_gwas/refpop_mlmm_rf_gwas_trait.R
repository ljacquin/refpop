# script meant to perform wiser and adjusted ls-means mlmm and random forest gwas
# for the refpop
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
if ("refpop_env" %in% conda_list()$name) {
  print("using refpop_env")
  use_condaenv("refpop_env")
}
# install other requirements from github if necessary
install_other_requirements <- F
if (install_other_requirements) {
  # reticulate::install_miniconda()
  conda_create("refpop_env")
  use_condaenv("refpop_env")
  library(devtools)
  devtools::install_github("ljacquin/KRMM")
  devtools::install_github("rstudio/tensorflow")
  library(tensorflow)
  install_tensorflow(envname = "refpop_env")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
  install.packages("umap")
}
use_tensorflow_ <- F
if (use_tensorflow_) {
  # leave tensorflow and keras for later use
  library(tensorflow)
  library(keras3)
  tensorflow::tf$random$set_seed(0)
  py_module_available("keras") # must return TRUE
  py_module_available("tensorflow") # must return TRUE
  py_discover_config("keras") # more info on the python env, tf and keras
}
library(umap)
library(MASS)
library(data.table)
library(stringr)
library(lme4)
library(tidyr)
library(FactoMineR)
library(doParallel)
library(doRNG)
library(robustbase)
library(foreach)
library(parallel)
library(missForest)
library(Matrix)
library(matrixcalc)
library(rgl)
library(Rfast)
library(cvTools)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(KRMM)
library(kernlab)
library(whitening)
library(glmnet)
library(ranger)
library(mixOmics)
library(future)
library(future.apply)
library(mlmm.gwas)
library(fastqq)
library(evir)

# set a seed for reproducibility
set.seed(123)

# get all mlmm.gwas from namespace and make them globally visible
list_mlmm_gwas_fcts <- ls(getNamespace("mlmm.gwas"), all.names = TRUE)[-1]
for (f_ in list_mlmm_gwas_fcts) {
  assign(f_, get(f_, envir = getNamespace("mlmm.gwas")),
    envir = .GlobalEnv
  )
}

# define computation mode, i.e. local or cluster
computation_mode <- "cluster"

# if comutations are local in rstudio, detect and set script path
# automatically using rstudioapi
if (identical(computation_mode, "local")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}

# source functions
source("../functions.R")

# set options
options(future.globals.maxSize = 60 * 1024^3)
options(expressions = 5e5)
options(warn = -1)

# define mlmm maximum number of steps for forward stepwise regression
max_step_mlmm_fwd_reg_ <- 10

# define alpha for corrected Bonferroni threshold
alpha_bonf_ <- 0.05

# define probabilities for peak over threshold (POT) approach using generalized
# pareto distribution (GPD)
prob_basis_tresh_ <- 0.99
prob_qgpd_tresh_ <- 0.999

# should mtry be optimized for ranger ?
tune_mtry_ranger_ <- T

# define number of recomputations for normalized variable importance
n_recompute_nvi <- 100

# set input paths
geno_dir_path <- "../../data/genotype_data/"
pheno_dir_path <- "../../data/phenotype_data/"

# set path for wiser phenotypes estimated using whitening
wiser_pheno_dir_path <- "../../data/phenotype_data/wiser_phenotype_estimates/"

# output result path for genotype graphics
output_gwas_results_path <- "../../results/gwas/"
output_gwas_graphics_path <- "../../results/graphics/gwas_graphics/"

# define list of all traits
list_all_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew"
)

# should wiser and genomic prediction be performed for a subset of traits ?
# this can be useful for managing computation resources, and most importantly
# for analyzing specific traits without modifying the original list of all traits
use_subset_list <- F

# define traits_ for wiser
if (use_subset_list) {
  traits_ <- c(
    "Color_over", "Powdery_mildew"
  )
} else {
  traits_ <- list_all_traits_
}

# get trait arguments
args <- commandArgs(trailingOnly = T)
trait_num <- as.integer(args[1])

# define trait_
trait_ <- traits_[trait_num]
print(paste0("trait: ", trait_))

# create result repository for trait_ if it does not exist
if (!dir.exists(paste0(output_gwas_graphics_path, trait_))) {
  dir.create(paste0(output_gwas_graphics_path, trait_))
  dir.create(paste0(output_gwas_graphics_path, trait_, "/wiser"))
  dir.create(paste0(output_gwas_graphics_path, trait_, "/ls_means"))
}

# get wiser and adjusted ls-means phenotypes, and genotype data (genomic data)

# get trait wiser phenotypes
wiser_trait_ <- readRDS(paste0(
  wiser_pheno_dir_path,
  paste0("wiser_obj_linear_kernel_", trait_)
))$wiser_phenotypes
colnames(wiser_trait_)[2] <- trait_

# get trait ls-means phenotypes
ls_means_trait_ <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "adjusted_ls_mean_phenotypes.csv"
)))[, c("Genotype", trait_)]

# get marker data
marker_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_data.csv"
)))

# get map data
map_df <- as.data.frame(readRDS(paste0(
  geno_dir_path,
  "phased_data/chr_map.RDS"
)))

# remove markers with maf inferior to 5% (default value) to prevent loss of
# power and false associations for under-represented alleles
marker_df <- filter_maf(marker_df)

# merge ls_means_trait_ and marker_df for integrity of analyses
merged_df <- merge(wiser_trait_, marker_df, by = "Genotype")
wiser_trait_ <- merged_df[, c("Genotype", trait_)]

# get common genotype marker data, rename rows and remove genotype and trait_
# columns
geno_names <- merged_df$Genotype
marker_df <- merged_df[, -match(c("Genotype", trait_), colnames(merged_df))]
rownames(marker_df) <- geno_names

# merge ls_means_trait_ and wiser_trait_ for integrity of analyses
merged_df <- merge(ls_means_trait_, wiser_trait_, by = "Genotype")
ls_means_trait_ <- merged_df[, c("Genotype", paste0(trait_, ".x"))]
colnames(ls_means_trait_)[2] <- trait_

# free some space
rm(merged_df)

# get common markers between map and marker data, and order them
marker_df_snp_names_ <- data.frame(
  "index" = 1:ncol(marker_df),
  "snp.name" = colnames(marker_df)
)
common_snp_names_ <- merge(
  map_df, marker_df_snp_names_,
  by = "snp.name"
)
common_snp_names_ <- common_snp_names_[order(
  common_snp_names_$chromosome,
  common_snp_names_$position
), ]

# reorder the markers (if necessary), and verify equality for integrity
marker_df <- marker_df[, common_snp_names_$snp.name]
map_df <- map_df[match(common_snp_names_$snp.name, map_df$snp.name), ]
map_df <- map_df[, c("snp.name", "chromosome", "position")]
identical(map_df$snp.name, colnames(marker_df))

# calculate corrected Bonferroni threshold
threshold_ <- alpha_bonf_ / ncol(marker_df)

# compute gram matrix
marker_matrix <- scale(apply(marker_df, 2, as.numeric), scale = F)
k_mat <- marker_matrix %*% t(marker_matrix)
rownames(marker_matrix) <- rownames(marker_df)
colnames(marker_matrix) <- colnames(marker_df)
colnames(k_mat) <- rownames(k_mat) <- rownames(marker_df)

# -- reformat response for wiser
wiser_y <- wiser_trait_[, trait_]
names(wiser_y) <- wiser_trait_$Genotype

# ---- mlmm gwas for wiser

# fit mlmm forward stepwise regression algorithm for wiser phenotypes
mlmm_obj_wiser_ <- mlmm_allmodels_(
  wiser_y, list(marker_df), list(k_mat),
  maxsteps = max_step_mlmm_fwd_reg_
)

# get mlmm algorithm fitted components (i.e. p-values and degrees of freedom)
mlmm_p_val_wiser_ <- mlmm_obj_wiser_$p_val
mlmm_deg_freedom_wiser_ <- mlmm_obj_wiser_$df_list

# get step associated to "optimal" model based on ebic minimization
res_eBIC <- as.data.frame(eBIC_allmodels(
  wiser_y,
  frommlmm_toebic(list(marker_matrix), mlmm_p_val_wiser_),
  list(k_mat), ncol(marker_matrix)
))
optimal_model_associated_step_ <- which.min(
  res_eBIC[, str_detect(colnames(res_eBIC), "eBIC")][-1]
)

# make a manhattan plot
manhattan_plot_done_ <- mlmm_manhattan_plot_(
  mlmm_p_val_wiser_,
  map = map_df,
  steps = optimal_model_associated_step_,
  output_path_ = paste0(output_gwas_graphics_path, trait_, "/wiser/"),
  file_name_ = paste0(trait_, "_wiser_mlmm_manhattan_plot.png"),
  threshold_ = -log10(threshold_),
  main_ =
    paste0(
      "MLMM GWAS of ", trait_,
      " WISER phenotypes using a corrected Bonferroni threshold for 0.05"
    )
)
manhattan_plot_done_

# get p-values associated with the step of the identified "optimal" model (i.e.,
# the model minimizing the extended-BIC). Note that these p-values do not correspond
# to the "optimal" model itself but to other models of the same form under H1 that
# were not significant at the step of the identified "optimal" model.
# The significant p-values specifically associated with the "optimal" model under
# H1 are those corresponding to the SNPs included in this model.
p_val_step_opt_model_ <- mlmm_p_val_wiser_[[
  optimal_model_associated_step_ + 1
]]
head(sort(p_val_step_opt_model_, decreasing = F))
if (sum(str_detect(names(p_val_step_opt_model_), "selec_")) > 0) {
  names(p_val_step_opt_model_) <- str_replace_all(
    names(p_val_step_opt_model_),
    "selec_",
    replacement = ""
  )
}

# get degrees of freedom of "optimal" model (i.e. minimizing extended-BIC). In
# other words, get df1 = q = p1 (nb. param. under H1) - p0 (nb. param. under H0)
# which equals 1 in case of the nullity test of q = 1 variable (SNP) under H0,
# and df2 = n - p1 which decreases by 1 at each step (i.e. p1 increases by 1 at
# each step).
deg_freed_step_opt_model_ <- sort(unique(mlmm_deg_freedom_wiser_[[
  optimal_model_associated_step_ + 1
]]))

# make mlmm gwas qqplot for trait WISER estimated phenotypes
qq_plot_done_ <- qqplot_f_distrib_mlmm_gwas_(
  p_val_obs_ = p_val_step_opt_model_,
  df1 = as.numeric(deg_freed_step_opt_model_[1]),
  df2 = as.numeric(deg_freed_step_opt_model_[2]),
  seed_ = 42,
  main_ = paste0(
    "MLMM GWAS Q-Q plot for ",
    trait_, " WISER estimated phenotypes"
  ),
  output_path_ = paste0(output_gwas_graphics_path, trait_, "/wiser/"),
  file_name_ = paste0(trait_, "_wiser_mlmm_qq_plot.png")
)
qq_plot_done_

# get significant snps according to defined corrected Bonferroni threshold
mlmm_signif_snps_wiser_ <- p_val_step_opt_model_[
  which(p_val_step_opt_model_ < threshold_)
]

if (length(mlmm_signif_snps_wiser_) > 0) {
  # reformat significant snp(s) data
  mlmm_signif_snps_wiser_ <- data.frame(
    "snp_name" = names(mlmm_signif_snps_wiser_),
    "p_value" = as.numeric(mlmm_signif_snps_wiser_),
    "chromosome_number" = map_df$chromosome[
      match(names(mlmm_signif_snps_wiser_), map_df$snp.name)
    ],
    "position_in_bp" = map_df$position[
      match(names(mlmm_signif_snps_wiser_), map_df$snp.name)
    ]
  )

  # write significant snp(s) infos
  fwrite(mlmm_signif_snps_wiser_, file = paste0(
    paste0(output_gwas_graphics_path, trait_, "/wiser/"),
    trait_, "_wiser_mlmm_significant_marker_info.csv"
  ))

  # get the snps with highest significancy by chromosome
  mlmm_highest_signif_snps_wiser_ <- mlmm_signif_snps_wiser_ %>%
    group_by(chromosome_number) %>%
    filter(p_value == min(p_value)) %>%
    ungroup()

  # make boxplots for trait phenotypes according to SNP allele combinations (
  # i.e. homozygote and heterozygote ) of genotypes for SNP markers with highest
  # score (i.e. p-value or normalized variable importance)
  for (snp_ in mlmm_highest_signif_snps_wiser_$snp_name) {
    # get marker score for boxplot
    score_ <- signif(
      as.numeric(
        mlmm_highest_signif_snps_wiser_$p_value[
          match(snp_, mlmm_highest_signif_snps_wiser_$snp_name)
        ]
      ), 3
    )
    df_ <- data.frame(
      "wiser_y" = wiser_y, "marker" = marker_matrix[, snp_]
    )
    r2_ <- summary(
      lm(as.formula("wiser_y ~ 1 + marker"), data = df_)
    )$r.squared
    r2_ <- signif(as.numeric(r2_), 3)
    geno_box_plot_ <- genotypes_boxplot_(
      X = marker_df,
      Y = wiser_y,
      ylab_ = paste0(trait_, " WISER phenotypes"),
      xlab_ = paste0(
        "Genotype alleles combination on chromosome ",
        map_df$chromosome[match(snp_, map_df$snp.name)]
      ),
      markers = snp_,
      r2_ = r2_,
      score_ = score_,
      score_text_ = "p-value",
      output_path_ = paste0(output_gwas_graphics_path, trait_, "/wiser/"),
      file_name_ = paste0(
        trait_, "_wiser_mlmm_highest_signif_marker_", snp_, "_chromo_",
        map_df$chromosome[match(snp_, map_df$snp.name)], ".png"
      )
    )
  }
}

# ---- rf gwas for wiser

# tune mtry before fitting ranger
if (tune_mtry_ranger_) {
  rf_obj_opt_mtry_ <- tune_mtry_ranger_rf(
    X = marker_matrix, Y = wiser_y,
    mtry_grid_ = floor(ncol(marker_matrix) / c(10, 5, 3)),
    num_trees_ = 500,
    pkgs_to_export_ = "ranger"
  )
  opt_mtry <- rf_obj_opt_mtry_$opt_mtry
} else {
  # default to classical value used in rf implementations (which works pretty
  # well in many cases)
  opt_mtry <- floor(ncol(marker_matrix) / 3)
}

# fit rf using ranger
rf_obj_wiser_ <- ranger(
  y = wiser_y,
  x = marker_matrix,
  mtry = opt_mtry,
  num.trees = 500,
  importance = "permutation"
)

# get normalized variable importance (nvi) for markers
var_imp_wiser_ <- importance(
  rf_obj_wiser_
)

# normalize importance values using min-max normalization
var_imp_wiser_ <- min_max_normalization(var_imp_wiser_)

# apply a peak over threshold (POT) approach, using the generalized pareto
# distribution (GPD), to model the tail of the normalized variable importance
# distribution which has extreme values
quant_basis_tresh_ <- quantile(var_imp_wiser_, prob_basis_tresh_)
gpd_fit <- gpd(var_imp_wiser_, threshold = quant_basis_tresh_)
xi_ <- gpd_fit$par.ests["xi"] # shape parameter
beta_ <- gpd_fit$par.ests["beta"] # scale parameter
quant_qgpd_tresh_ <- as.numeric(
  qgpd(prob_qgpd_tresh_, xi = xi_, beta = beta_)
)

if (length(var_imp_wiser_) > 0) {
  # get the subset of variables (i.e. markers) detected as important
  sub_var_imp_wiser_ <- var_imp_wiser_[
    var_imp_wiser_ > quant_basis_tresh_
  ]

  if (length(sub_var_imp_wiser_) > 0) {
    sub_mat_ <- marker_matrix[, names(sub_var_imp_wiser_)]

    # recompute normalized variable importance several times for subset of variables
    # detected as important
    list_sub_var_imp_wiser <- vector("list", n_recompute_nvi)
    for (recomp_idx_ in 1:n_recompute_nvi) {
      list_sub_var_imp_wiser[[recomp_idx_]] <- min_max_normalization(
        importance(
          ranger(
            y = wiser_y,
            x = sub_mat_,
            mtry = floor(ncol(sub_mat_) / 3),
            num.trees = 500,
            importance = "permutation"
          )
        )
      )
    }
    mean_sub_var_imp_wiser <- Reduce(`+`, list_sub_var_imp_wiser) /
      length(list_sub_var_imp_wiser)

    # reassign mean nvi into original vector of nvi for variables detected as
    # important
    var_imp_wiser_[
      match(names(mean_sub_var_imp_wiser), names(var_imp_wiser_))
    ] <- mean_sub_var_imp_wiser
  }
}

# create rf result data frame specific to manhattan ggplot
rf_df_results_ <- data.frame(
  "snp" = map_df$snp.name,
  "chr" = map_df$chromosome,
  "bp" = map_df$position,
  "score" = var_imp_wiser_
)

# make a manhattan ggplot
manhattan_ggplot_plotly_done_ <- manhattan_ggplot_plotly_(
  df_ = rf_df_results_,
  trait_ = trait_,
  qgpd_threshold = quant_qgpd_tresh_,
  x_lab_ = "SNP index",
  y_lab_ = "Normalized variable importance (NVI)",
  title_ = paste0(
    "Random forest GWAS of ", trait_,
    " WISER phenotypes using a 0.001 threshold for a generalized Pareto distribution (GPD)"
  ),
  title_size_ = 16,
  axis_title_size_ = 14,
  output_path_ = paste0(output_gwas_graphics_path, trait_, "/wiser/"),
  file_name_ = paste0(trait_, "_wiser_rf_manhattan_plot.png"),
  make_ggplotly_ = F
)
manhattan_ggplot_plotly_done_

# get significant snps according to generalized pareto distribution's
# 0.1% threshold
rf_signif_snps_wiser_ <- var_imp_wiser_[
  var_imp_wiser_ > quant_qgpd_tresh_
]

if (length(rf_signif_snps_wiser_) > 0) {
  # if any significant snps, convert their infos into a data frame containing
  # snp name, score and chromosome number
  rf_signif_snps_wiser_df_ <- data.frame(
    "snp_name" = names(rf_signif_snps_wiser_),
    "normalized_variable_importance" = rf_signif_snps_wiser_,
    "chromosome_number" = map_df$chromosome[
      match(names(rf_signif_snps_wiser_), map_df$snp.name)
    ],
    "position_in_bp" = map_df$position[
      match(names(rf_signif_snps_wiser_), map_df$snp.name)
    ]
  )

  # get and write significant markers infos
  fwrite(rf_signif_snps_wiser_df_, file = paste0(
    paste0(output_gwas_graphics_path, trait_, "/wiser/"),
    trait_, "_wiser_rf_significant_marker_info.csv"
  ))

  # get the snps with highest importance by chromosome
  rf_highest_signif_snps_wiser_df_ <- rf_signif_snps_wiser_df_ %>%
    group_by(chromosome_number) %>%
    filter(normalized_variable_importance ==
      max(normalized_variable_importance)) %>%
    ungroup()

  # make boxplots for trait phenotypes according to SNP allele combinations (
  # i.e. homozygote and heterozygote ) of genotypes for SNP markers with highest
  # score (i.e. p-value or normalized variable importance)
  for (snp_ in rf_highest_signif_snps_wiser_df_$snp_name) {
    # get marker score for boxplot
    score_ <- signif(
      as.numeric(
        rf_highest_signif_snps_wiser_df_$normalized_variable_importance[
          match(snp_, rf_highest_signif_snps_wiser_df_$snp_name)
        ]
      ), 3
    )
    df_ <- data.frame(
      "wiser_y" = wiser_y, "marker" = marker_matrix[, snp_]
    )
    r2_ <- ranger(
      as.formula("wiser_y ~ 1 + marker"),
      data = df_
    )$r.squared
    r2_ <- signif(as.numeric(r2_), 3)
    geno_box_plot_ <- genotypes_boxplot_(
      X = marker_df,
      Y = wiser_y,
      ylab_ = paste0(trait_, " WISER phenotypes"),
      xlab_ = paste0(
        "Genotype alleles combination on chromosome ",
        map_df$chromosome[match(snp_, map_df$snp.name)]
      ),
      markers = snp_,
      r2_ = r2_,
      score_ = score_,
      score_text_ = "NVI",
      output_path_ = paste0(output_gwas_graphics_path, trait_, "/wiser/"),
      file_name_ = paste0(
        trait_, "_wiser_rf_highest_signif_marker_", snp_, "_chromo_",
        map_df$chromosome[match(snp_, map_df$snp.name)], ".png"
      )
    )
  }
}

# -- reformat response for ls-means
ls_means_y <- ls_means_trait_[, trait_]
names(ls_means_y) <- ls_means_trait_$Genotype

# ---- mlmm gwas for ls-means

# fit mlmm forward stepwise regression algorithm for ls-means phenotypes
mlmm_obj_ls_means_ <- mlmm_allmodels_(
  ls_means_y, list(marker_df), list(k_mat),
  maxsteps = max_step_mlmm_fwd_reg_
)

# get mlmm algorithm fitted components (i.e. p-values and degrees of freedom)
mlmm_p_val_ls_means_ <- mlmm_obj_ls_means_$p_val
mlmm_deg_freedom_ls_means_ <- mlmm_obj_ls_means_$df_list

# get step associated to "optimal" model based on ebic minimization
res_eBIC <- as.data.frame(eBIC_allmodels(
  ls_means_y,
  frommlmm_toebic(list(marker_matrix), mlmm_p_val_ls_means_),
  list(k_mat), ncol(marker_matrix)
))
optimal_model_associated_step_ <- which.min(
  res_eBIC[, str_detect(colnames(res_eBIC), "eBIC")][-1]
)

# make a manhattan plot
manhattan_plot_done_ <- mlmm_manhattan_plot_(
  mlmm_p_val_ls_means_,
  map = map_df,
  steps = optimal_model_associated_step_,
  output_path_ = paste0(output_gwas_graphics_path, trait_, "/ls_means/"),
  file_name_ = paste0(trait_, "_ls_means_mlmm_manhattan_plot.png"),
  threshold_ = -log10(threshold_),
  main_ =
    paste0(
      "MLMM GWAS of ", trait_,
      " LS-means phenotypes using a corrected Bonferroni threshold for 0.05"
    )
)
manhattan_plot_done_

# get p-values associated with the step of the identified "optimal" model (i.e.,
# the model minimizing the extended-BIC). Note that these p-values do not correspond
# to the "optimal" model itself but to other models of the same form under H1 that
# were not significant at the step of the identified "optimal" model.
# The significant p-values specifically associated with the "optimal" model under
# H1 are those corresponding to the SNPs included in this model.
p_val_step_opt_model_ <- mlmm_p_val_ls_means_[[
  optimal_model_associated_step_ + 1
]]
head(sort(p_val_step_opt_model_, decreasing = F))

# if any, reformat snp names which do not respect the nomenclature to the
# expected standard
if (sum(str_detect(names(p_val_step_opt_model_), "selec_")) > 0) {
  names(p_val_step_opt_model_) <- str_replace_all(
    names(p_val_step_opt_model_),
    "selec_",
    replacement = ""
  )
}

# get degrees of freedom of "optimal" model (i.e. minimizing extended-BIC). In
# other words, get df1 = q = p1 (nb. param. under H1) - p0 (nb. param. under H0)
# which equals 1 in case of the nullity test of q = 1 variable (SNP) under H0,
# and df2 = n - p1 which decreases by 1 at each step (i.e. p1 increases by 1 at
# each step).
deg_freed_step_opt_model_ <- sort(unique(mlmm_deg_freedom_ls_means_[[
  optimal_model_associated_step_ + 1
]]))

# make mlmm gwas qqplot for trait LS-means estimated phenotypes
qq_plot_done_ <- qqplot_f_distrib_mlmm_gwas_(
  p_val_obs_ = p_val_step_opt_model_,
  df1 = as.numeric(deg_freed_step_opt_model_[1]),
  df2 = as.numeric(deg_freed_step_opt_model_[2]),
  seed_ = 42,
  main_ = paste0(
    "MLMM GWAS Q-Q plot for ",
    trait_, " LS-means estimated phenotypes"
  ),
  output_path_ = paste0(output_gwas_graphics_path, trait_, "/ls_means/"),
  file_name_ = paste0(trait_, "_ls_means_mlmm_qq_plot.png")
)
qq_plot_done_


# get significant snps according to defined corrected Bonferroni threshold
mlmm_signif_snps_ls_means_ <- p_val_step_opt_model_[
  which(p_val_step_opt_model_ < threshold_)
]

if (length(mlmm_signif_snps_ls_means_) > 0) {
  # reformat significant snp(s) data
  mlmm_signif_snps_ls_means_ <- data.frame(
    "snp_name" = names(mlmm_signif_snps_ls_means_),
    "p_value" = as.numeric(mlmm_signif_snps_ls_means_),
    "chromosome_number" = map_df$chromosome[
      match(names(mlmm_signif_snps_ls_means_), map_df$snp.name)
    ],
    "position_in_bp" = map_df$position[
      match(names(mlmm_signif_snps_ls_means_), map_df$snp.name)
    ]
  )

  # write significant snp(s) infos
  fwrite(mlmm_signif_snps_ls_means_, file = paste0(
    paste0(output_gwas_graphics_path, trait_, "/ls_means/"),
    trait_, "_ls_means_mlmm_significant_marker_info.csv"
  ))

  # get the snps with highest significancy by chromosome
  mlmm_highest_signif_snps_ls_means_ <- mlmm_signif_snps_ls_means_ %>%
    group_by(chromosome_number) %>%
    filter(p_value == min(p_value)) %>%
    ungroup()

  # make boxplots for trait phenotypes according to SNP allele combinations (
  # i.e. homozygote and heterozygote ) of genotypes for SNP markers with highest
  # score (i.e. p-value or normalized variable importance)
  for (snp_ in mlmm_highest_signif_snps_ls_means_$snp_name) {
    # get marker score for boxplot
    score_ <- signif(
      as.numeric(
        mlmm_highest_signif_snps_ls_means_$p_value[
          match(snp_, mlmm_highest_signif_snps_ls_means_$snp_name)
        ]
      ), 3
    )
    df_ <- data.frame(
      "ls_means_y" = ls_means_y, "marker" = marker_matrix[, snp_]
    )
    r2_ <- summary(
      lm(as.formula("ls_means_y ~ 1 + marker"), data = df_)
    )$r.squared
    r2_ <- signif(as.numeric(r2_), 3)
    geno_box_plot_ <- genotypes_boxplot_(
      X = marker_df,
      Y = ls_means_y,
      ylab_ = paste0(trait_, " LS-means phenotypes"),
      xlab_ = paste0(
        "Genotype alleles combination on chromosome ",
        map_df$chromosome[match(snp_, map_df$snp.name)]
      ),
      markers = snp_,
      r2_ = r2_,
      score_ = score_,
      score_text_ = "p-value",
      output_path_ = paste0(output_gwas_graphics_path, trait_, "/ls_means/"),
      file_name_ = paste0(
        trait_, "_ls_means_mlmm_highest_signif_marker_", snp_, "_chromo_",
        map_df$chromosome[match(snp_, map_df$snp.name)], ".png"
      )
    )
  }
}

# ---- rf gwas for ls-means

# tune mtry before fitting ranger
if (tune_mtry_ranger_) {
  rf_obj_opt_mtry_ <- tune_mtry_ranger_rf(
    X = marker_matrix, Y = ls_means_y,
    mtry_grid_ = floor(ncol(marker_matrix) / c(10, 5, 3)),
    num_trees_ = 500,
    pkgs_to_export_ = "ranger"
  )
  opt_mtry <- rf_obj_opt_mtry_$opt_mtry
} else {
  # default to classical value used in rf implementations (which works pretty
  # well in many cases)
  opt_mtry <- floor(ncol(marker_matrix) / 3)
}

# fit rf using ranger
rf_obj_ls_means_ <- ranger(
  y = ls_means_y,
  x = marker_matrix,
  mtry = opt_mtry,
  num.trees = 500,
  importance = "permutation"
)

# get variable importance for markers
var_imp_ls_means_ <- importance(
  rf_obj_ls_means_
)

# normalize importance values using min-max normalization
var_imp_ls_means_ <- min_max_normalization(var_imp_ls_means_)

# apply a peak over threshold (POT) approach, using the generalized pareto
# distribution (GPD), to model the tail of the variable importance distribution
# which has extreme values
quant_basis_tresh_ <- quantile(var_imp_ls_means_, prob_basis_tresh_)
gpd_fit <- gpd(var_imp_ls_means_, threshold = quant_basis_tresh_)
xi_ <- gpd_fit$par.ests["xi"] # shape parameter
beta_ <- gpd_fit$par.ests["beta"] # scale parameter
quant_qgpd_tresh_ <- as.numeric(
  qgpd(prob_qgpd_tresh_, xi = xi_, beta = beta_)
)

if (length(var_imp_ls_means_) > 0) {
  # get the subset of variables (i.e. markers) detected as important
  sub_var_imp_ls_means_ <- var_imp_ls_means_[
    var_imp_ls_means_ > quant_basis_tresh_
  ]

  if (length(sub_var_imp_ls_means_) > 0) {
    sub_mat_ <- marker_matrix[, names(sub_var_imp_ls_means_)]

    # recompute normalized variable importance several times for subset of variables
    # detected as important
    list_sub_var_imp_ls_means <- vector("list", n_recompute_nvi)
    for (recomp_idx_ in 1:n_recompute_nvi) {
      list_sub_var_imp_ls_means[[recomp_idx_]] <- min_max_normalization(
        importance(
          ranger(
            y = ls_means_y,
            x = sub_mat_,
            mtry = floor(ncol(sub_mat_) / 3),
            num.trees = 500,
            importance = "permutation"
          )
        )
      )
    }
    mean_sub_var_imp_ls_means <- Reduce(`+`, list_sub_var_imp_ls_means) /
      length(list_sub_var_imp_ls_means)

    # reassign mean nvi into original vector of nvi for variables detected as
    # important
    var_imp_ls_means_[
      match(names(mean_sub_var_imp_ls_means), names(var_imp_ls_means_))
    ] <- mean_sub_var_imp_ls_means
  }
}

# create rf result data frame specific to manhattan ggplot
rf_df_results_ <- data.frame(
  "snp" = map_df$snp.name,
  "chr" = map_df$chromosome,
  "bp" = map_df$position,
  "score" = var_imp_ls_means_
)

# make a manhattan ggplot
manhattan_ggplot_plotly_done_ <- manhattan_ggplot_plotly_(
  df_ = rf_df_results_,
  trait_ = trait_,
  qgpd_threshold = quant_qgpd_tresh_,
  x_lab_ = "SNP index",
  y_lab_ = "Normalized variable importance (NVI)",
  title_ = paste0(
    "Random forest GWAS of ", trait_,
    " LS-means phenotypes using a 0.001 threshold for a generalized Pareto distribution (GPD)"
  ),
  title_size_ = 16,
  axis_title_size_ = 14,
  output_path_ = paste0(output_gwas_graphics_path, trait_, "/ls_means/"),
  file_name_ = paste0(trait_, "_ls_means_rf_manhattan_plot.png"),
  make_ggplotly_ = F
)
manhattan_ggplot_plotly_done_

# get significant snps according to generalized pareto distribution's
# 0.1% threshold
rf_signif_snps_ls_means_ <- var_imp_ls_means_[
  var_imp_ls_means_ > quant_qgpd_tresh_
]

if (length(rf_signif_snps_ls_means_) > 0) {
  # if any significant snps, convert their infos into a data frame containing
  # snp name, score and chromosome number
  rf_signif_snps_ls_means_df_ <- data.frame(
    "snp_name" = names(rf_signif_snps_ls_means_),
    "normalized_variable_importance" = rf_signif_snps_ls_means_,
    "chromosome_number" = map_df$chromosome[
      match(names(rf_signif_snps_ls_means_), map_df$snp.name)
    ],
    "position_in_bp" = map_df$position[
      match(names(rf_signif_snps_ls_means_), map_df$snp.name)
    ]
  )

  # get and write significant markers infos
  fwrite(rf_signif_snps_ls_means_df_, file = paste0(
    paste0(output_gwas_graphics_path, trait_, "/ls_means/"),
    trait_, "_ls_means_rf_significant_marker_info.csv"
  ))

  # get the snps with highest importance by chromosome
  rf_highest_signif_snps_ls_means_df_ <- rf_signif_snps_ls_means_df_ %>%
    group_by(chromosome_number) %>%
    filter(normalized_variable_importance ==
      max(normalized_variable_importance)) %>%
    ungroup()

  # make boxplots for trait phenotypes according to SNP allele combinations (
  # i.e. homozygote and heterozygote ) of genotypes for SNP markers with highest
  # score (i.e. p-value or normalized variable importance)
  for (snp_ in rf_highest_signif_snps_ls_means_df_$snp_name) {
    # get marker score for boxplot
    score_ <- signif(
      as.numeric(
        rf_highest_signif_snps_ls_means_df_$normalized_variable_importance[
          match(snp_, rf_highest_signif_snps_ls_means_df_$snp_name)
        ]
      ), 3
    )
    df_ <- data.frame(
      "ls_means_y" = ls_means_y, "marker" = marker_matrix[, snp_]
    )
    r2_ <- ranger(
      as.formula("ls_means_y ~ 1 + marker"),
      data = df_
    )$r.squared
    r2_ <- signif(as.numeric(r2_), 3)
    geno_box_plot_ <- genotypes_boxplot_(
      X = marker_df,
      Y = ls_means_y,
      ylab_ = paste0(trait_, " LS-means phenotypes"),
      xlab_ = paste0(
        "Genotype alleles combination on chromosome ",
        map_df$chromosome[match(snp_, map_df$snp.name)]
      ),
      markers = snp_,
      r2_ = r2_,
      score_ = score_,
      score_text_ = "NVI",
      output_path_ = paste0(output_gwas_graphics_path, trait_, "/ls_means/"),
      file_name_ = paste0(
        trait_, "_ls_means_rf_highest_signif_marker_", snp_, "_chromo_",
        map_df$chromosome[match(snp_, map_df$snp.name)], ".png"
      )
    )
  }
}

# get all significant snps identified using wiser and ls-means, and the two analysis
# models used (i.e. mlmm and rf)
all_signif_snps_ <- c()

if (!is.null(nrow(mlmm_signif_snps_wiser_)) &&
  nrow(mlmm_signif_snps_wiser_) > 0) {
  all_signif_snps_ <- c(
    all_signif_snps_,
    mlmm_signif_snps_wiser_$snp_name
  )
}
if (!is.null(nrow(mlmm_signif_snps_ls_means_)) &&
  nrow(mlmm_signif_snps_ls_means_) > 0) {
  all_signif_snps_ <- c(
    all_signif_snps_,
    mlmm_signif_snps_ls_means_$snp_name
  )
}
if (!is.null(rf_signif_snps_wiser_) &&
  length(rf_signif_snps_wiser_) > 0) {
  all_signif_snps_ <- c(
    all_signif_snps_,
    names(rf_signif_snps_wiser_)
  )
}
if (!is.null(rf_signif_snps_ls_means_) &&
  length(rf_signif_snps_ls_means_) > 0) {
  all_signif_snps_ <- c(
    all_signif_snps_,
    names(rf_signif_snps_ls_means_)
  )
}
# if any, get counts of significantly detected snps and assign their statistics according
# to the different models (i.e. mlmm and rf) and data used (i.e. ls-means and
# wiser)
if (length(all_signif_snps_) > 0) {
  count_signif_snps_ <- table(all_signif_snps_)

  stat_signif_snps_df_ <- data.frame(
    "snp_name" = names(count_signif_snps_),
    "count" = as.numeric(count_signif_snps_),
    "chromosome_number" = map_df$chromosome[
      match(names(count_signif_snps_), map_df$snp.name)
    ],
    "position_in_bp" = map_df$position[
      match(names(count_signif_snps_), map_df$snp.name)
    ],
    "ls_means_mlmm_neg_log10_p" = rep(NA, length(count_signif_snps_)),
    "ls_means_lm_r2" = rep(NA, length(count_signif_snps_)),
    "ls_means_rf_norm_var_import" = rep(NA, length(count_signif_snps_)),
    "ls_means_rf_r2" = rep(NA, length(count_signif_snps_)),
    "wiser_mlmm_neg_log10_p" = rep(NA, length(count_signif_snps_)),
    "wiser_lm_r2" = rep(NA, length(count_signif_snps_)),
    "wiser_rf_norm_var_import" = rep(NA, length(count_signif_snps_)),
    "wiser_rf_r2" = rep(NA, length(count_signif_snps_))
  )

  # assign ls-means p-values and normalized variable importance
  if (!is.null(nrow(mlmm_signif_snps_ls_means_)) &&
    nrow(mlmm_signif_snps_ls_means_) > 0) {
    idx_ <- match(
      mlmm_signif_snps_ls_means_$snp_name,
      names(count_signif_snps_)
    )
    stat_signif_snps_df_$ls_means_mlmm_neg_log10_p[
      idx_
    ] <- signif(
      -log10(as.numeric(mlmm_signif_snps_ls_means_$p_value)),
      3
    )
  }
  if (!is.null(rf_signif_snps_ls_means_) &&
    length(rf_signif_snps_ls_means_) > 0) {
    idx_ <- match(
      rf_signif_snps_ls_means_df_$snp_name,
      names(count_signif_snps_)
    )
    stat_signif_snps_df_$ls_means_rf_norm_var_import[
      idx_
    ] <- signif(
      rf_signif_snps_ls_means_df_$normalized_variable_importance,
      3
    )
  }

  # assign wiser p-values and normalized variable importance
  if (!is.null(nrow(mlmm_signif_snps_wiser_)) &&
    nrow(mlmm_signif_snps_wiser_) > 0) {
    idx_ <- match(
      mlmm_signif_snps_wiser_$snp_name,
      names(count_signif_snps_)
    )
    stat_signif_snps_df_$wiser_mlmm_neg_log10_p[
      idx_
    ] <- signif(
      -log10(as.numeric(mlmm_signif_snps_wiser_$p_value)),
      3
    )
  }
  if (!is.null(rf_signif_snps_wiser_) &&
    length(rf_signif_snps_wiser_) > 0) {
    idx_ <- match(
      rf_signif_snps_wiser_df_$snp_name,
      names(count_signif_snps_)
    )
    stat_signif_snps_df_$wiser_rf_norm_var_import[
      idx_
    ] <- signif(
      rf_signif_snps_wiser_df_$normalized_variable_importance,
      3
    )
  }

  # compute the coefficient of determination R2 for all the significant snps
  for (i in 1:nrow(stat_signif_snps_df_)) {
    if (!is.na(stat_signif_snps_df_$ls_means_mlmm_neg_log10_p[i])) {
      snp_ <- stat_signif_snps_df_$snp_name[i]
      df_ <- data.frame(
        "ls_means_y" = ls_means_y, "marker" = marker_matrix[, snp_]
      )
      r2_ <- summary(
        lm(as.formula("ls_means_y ~ 1 + marker"), data = df_)
      )$r.squared
      stat_signif_snps_df_$ls_means_lm_r2[i] <- signif(r2_, 3)
    }

    if (!is.na(stat_signif_snps_df_$ls_means_rf_norm_var_import[i])) {
      snp_ <- stat_signif_snps_df_$snp_name[i]
      df_ <- data.frame(
        "ls_means_y" = ls_means_y, "marker" = marker_matrix[, snp_]
      )
      r2_ <- ranger(
        as.formula("ls_means_y ~ 1 + marker"),
        data = df_
      )$r.squared
      stat_signif_snps_df_$ls_means_rf_r2[i] <- signif(r2_, 3)
    }

    if (!is.na(stat_signif_snps_df_$wiser_mlmm_neg_log10_p[i])) {
      snp_ <- stat_signif_snps_df_$snp_name[i]
      df_ <- data.frame(
        "wiser_y" = wiser_y, "marker" = marker_matrix[, snp_]
      )
      r2_ <- summary(
        lm(as.formula("wiser_y ~ 1 + marker"), data = df_)
      )$r.squared
      stat_signif_snps_df_$wiser_lm_r2[i] <- signif(r2_, 3)
    }

    if (!is.na(stat_signif_snps_df_$wiser_rf_norm_var_import[i])) {
      snp_ <- stat_signif_snps_df_$snp_name[i]
      df_ <- data.frame(
        "wiser_y" = wiser_y, "marker" = marker_matrix[, snp_]
      )
      r2_ <- ranger(
        as.formula("wiser_y ~ 1 + marker"),
        data = df_
      )$r.squared
      stat_signif_snps_df_$wiser_rf_r2[i] <- signif(r2_, 3)
    }
  }

  # reorder data frame according to chromosome number and position in bp
  stat_signif_snps_df_ <- stat_signif_snps_df_ %>%
    arrange(chromosome_number, position_in_bp)

  # write complete statistics for significant snps detected for trait
  fwrite(stat_signif_snps_df_, file = paste0(
    output_gwas_results_path,
    trait_, "_gwas_results.csv"
  ))
}
