# script meant to perform genomic prediction and analyses for refpop
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(MASS)
library(data.table)
library(stringr)
library(lme4)
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
library(reticulate)
library(devtools)
if ("refpop_env" %in% conda_list()$name) {
  use_condaenv("refpop_env")
}
# install other requirements from github if necessary
install_other_requirements <- F
if (install_other_requirements) {
  # reticulate::install_miniconda()
  conda_create("refpop_env")
  use_condaenv("refpop_env")
  devtools::install_github("ljacquin/KRMM")
  devtools::install_github("rstudio/tensorflow")
  library(tensorflow)
  install_tensorflow(envname = "refpop_env")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
  install.packages("umap")
}
library(KRMM)
library(kernlab)
library(glmnet)
library(ranger)
library(tensorflow)
library(keras3)
library(umap)
tensorflow::tf$random$set_seed(0)
py_module_available("keras") # must return TRUE
py_module_available("tensorflow") # must return TRUE
py_discover_config("keras") # more info on the python env, tf and keras

# define computation mode, i.e. local or cluster
computation_mode <- "local"

# if comutations are local in rstudio, detect and set script path
# automatically using rstudioapi
if (identical(computation_mode, "local")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}

# source functions
source("../functions.R")

# set options to increase memory and suppress warnings
options(expressions = 5e5)
options(warn = -1)

# define function(s) and package(s) to export for parallelization
pkgs_to_export_ <- c(
  "ranger",
  "kernlab",
  "KRMM",
  "glmnet",
  "foreach",
  "cvTools"
)
# set input paths
geno_dir_path <- "../../data/genotype_data/"
pheno_dir_path <- "../../data/phenotype_data/"

# output result path for genotype graphics
output_pred_results_path <- "../../results/genomic_prediction/"
output_pred_graphics_path <- "../../results/graphics/genomic_prediction_graphics/"

# define selected_traits_
selected_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Weight_sample", "Sample_size"
)

# define trait_
trait_ <- "Scab"

# define shift seed value by
mult_seed_by_ <- 100

# set k for K-folds cv
k_folds_ <- 5

# define number of shuffles
n_shuff_ <- 1

# color palette for families (28 counts)
color_palette_family <- c(
  "black",
  "lightblue",
  colorRampPalette(c("blue", "deepskyblue"))(10),
  colorRampPalette(c("orange", "orange2"))(2),
  colorRampPalette(c("aquamarine"))(1),
  colorRampPalette(c("magenta", "magenta2"))(2),
  colorRampPalette(c("firebrick3", "firebrick4"))(2),
  colorRampPalette(c("green", "darkgreen"))(5),
  colorRampPalette(c("darkorchid4"))(1),
  colorRampPalette(c("gold1", "gold2"))(2),
  colorRampPalette(c("lightcoral"))(1)
)

# color palette for origins (11 counts)
color_palette_origin <- c(
  "magenta",
  "lightblue",
  "blue",
  "orange",
  "aquamarine",
  "red",
  "green",
  "darkorchid4",
  "gold2",
  "lightcoral",
  "yellow"
)

# get spats and adjusted ls-means phenotype and genotype data,
# and family and origin info

spats_pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path, "spats_per_env_adjusted_phenotypes/",
  trait_, "_spats_adjusted_phenotypes_long_format.csv"
)))

pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "adjusted_ls_means_phenotypes.csv"
)))

geno_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_data.csv"
)))

geno_fam_orig_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_family_origin_information.csv"
)))

set.seed(1)
pheno_df <- pheno_df[sample(1:100, replace = F), ]

# remove monomorphic markers
geno_df <- remove_monomorphic_markers(geno_df)
monomorphic_markers_list_ <- geno_df$monomorphic_markers
geno_df <- geno_df$filtered_df

# set genotypes in the same order between pheno and geno data
# and sample markers according to a uniform distribution
set.seed(123)
snp_sample_size_ <- 50e3

idx_snp_sample_size_ <- sample(2:ncol(geno_df),
  size = snp_sample_size_, replace = F
)

geno_df <- geno_df[, c(
  match("Genotype", colnames(geno_df)),
  idx_snp_sample_size_
)]

# merge pheno_df and geno_df for integrity of analyses and slice the merged df
merged_df <- merge(pheno_df, geno_df, by = "Genotype")
pheno_df <- merged_df[, c("Genotype", selected_traits_)]
geno_df <- merged_df[, -match(
  c("Genotype", selected_traits_),
  colnames(merged_df)
)]

# remove na for analyzed trait_ and corresponding rows for marker data
idx_na_trait_ <- which(is.na(pheno_df[, trait_]))
if (length(idx_na_trait_) > 0) {
  pheno_df <- pheno_df[-idx_na_trait_, ]
  geno_df <- geno_df[-idx_na_trait_, ]
}
rownames(geno_df) <- pheno_df$Genotype

# get number of genotypes
n <- nrow(geno_df)

# compute the general mixed model (Rao & Keffe) Y = X*Beta + E = X*Beta + Zu + E
# and whiten the residuals of E where Var(E) = Σ = sigma2_u*ZKZ' sigma2_e*In',
# i.e. compute  Ytilde = Xtilde*Beta + Etilde where tilde terms are obtained
# by multiplying by L^{-1} after Cholesky decomp. of Σ = LL'
ebv_post_white_ <- compute_ebv_post_whitening(
  trait_, spats_pheno_df, geno_df,
  sigma2u = 1, sigma2e = 1
)
fitted_spats_pheno_df <- ebv_post_white_$fitted_spats_pheno_df
x_mat <- ebv_post_white_$x_mat
z_mat <- ebv_post_white_$z_mat
k_mat <- ebv_post_white_$k_mat
u_hat <- ebv_post_white_$u_hat
beta_hat <- ebv_post_white_$beta_hat
beta_hat
plot(density(u_hat))

sigma2_upper_bound <- var(fitted_spats_pheno_df[, trait_])
var_comp_estim_abc_ <- abc_variance_component_estimation(
  trait_, fitted_spats_pheno_df,
  x_mat, z_mat, k_mat,
  beta_hat,
  prior_sig2_u = c(1e-2, sigma2_upper_bound),
  prior_sig2_e = c(1e-2, sigma2_upper_bound),
  n_sim_ = 100,
  quantile_threshold = 0.05
)
var_comp_estim_abc_

ebv_post_white_ <- compute_ebv_post_whitening(
  trait_, spats_pheno_df, geno_df,
  sigma2u = var_comp_estim_abc_$sigma2_u_hat_mean,
  sigma2e = var_comp_estim_abc_$sigma2_e_hat_mean
)
u_hat <- ebv_post_white_$u_hat
plot(density(u_hat))

# register parallel backend
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# create folds for k-fold cross-validation
Folds <- cvFolds(n, K = k_folds_, type = "consecutive")

df_result_ <- foreach(
  shuff_ = 1:n_shuff_,
  .packages = pkgs_to_export_,
  .combine = rbind
) %dopar% {
  # set seed, define a new set of indices,
  set.seed(shuff_ * mult_seed_by_)
  idx_ <- sample(1:n, size = n, replace = FALSE)
  fold_results <- foreach(
    fold_ = 1:k_folds_,
    .packages = pkgs_to_export_,
    .combine = rbind
  ) %dopar% {
    # define indices for validation and training based on k-folds cv
    idx_val_fold <- which(Folds$which == fold_)
    idx_val <- idx_[idx_val_fold]
    idx_train <- idx_[-idx_val_fold]

    # initialize vector of results for the current fold
    fold_result <- c(
      "GBLUP_ebv_post_whitening" = NA, "RKHS_ebv_post_whitening" = NA,
      "GBLUP_ls_means" = NA, "RKHS_ls_means" = NA
    )
    # predict estimated breeding values (ebv), using ebv during training
    # train and predict with GBLUP (linear kernel krmm)
    linear_krmm_model <- krmm(
      Y = u_hat[idx_train],
      Matrix_covariates = geno_df[idx_train, ],
      method = "GBLUP"
    )
    f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
      Matrix_covariates = geno_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result["GBLUP_ebv_post_whitening"] <- cor(
      f_hat_val_linear_krmm,
      u_hat[idx_val]
    )
    # train and predict with RKHS (non-linear Gaussian kernel krmm)
    gaussian_krmm_model <- krmm(
      Y = u_hat[idx_train],
      Matrix_covariates = geno_df[idx_train, ],
      method = "RKHS", kernel = "Gaussian",
      rate_decay_kernel = 0.1
    )
    f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
      Matrix_covariates = geno_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result["RKHS_ebv_post_whitening"] <- cor(
      f_hat_val_gaussian_krmm,
      u_hat[idx_val]
    )
    # predict ls-means, using ls-means during training
    # train and predict with GBLUP (linear kernel krmm)
    linear_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = geno_df[idx_train, ],
      method = "GBLUP"
    )
    f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
      Matrix_covariates = geno_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result["GBLUP_ls_means"] <- cor(
      f_hat_val_linear_krmm,
      pheno_df[idx_val, trait_]
    )
    # train and predict with RKHS (non-linear Gaussian kernel krmm)
    gaussian_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = geno_df[idx_train, ],
      method = "RKHS", kernel = "Gaussian",
      rate_decay_kernel = 0.1
    )
    f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
      Matrix_covariates = geno_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result["RKHS_ls_means"] <- cor(
      f_hat_val_gaussian_krmm,
      pheno_df[idx_val, trait_]
    )
    fold_result
  }
  fold_results
}
# stop the parallel backend
stopCluster(cl)
registerDoSEQ()

# create directory for trait_ graphics if it does not exist
if (!dir.exists(paste0(output_pred_graphics_path, trait_, "/"))) {
  dir.create(paste0(output_pred_graphics_path, trait_, "/"))
}

# get methods names
df_result_ <- as.data.frame(apply(df_result_, 2, as.numeric))
method_names <- colnames(df_result_)

# initialize plot_ly boxplot graphic
boxplots_pa_ <- plot_ly()

# add boxplots
for (method_ in method_names) {
  boxplots_pa_ <- add_boxplot(
    boxplots_pa_,
    y = df_result_[[method_]],
    name = method_,
    boxpoints = "all",
    jitter = 0.3,
    pointpos = -1.8
  )
}

# add layout
boxplots_pa_ <- boxplots_pa_ %>%
  layout(
    title = paste0(
      "Genomic prediction PA distributions of methods for ",
      trait_, ", based on ", snp_sample_size_, " SNPs across ",
      n_shuff_, " shuffling scenarios for ", k_folds_, "-folds CV"
    ),
    yaxis = list(title = "Predictive ability (PA)", range = c(0, 1)),
    legend = list(title = list(text = "Prediction method"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path, trait_, "/post_whitening_predictive_ability_",
  trait_, "_", snp_sample_size_, "_SNP_", k_folds_, "_folds_CV_single_step.html"
))

# add stats and save predictive ability results
df_result_[, 1:4] <- signif(apply(df_result_[, 1:4], 2, as.numeric), 2)
rownames(df_result_) <- paste0("pa_scenario_", 1:nrow(df_result_))

df_stat <- as.data.frame(rbind(
  apply(df_result_[, 1:4], 2, mean),
  apply(df_result_[, 1:4], 2, sd)
))
df_stat <- signif(apply(df_stat, 2, as.numeric), 2)
rownames(df_stat) <- c("pa_mean", "pa_sd")
df_stat <- as.data.frame(df_stat)

df_result_ <- rbind(df_result_, df_stat)

fwrite(df_result_,
  file = paste0(
    output_pred_results_path,
    "post_whitening_genomic_pred_results_", ncol(geno_df), "_SNP_",
    trait_, "_", k_folds_, "_folds_CV_single_step.csv"
  ), row.names = T
)
