# script meant to perform genomic prediction and analyses for refpop
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
use_tensorflow_or_umap <- F
if (use_tensorflow_or_umap) {
  # leave tensorflow and keras for later use
  library(tensorflow)
  library(keras3)
  library(umap)
  tensorflow::tf$random$set_seed(0)
  py_module_available("keras") # must return TRUE
  py_module_available("tensorflow") # must return TRUE
  py_discover_config("keras") # more info on the python env, tf and keras
}
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
options(future.globals.maxSize = 6 * 1024^3)
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

# set path for identified environment and phenotypic outliers
outlier_dir_path <- "../../results/phenotype_outlier_detection/"

# set path for wiser phenotypes estimated using whitening
wiser_pheno_dir_path <- "../../data/phenotype_data/wiser_phenotypes_and_estimates/"

# output result path for genotype graphics
output_pred_results_path <- "../../results/genomic_prediction/"
output_pred_graphics_path <- "../../results/graphics/genomic_prediction_graphics/"

# define kernels for wiser
kernels_ <- c("linear", "identity")

# define traits_ for wiser
traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Weight_sample"
)
skipped_vars_ <- c("Sample_size")

# get kernel and trait arguments
args <- commandArgs(trailingOnly = TRUE)
kernel_num <- as.integer(args[1])
trait_num <- as.integer(args[2])

# kernel type, i.e. "linear" or "identity" for genomic covariance matrix
# (i.e. Gram matrix). NB. "identity" is not recommended due to hypothesis of
# independence between genotypes which is highly unlikely
kernel_ <- kernels_[kernel_num]
print(paste0("kernel: ", kernel_))

# define trait_
trait_ <- traits_[trait_num]
print(paste0("trait: ", trait_))

# define shift seed value by
mult_seed_by_ <- 100

# set k for K-folds cv
k_folds_ <- 5

# define number of shuffles
n_shuff_ <- 20

# define number of snp for sampling
# according to a uniform distribution
snp_sample_size_ <- 50e3

# get raw and adjusted ls-means phenotype and genotype data
raw_pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path, "phenotype_raw_data_no_outliers.csv"
)))

pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "adjusted_ls_means_phenotypes.csv"
)))

omic_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_data.csv"
)))

# get environments based on (estimable) heritabilities which are not identified
# as outliers (i.e. kept environments after SpATS(), cf.
# refpop_1_spat_hetero_correct_per_env_trait_and_h2_estim.R) and for which trait
# has available data
sel_env_ <- as.data.frame(
  fread(
    paste0(
      outlier_dir_path,
      "envir_per_trait_retained_based_on_non_miss_data_and_estimable_h2_non_outliers_for_adj_phenotypes.csv"
    )
  )
)
sel_env_ <- str_split(sel_env_[sel_env_$trait == trait_, ], pattern = ", ")[[2]]
raw_pheno_df <- raw_pheno_df[raw_pheno_df$Envir %in% sel_env_, ]

# remove monomorphic markers
omic_df <- remove_monomorphic_markers(omic_df)
monomorphic_markers_list_ <- omic_df$monomorphic_markers
omic_df <- omic_df$filtered_df

# sample markers according to a uniform distribution
set.seed(123)
idx_snp_sample_size_ <- sample(2:ncol(omic_df),
  size = snp_sample_size_, replace = F
)
omic_df <- omic_df[, c(
  match("Genotype", colnames(omic_df)),
  idx_snp_sample_size_
)]

# merge pheno_df and omic_df for integrity of analyses
merged_df <- merge(pheno_df, omic_df, by = "Genotype")
pheno_df <- merged_df[, c("Genotype", traits_, skipped_vars_)]
omic_df <- merged_df[, -match(
  c("Genotype", traits_, skipped_vars_),
  colnames(merged_df)
)]

# remove na for analyzed trait_ and corresponding rows for marker data
idx_na_trait_ <- which(is.na(pheno_df[, trait_]))
if (length(idx_na_trait_) > 0) {
  pheno_df <- pheno_df[-idx_na_trait_, ]
  omic_df <- omic_df[-idx_na_trait_, ]
}
rownames(omic_df) <- pheno_df$Genotype

# get number of genotypes
n <- length(pheno_df$Genotype)

# get optimal whitening method and regularization parameter using k-folds CV
opt_white_reg_par <- optimize_whitening_and_regularization(
  omic_df, raw_pheno_df, trait_,
  fixed_effects_vars = c(
    "Envir", "Country", "Year",
    "Row", "Position", "Management"
  ),
  fixed_effects_vars_computed_as_factor = c(
    "Envir", "Country", "Year",
    "Row", "Position", "Management"
  ),
  site_var = "Country",
  fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
  random_effects_vars = "Genotype",
  whitening_method_grid = c("ZCA-cor", "PCA-cor", "Cholesky"),
  alpha_frob_grid = c(0.01, 0.1),
  k_folds = 5
)
print(opt_white_reg_par)
opt_whitening_method_ <- as.character(opt_white_reg_par$opt_whitening_method)
opt_alpha_frob_ <- opt_white_reg_par$opt_alpha_frob
rm(opt_white_reg_par)

# compute wiser phenotypes,
# since computations are long, save results for later use

# if exist, read corrected phenotype data for trait associated to the
# defined kernel
if (file.exists(paste0(
  wiser_pheno_dir_path,
  "wiser_phenotypes_and_estimates_", kernel_,
  "_kernel_", trait_, ".csv"
))) {
  # load corrected phenotypes if file exists
  wiser_pheno_df <- as.data.frame(fread(paste0(
    wiser_pheno_dir_path,
    "wiser_phenotypes_and_estimates_", kernel_,
    "_kernel_", trait_, ".csv"
  )))
  v_hat <- wiser_pheno_df$v_hat
} else {
  # estimate wiser phenotype
  start_time_ <- Sys.time()
  wiser_obj <- estimate_wiser_phenotype(omic_df, raw_pheno_df, trait_,
    fixed_effects_vars = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    fixed_effects_vars_computed_as_factor = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    site_var = "Country",
    fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
    random_effects_vars = "Genotype",
    kernel_type = kernel_,
    whitening_method = opt_whitening_method_,
    alpha_frob_ = opt_alpha_frob_
  )
  end_time_ <- Sys.time()
  time_taken_ <- end_time_ - start_time_
  print(paste0(
    "Execution time for wiser computation of ",
    trait_, " with ", kernel_, " kernel : ",
    signif(time_taken_, 3)
  ))

  # get estimated wiser phenotype
  v_hat <- wiser_obj$wiser_phenotypes$v_hat

  # save wiser phenotype for kernel and trait
  fwrite(wiser_obj$wiser_phenotypes, paste0(
    wiser_pheno_dir_path,
    "wiser_phenotypes_and_estimates_", kernel_,
    "_kernel_", trait_, ".csv"
  ), row.names = F, col.names = T)

  # save wiser object for kernel and trait
  saveRDS(wiser_obj, paste0(
    wiser_pheno_dir_path,
    "wiser_obj_", kernel_,
    "_kernel_", trait_
  ))

  # verify the property that x_mat_tilde is indeed orthogonal
  # (i.e. decorrelated) to xi_hat
  x_mat_tilde_mult_xi_test <- t(
    wiser_obj$wiser_x_mat_tilde
  ) %*% wiser_obj$wiser_xi_hat

  # for readability purposes print only the first 5 elements,
  # the mean and median ofx_mat_tilde_mult_xi_test
  print(paste0(
    "x_mat_tilde_mult_xi_test[1:10]: ",
    paste0(head(x_mat_tilde_mult_xi_test, 5), collapse = ", ")
  ))
  print(paste0(
    "mean of x_mat_tilde_mult_xi_test: ",
    mean(x_mat_tilde_mult_xi_test)
  ))
  print(paste0(
    "median of x_mat_tilde_mult_xi_test: ",
    median(x_mat_tilde_mult_xi_test)
  ))
}

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
    if (kernel_ == "linear") {
      fold_result <- c(
        "RF_wiser_linear" = NA,
        "SVR_wiser_linear" = NA,
        "GBLUP_wiser_linear" = NA,
        "RKHS_wiser_linear" = NA,
        "LASSO_wiser_linear" = NA,
        "RF_ls_means" = NA,
        "SVR_ls_means" = NA,
        "GBLUP_ls_means" = NA,
        "RKHS_ls_means" = NA,
        "LASSO_ls_means" = NA
      )
    } else {
      fold_result <- c(
        "RF_wiser_identity" = NA,
        "SVR_wiser_identity" = NA,
        "GBLUP_wiser_identity" = NA,
        "RKHS_wiser_identity" = NA,
        "LASSO_wiser_identity" = NA,
        "RF_ls_means" = NA,
        "SVR_ls_means" = NA,
        "GBLUP_ls_means" = NA,
        "RKHS_ls_means" = NA,
        "LASSO_ls_means" = NA
      )
    }

    # training and prediction based on computed phenotypes, i.e. v_hat

    # train and predict with Random Forest
    rf_model <- ranger(
      y = v_hat[idx_train],
      x = omic_df[idx_train, ],
      mtry = ncol(omic_df) / 3,
      num.trees = 1000
    )
    f_hat_val_rf <- predict(
      rf_model,
      omic_df[idx_val, ]
    )
    fold_result[1] <- cor(
      f_hat_val_rf$predictions,
      v_hat[idx_val]
    )

    # train and predict with SVR
    # a correct value for c_par according to Cherkassy and Ma (2004).
    # Neural networks 17, 113-126 is defined as follows
    c_par <- max(
      abs(mean(v_hat[idx_train])
      + 3 * sd(v_hat[idx_train])),
      abs(mean(v_hat[idx_train])
      - 3 * sd(v_hat[idx_train]))
    )
    gaussian_svr_model <- ksvm(
      x = as.matrix(omic_df[idx_train, ]),
      y = v_hat[idx_train],
      scaled = FALSE, type = "eps-svr",
      kernel = "rbfdot",
      kpar = "automatic", C = c_par, epsilon = 0.1
    )
    f_hat_val_gaussian_svr <- predict(
      gaussian_svr_model,
      as.matrix(omic_df[idx_val, ])
    )
    fold_result[2] <- cor(
      f_hat_val_gaussian_svr,
      v_hat[idx_val]
    )

    # train and predict with GBLUP (linear kernel krmm)
    linear_krmm_model <- krmm(
      Y = v_hat[idx_train],
      Matrix_covariates = omic_df[idx_train, ],
      method = "GBLUP"
    )
    f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
      Matrix_covariates = omic_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result[3] <- cor(
      f_hat_val_linear_krmm,
      v_hat[idx_val]
    )
    # train and predict with RKHS (non-linear Gaussian kernel krmm)
    gaussian_krmm_model <- krmm(
      Y = v_hat[idx_train],
      Matrix_covariates = omic_df[idx_train, ],
      method = "RKHS", kernel = "Gaussian",
      rate_decay_kernel = 0.1
    )
    f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
      Matrix_covariates = omic_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result[4] <- cor(
      f_hat_val_gaussian_krmm,
      v_hat[idx_val]
    )

    # train and predict with LASSO
    cv_fit_lasso_model <- cv.glmnet(
      intercept = TRUE, y = v_hat[idx_train],
      x = as.matrix(omic_df[idx_train, ]),
      type.measure = "mse", alpha = 1.0, nfold = 10,
      parallel = TRUE
    )
    f_hat_val_lasso <- predict(cv_fit_lasso_model,
      newx = as.matrix(omic_df[idx_val, ]),
      s = "lambda.min"
    )
    fold_result[5] <- cor(
      f_hat_val_lasso,
      v_hat[idx_val]
    )

    # training and prediction based on ls-means

    # train and predict with Random Forest
    rf_model <- ranger(
      y = pheno_df[idx_train, trait_],
      x = omic_df[idx_train, ],
      mtry = ncol(omic_df) / 3,
      num.trees = 1000
    )
    f_hat_val_rf <- predict(
      rf_model,
      omic_df[idx_val, ]
    )
    fold_result[6] <- cor(
      f_hat_val_rf$predictions,
      pheno_df[idx_val, trait_]
    )

    # train and predict with SVR
    # a correct value for c_par according to Cherkassy and Ma (2004).
    # Neural networks 17, 113-126 is defined as follows
    c_par <- max(
      abs(mean(pheno_df[idx_train, trait_])
      + 3 * sd(pheno_df[idx_train, trait_])),
      abs(mean(pheno_df[idx_train, trait_])
      - 3 * sd(pheno_df[idx_train, trait_]))
    )
    gaussian_svr_model <- ksvm(
      x = as.matrix(omic_df[idx_train, ]),
      y = pheno_df[idx_train, trait_],
      scaled = FALSE, type = "eps-svr",
      kernel = "rbfdot",
      kpar = "automatic", C = c_par, epsilon = 0.1
    )
    f_hat_val_gaussian_svr <- predict(
      gaussian_svr_model,
      as.matrix(omic_df[idx_val, ])
    )
    fold_result[7] <- cor(
      f_hat_val_gaussian_svr,
      pheno_df[idx_val, trait_]
    )

    # train and predict with GBLUP (linear kernel krmm)
    linear_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = omic_df[idx_train, ],
      method = "GBLUP"
    )
    f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
      Matrix_covariates = omic_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result[8] <- cor(
      f_hat_val_linear_krmm,
      pheno_df[idx_val, trait_]
    )

    # train and predict with RKHS (non-linear Gaussian kernel krmm)
    gaussian_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = omic_df[idx_train, ],
      method = "RKHS", kernel = "Gaussian",
      rate_decay_kernel = 0.1
    )
    f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
      Matrix_covariates = omic_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result[9] <- cor(
      f_hat_val_gaussian_krmm,
      pheno_df[idx_val, trait_]
    )

    # train and predict with LASSO
    cv_fit_lasso_model <- cv.glmnet(
      intercept = TRUE, y = pheno_df[idx_train, trait_],
      x = as.matrix(omic_df[idx_train, ]),
      type.measure = "mse", alpha = 1.0, nfold = 10,
      parallel = TRUE
    )
    f_hat_val_lasso <- predict(cv_fit_lasso_model,
      newx = as.matrix(omic_df[idx_val, ]),
      s = "lambda.min"
    )
    fold_result[10] <- cor(
      f_hat_val_lasso,
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
    yaxis = list(
      title = "Predictive ability (PA)",
      range = c(0, 1)
    ),
    legend = list(title = list(text = "Prediction method"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path, trait_, "/wiser_phenotype_predictive_ability_",
  trait_, "_", kernel_, "_kernel_", snp_sample_size_, "_SNP_",
  k_folds_, "_folds_CV.html"
))

# add stats and save predictive ability results
df_result_[, 1:10] <- signif(apply(df_result_[, 1:10], 2, as.numeric), 2)
rownames(df_result_) <- paste0("pa_scenario_", 1:nrow(df_result_))

df_stat <- as.data.frame(rbind(
  apply(df_result_[, 1:10], 2, mean),
  apply(df_result_[, 1:10], 2, sd)
))
df_stat <- signif(apply(df_stat, 2, as.numeric), 2)
rownames(df_stat) <- c("pa_mean", "pa_sd")
df_stat <- as.data.frame(df_stat)

df_result_ <- rbind(df_result_, df_stat)

fwrite(df_result_,
  file = paste0(
    output_pred_results_path,
    "wiser_phenotype_genomic_pred_results_", ncol(omic_df), "_SNP_",
    trait_, "_", kernel_, "_kernel_", k_folds_, "_folds_CV.csv"
  ), row.names = T
)

print(df_result_)
