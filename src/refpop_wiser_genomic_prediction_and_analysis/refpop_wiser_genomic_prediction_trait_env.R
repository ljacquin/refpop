# script meant to perform wiser and ls-means genomic prediction and analyses for
# a refpop environment
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
library(emmeans)

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

# set color gradients and color vector for predictive abilities
blue_gradient <- c("#90B3E0", "#3D9BC5", "#005AB5", "#00407A", "#002A66")
yellow_orange_gradient <- colorRampPalette(c("#FFEA00", "#FF7A00"))(5)
pa_colors_ <- c(blue_gradient, yellow_orange_gradient)

# set color vector for computed genomic heritabilities (h2)
h2_colors_ <- c(blue_gradient[3], yellow_orange_gradient[3])

# define number of cores
nb_cores_ <- 12

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

# set path for wiser phenotypes estimated using whitening
wiser_pheno_dir_path <- "../../data/phenotype_data/wiser_phenotype_estimates_env/"

# output result path for genotype graphics
output_pred_results_path <- "../../results/genomic_prediction_env/"
output_pred_graphics_path <- "../../results/graphics/genomic_prediction_graphics_env/"

# define kernels for wiser
kernels_ <- c("linear", "identity")

# define traits_
traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Weight_sample"
)

# define selected year and management type for analysis
sel_year_management <- "2023_1"

# get kernel and trait arguments
args <- commandArgs(trailingOnly = T)
env_num <- as.integer(args[1])
kernel_num <- as.integer(args[2])
trait_num <- as.integer(args[3])

# # kernel type, i.e. "linear" or "identity" for genomic covariance matrix
# # (i.e. Gram matrix). NB. "identity" is not recommended due to hypothesis of
# # independence between genotypes which is highly unlikely
kernel_ <- kernels_[kernel_num]

# # define trait_
trait_ <- traits_[trait_num]

# define shift seed value by
mult_seed_by_ <- 100

# set k for K-folds cv
k_folds_ <- 5

# define number of shuffles
n_shuff_ <- 20

# define number of snp for sampling
# according to a uniform distribution
snp_sample_size_ <- 50e3

# get phenotype data
raw_pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path, "phenotype_data.csv"
)))

# get fist environment for selected year
list_env_sel_year_manage <- raw_pheno_df$Envir[
  str_detect(raw_pheno_df$Envir, pattern = sel_year_management)
]
env_ <- list_env_sel_year_manage[env_num]

# print combination of kernel, trait and analyzed environment
print(paste0("kernel: ", kernel_))
print(paste0("trait: ", trait_))
print(paste0("environment: ", env_))

# get phenotype for selected environments and genotype data
raw_pheno_df <- raw_pheno_df[raw_pheno_df$Envir == env_, ]
raw_pheno_df <- raw_pheno_df[!is.na(raw_pheno_df[, trait_]), ]

omic_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_data.csv"
)))

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

# compute trait ls-means per genotype for selected environnement (i.e. env_)
lm_ <- lm(
  as.formula(paste0(
    trait_,
    " ~ Genotype + Row + Position"
  )),
  data = raw_pheno_df
)
ls_means_df <- as.data.frame(
  emmeans(lm_, ~Genotype)
)[, c("Genotype", "emmean")]
colnames(ls_means_df)[2] <- trait_

# merge pheno_df and omic_df for integrity of analyses and slice the merged df
merged_df <- merge(ls_means_df, omic_df, by = "Genotype")
pheno_df <- merged_df[, c("Genotype", trait_)]
omic_df <- merged_df[, -match(
  c("Genotype", trait_),
  colnames(merged_df)
)]

# assign genotype names to omic_df rows and get number of genotypes
rownames(omic_df) <- merged_df$Genotype
n <- length(merged_df$Genotype)
rm(merged_df)

# compute wiser phenotype, corrected for fixed effects which takes into account
# genetic covariance between genotypes, using a whitening algorithm and
# approximate bayesian computation (ABC)
# since computations are long, save results for later use

# if exist, read corrected phenotype data for trait associated to the
# defined kernel
if (file.exists(paste0(
  wiser_pheno_dir_path,
  "wiser_phenotype_estimates_env_", kernel_,
  "_kernel_", trait_, "_", env_, ".csv"
))) {
  # load corrected phenotypes if file exists
  wiser_obj <- readRDS(paste0(
    wiser_pheno_dir_path,
    "wiser_obj_", kernel_,
    "_kernel_", trait_, "_", env_
  ))
  wiser_pheno_df <- wiser_obj$wiser_phenotypes
  omic_df <- wiser_obj$wiser_omic_data
  rm(wiser_obj)
} else {
  # get optimal whitening method and regularization parameter using k-folds CV
  opt_white_reg_par <- optimize_whitening_and_regularization(
    omic_df, raw_pheno_df, trait_,
    fixed_effects_vars = c(
      "Row", "Position"
    ),
    fixed_effects_vars_computed_as_factor = c(
      "Row", "Position"
    ),
    envir_var = NULL,
    fixed_effects_vars_computed_as_factor_by_envir = c("Row", "Position"),
    random_effects_vars = "Genotype",
    whitening_method_grid = c("ZCA-cor", "PCA-cor", "Cholesky"),
    alpha_grid = c(0.01, 0.1),
    k_folds = 5
  )
  print(opt_white_reg_par)
  opt_whitening_method_ <- as.character(opt_white_reg_par$opt_whitening_method)
  opt_alpha_ <- as.numeric(opt_white_reg_par$opt_alpha_)
  rm(opt_white_reg_par)

  # estimate wiser phenotype
  start_time_ <- Sys.time()
  wiser_obj <- estimate_wiser_phenotype(omic_df, raw_pheno_df, trait_,
    fixed_effects_vars = c(
      "Row", "Position"
    ),
    fixed_effects_vars_computed_as_factor = c(
      "Row", "Position"
    ),
    envir_var = NULL,
    fixed_effects_vars_computed_as_factor_by_envir = c("Row", "Position"),
    random_effects_vars = "Genotype",
    kernel_type = kernel_,
    whitening_method = opt_whitening_method_,
    alpha_ = opt_alpha_
  )
  end_time_ <- Sys.time()
  time_taken_ <- end_time_ - start_time_
  print(paste0(
    "Execution time for wiser computation of ",
    trait_, " with ", kernel_, " kernel : ",
    signif(time_taken_, 3)
  ))

  # get estimated wiser phenotype, and associated marker data
  wiser_pheno_df <- wiser_obj$wiser_phenotypes
  omic_df <- wiser_obj$wiser_omic_data

  # save wiser phenotype for kernel and trait
  fwrite(wiser_obj$wiser_phenotypes, paste0(
    wiser_pheno_dir_path,
    "wiser_phenotype_estimates_env_", kernel_,
    "_kernel_", trait_, "_", env_, ".csv"
  ), row.names = F, col.names = T)

  # save wiser object for kernel and trait
  saveRDS(wiser_obj, paste0(
    wiser_pheno_dir_path,
    "wiser_obj_", kernel_,
    "_kernel_", trait_, "_", env_
  ))

  # test that x_mat_tilde is indeed orthogonal (i.e. decorrelated) to xi_hat
  x_mat_tilde_mult_xi_test <- t(
    wiser_obj$wiser_x_mat_tilde
  ) %*% wiser_obj$wiser_xi_hat
  # for readability purposes print only the first 10 elements,
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

# merge ls_means and wiser phenotypes by genotype for integrity of analyses
merged_df <- merge(ls_means_df, wiser_pheno_df, by = "Genotype")
ls_means_trait_ <- merged_df[, trait_]
v_hat <- merged_df[, "v_hat"]
omic_df <- omic_df[rownames(omic_df) %in% merged_df$Genotype, ]

# create directory for trait_ graphics if it does not exist
if (!dir.exists(paste0(output_pred_graphics_path, trait_, "/"))) {
  dir.create(paste0(output_pred_graphics_path, trait_, "/"))
}

# make a scatter plot between ls-means and wiser phenotypes with a simple linear
# regression fit
ls_means_wiser_plot_ <- create_scatter_plot_with_linear_fit(
  merged_df, trait_, env_
)

# save ls_means_wiser_plot_
saveWidget(ls_means_wiser_plot_, file = paste0(
  output_pred_graphics_path, trait_, "/ls_means_wiser_pheno_plot_",
  trait_, "_", kernel_, "_kernel_", snp_sample_size_, "_SNP_",
  env_, ".html"
))

# register parallel backend
cl <- makeCluster(nb_cores_)
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
  idx_ <- sample(1:n, size = n, replace = F)
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
      "RF_wiser_pa" = NA,
      "SVR_wiser_pa" = NA,
      "GBLUP_wiser_pa" = NA,
      "GBLUP_wiser_h2" = NA,
      "RKHS_wiser_pa" = NA,
      "LASSO_wiser_pa" = NA,
      "RF_ls_means_pa" = NA,
      "SVR_ls_means_pa" = NA,
      "GBLUP_ls_means_pa" = NA,
      "GBLUP_ls_means_h2" = NA,
      "RKHS_ls_means_pa" = NA,
      "LASSO_ls_means_pa" = NA
    )

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
    fold_result["RF_wiser_pa"] <- cor(
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
      scaled = F, type = "eps-svr",
      kernel = "rbfdot",
      kpar = "automatic", C = c_par, epsilon = 0.1
    )
    f_hat_val_gaussian_svr <- predict(
      gaussian_svr_model,
      as.matrix(omic_df[idx_val, ])
    )
    fold_result["SVR_wiser_pa"] <- cor(
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
    fold_result["GBLUP_wiser_pa"] <- cor(
      f_hat_val_linear_krmm,
      v_hat[idx_val]
    )
    fold_result["GBLUP_wiser_h2"] <- compute_genomic_h2(
      linear_krmm_model$sigma2K_hat,
      linear_krmm_model$sigma2E_hat
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
    fold_result["RKHS_wiser_pa"] <- cor(
      f_hat_val_gaussian_krmm,
      v_hat[idx_val]
    )

    # train and predict with LASSO
    cv_fit_lasso_model <- cv.glmnet(
      intercept = T, y = v_hat[idx_train],
      x = as.matrix(omic_df[idx_train, ]),
      type.measure = "mse", alpha = 1.0, nfold = 10,
      parallel = T
    )
    f_hat_val_lasso <- predict(cv_fit_lasso_model,
      newx = as.matrix(omic_df[idx_val, ]),
      s = "lambda.min"
    )
    fold_result["LASSO_wiser_pa"] <- cor(
      f_hat_val_lasso,
      v_hat[idx_val]
    )

    # training and prediction based on ls-means

    # train and predict with Random Forest
    rf_model <- ranger(
      y = ls_means_trait_[idx_train],
      x = omic_df[idx_train, ],
      mtry = ncol(omic_df) / 3,
      num.trees = 1000
    )
    f_hat_val_rf <- predict(
      rf_model,
      omic_df[idx_val, ]
    )
    fold_result["RF_ls_means_pa"] <- cor(
      f_hat_val_rf$predictions,
      ls_means_trait_[idx_val]
    )

    # train and predict with SVR
    # a correct value for c_par according to Cherkassy and Ma (2004).
    # Neural networks 17, 113-126 is defined as follows
    c_par <- max(
      abs(mean(ls_means_trait_[idx_train])
      + 3 * sd(ls_means_trait_[idx_train])),
      abs(mean(ls_means_trait_[idx_train])
      - 3 * sd(ls_means_trait_[idx_train]))
    )
    gaussian_svr_model <- ksvm(
      x = as.matrix(omic_df[idx_train, ]),
      y = ls_means_trait_[idx_train],
      scaled = F, type = "eps-svr",
      kernel = "rbfdot",
      kpar = "automatic", C = c_par, epsilon = 0.1
    )
    f_hat_val_gaussian_svr <- predict(
      gaussian_svr_model,
      as.matrix(omic_df[idx_val, ])
    )
    fold_result["SVR_ls_means_pa"] <- cor(
      f_hat_val_gaussian_svr,
      ls_means_trait_[idx_val]
    )

    # train and predict with GBLUP (linear kernel krmm)
    linear_krmm_model <- krmm(
      Y = ls_means_trait_[idx_train],
      Matrix_covariates = omic_df[idx_train, ],
      method = "GBLUP"
    )
    f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
      Matrix_covariates = omic_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result["GBLUP_ls_means_pa"] <- cor(
      f_hat_val_linear_krmm,
      ls_means_trait_[idx_val]
    )
    fold_result["GBLUP_ls_means_h2"] <- compute_genomic_h2(
      linear_krmm_model$sigma2K_hat,
      linear_krmm_model$sigma2E_hat
    )

    # train and predict with RKHS (non-linear Gaussian kernel krmm)
    gaussian_krmm_model <- krmm(
      Y = ls_means_trait_[idx_train],
      Matrix_covariates = omic_df[idx_train, ],
      method = "RKHS", kernel = "Gaussian",
      rate_decay_kernel = 0.1
    )
    f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
      Matrix_covariates = omic_df[idx_val, ],
      add_fixed_effects = T
    )
    fold_result["RKHS_ls_means_pa"] <- cor(
      f_hat_val_gaussian_krmm,
      ls_means_trait_[idx_val]
    )

    # train and predict with LASSO
    cv_fit_lasso_model <- cv.glmnet(
      intercept = T, y = ls_means_trait_[idx_train],
      x = as.matrix(omic_df[idx_train, ]),
      type.measure = "mse", alpha = 1.0, nfold = 10,
      parallel = T
    )
    f_hat_val_lasso <- predict(cv_fit_lasso_model,
      newx = as.matrix(omic_df[idx_val, ]),
      s = "lambda.min"
    )
    fold_result["LASSO_ls_means_pa"] <- cor(
      f_hat_val_lasso,
      ls_means_trait_[idx_val]
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

# get methods pa
df_pa_ <- as.data.frame(df_result_[, str_detect(
  colnames(df_result_),
  pattern = "_pa"
)])
colnames(df_pa_) <- str_replace_all(
  colnames(df_pa_),
  pattern = "_pa",
  replacement = ""
)
method_names <- colnames(df_pa_)
df_pa_ <- as.data.frame(apply(df_pa_, 2, as.numeric))

# initialize plot_ly for violin + boxplots
violin_box_plots_pa_ <- plot_ly()

# add violin plots
for (i in seq_along(method_names)) {
  violin_box_plots_pa_ <- add_trace(
    violin_box_plots_pa_,
    type = "violin", # specify violin plot
    y = df_pa_[[method_names[i]]],
    name = method_names[i],
    points = F, # show all points
    jitter = 0.3, # add jitter for spread
    pointpos = -1.8, # adjust point position relative to the violin plot
    marker = list(color = pa_colors_[i]),
    fillcolor = pa_colors_[i],
    line = list(color = pa_colors_[i]),
    meanline = list(visible = F), # do not show mean line
    scalemode = "width", # keep width constant for comparison
    opacity = 0.6 # slight transparency to see the boxplot behind it
  )
}

# add boxplots on top of the violin plots but hide from legend
for (i in seq_along(method_names)) {
  violin_box_plots_pa_ <- add_boxplot(
    violin_box_plots_pa_,
    y = df_pa_[[method_names[i]]],
    name = method_names[i],
    marker = list(color = pa_colors_[i]),
    line = list(color = "black", width = 2), # black line for the boxplot
    fillcolor = "rgba(255,255,255,0)", # transparent fill to see the violin plot underneath
    width = 0.2, # narrower boxplot to fit inside the violin
    notchwidth = 0.4, # add notch for median
    showlegend = F # hide boxplot from legend
  )
}

# add layout
violin_box_plots_pa_ <- violin_box_plots_pa_ %>%
  layout(
    title = paste0(
      "Genomic PA distributions of models for ",
      trait_, ", based on ", ncol(omic_df), " SNP across ",
      n_shuff_, " shuffling scenarios for ", k_folds_, "-folds cv"
    ),
    yaxis = list(
      title = "Predictive ability (PA)",
      range = c(0, 1)
    ),
    legend = list(title = list(text = "Prediction model"))
  )

# save violin_box_plots_pa_ graphics
saveWidget(violin_box_plots_pa_, file = paste0(
  output_pred_graphics_path, trait_, "/genomic_pred_results_",
  trait_, "_", kernel_, "_kernel_", ncol(omic_df), "_SNP_",
  k_folds_, "_folds_CV_", env_, ".html"
))

# add stats and save predictive ability results
df_pa_ <- signif(apply(df_pa_, 2, as.numeric), 2)
rownames(df_pa_) <- paste0("pa_scenario_", 1:nrow(df_pa_))

df_stat <- as.data.frame(rbind(
  apply(df_pa_, 2, mean),
  apply(df_pa_, 2, sd)
))
df_stat <- signif(apply(df_stat, 2, as.numeric), 2)
rownames(df_stat) <- c("pa_mean", "pa_sd")
df_stat <- as.data.frame(df_stat)

df_pa_ <- rbind(df_pa_, df_stat)

fwrite(df_pa_,
  file = paste0(
    output_pred_results_path,
    "genomic_pred_results_", ncol(omic_df), "_SNP_",
    trait_, "_", kernel_, "_kernel_",
    k_folds_, "_folds_CV_", env_, ".csv"
  ), row.names = T
)
print(df_pa_)

# get computed h2 for gblup
df_h2_ <- as.data.frame(df_result_[, str_detect(
  colnames(df_result_),
  pattern = "_h2"
)])
colnames(df_h2_) <- str_replace_all(
  colnames(df_h2_),
  pattern = "_h2",
  replacement = ""
)
method_names <- colnames(df_h2_)
df_h2_ <- as.data.frame(apply(df_h2_, 2, as.numeric))

# initialize plot_ly for violin + boxplots
violin_box_plots_h2_ <- plot_ly()

# add violin plots without points
for (i in seq_along(method_names)) {
  violin_box_plots_h2_ <- add_trace(
    violin_box_plots_h2_,
    type = "violin", # specify violin plot
    y = df_h2_[[method_names[i]]],
    name = method_names[i],
    points = F, # remove points
    jitter = 0.3, # add jitter for spread (not used since points are removed)
    pointpos = -1.8, # adjust point position relative to the violin plot (not used here)
    marker = list(color = h2_colors_[i]),
    fillcolor = h2_colors_[i],
    line = list(color = h2_colors_[i]),
    meanline = list(visible = F), # do not show mean line
    scalemode = "width", # keep width constant for comparison
    opacity = 0.6 # slight transparency to see the boxplot behind it
  )
}

# add boxplots on top of the violin plots but hide from legend
for (i in seq_along(method_names)) {
  violin_box_plots_h2_ <- add_boxplot(
    violin_box_plots_h2_,
    y = df_h2_[[method_names[i]]],
    name = method_names[i],
    marker = list(color = h2_colors_[i]),
    line = list(color = "black", width = 2), # black line for the boxplot
    fillcolor = "rgba(255,255,255,0)", # transparent fill to see the violin plot underneath
    width = 0.2, # narrower boxplot to fit inside the violin
    notchwidth = 0.4, # add notch for median
    showlegend = F # hide boxplot from legend
  )
}

# add layout
violin_box_plots_h2_ <- violin_box_plots_h2_ %>%
  layout(
    title = paste0(
      "Genomic heritability (h2) distributions estimated from GBLUP prediction models for ",
      trait_, ", based on ", ncol(omic_df), " SNP across ",
      n_shuff_, " shuffling scenarios for ", k_folds_, "-folds cv"
    ),
    yaxis = list(
      title = "Genomic heritability (h2)",
      range = c(0, 1)
    ),
    legend = list(title = list(text = "Prediction model"))
  )

# save violin_box_plots_h2_ graphics
saveWidget(violin_box_plots_h2_, file = paste0(
  output_pred_graphics_path, trait_, "/genomic_h2_results_",
  trait_, "_", kernel_, "_kernel_", ncol(omic_df), "_SNP_",
  k_folds_, "_folds_CV_", env_, ".html"
))

# add stats and save h2 results
df_h2_ <- signif(apply(df_h2_, 2, as.numeric), 2)
rownames(df_h2_) <- paste0("pa_scenario_", 1:nrow(df_h2_))

df_stat <- as.data.frame(rbind(
  apply(df_h2_, 2, mean),
  apply(df_h2_, 2, sd)
))
df_stat <- signif(apply(df_stat, 2, as.numeric), 2)
rownames(df_stat) <- c("pa_mean", "pa_sd")
df_stat <- as.data.frame(df_stat)

df_h2_ <- rbind(df_h2_, df_stat)

fwrite(df_h2_,
  file = paste0(
    output_pred_results_path,
    "genomic_h2_results_", ncol(omic_df), "_SNP_",
    trait_, "_", kernel_, "_kernel_",
    k_folds_, "_folds_CV_", env_, ".csv"
  ), row.names = T
)
print(df_h2_)
