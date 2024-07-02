# script meant to perform genomic prediction and analyses for refpop
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(MASS)
library(data.table)
library(stringr)
library(FactoMineR)
library(doParallel)
library(doRNG)
library(robustbase)
library(foreach)
library(parallel)
library(missForest)
library(Matrix)
library(rgl)
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
computation_mode <- "cluster"

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
trait_ <- "Trunk_diameter"

# define shift seed value by
mult_seed_by_ <- 100

# set k for K-folds cv
k_folds_ <- 5

# define number of shuffles
n_shuff_ <- 20

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

# get phenotype and genotype data, and family and origin info
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
pheno_df <- merged_df[, selected_traits_]
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

# get number of phenotypes
n <- nrow(pheno_df)

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
      "RF" = NA, "SVR" = NA, "RKHS" = NA, "GBLUP" = NA, "LASSO" = NA,
      "SVR_support_vectors" = NA
    )

    # train and predict with Random Forest
    rf_model <- ranger(
      y = pheno_df[idx_train, trait_],
      x = geno_df[idx_train, ],
      mtry = ncol(geno_df) / 3,
      num.trees = 1000
    )
    f_hat_val_rf <- predict(
      rf_model,
      geno_df[idx_val, ]
    )
    fold_result["RF"] <- cor(
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
      x = as.matrix(geno_df[idx_train, ]),
      y = pheno_df[idx_train, trait_],
      scaled = FALSE, type = "eps-svr",
      kernel = "rbfdot",
      kpar = "automatic", C = c_par, epsilon = 0.1
    )
    idx_sv_ <- SVindex(gaussian_svr_model)
    sv_ <- pheno_df[idx_train, "Genotype"][idx_sv_]
    f_hat_val_gaussian_svr <- predict(
      gaussian_svr_model,
      as.matrix(geno_df[idx_val, ])
    )
    fold_result["SVR"] <- cor(
      f_hat_val_gaussian_svr,
      pheno_df[idx_val, trait_]
    )
    fold_result["SVR_support_vectors"] <- paste0(sv_, collapse = ", ")

    # train and predict with GBLUP (linear kernel krmm)
    linear_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = geno_df[idx_train, ],
      method = "GBLUP"
    )
    f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
      Matrix_covariates = geno_df[idx_val, ],
      add_flxed_effects = T
    )
    fold_result["GBLUP"] <- cor(
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
      add_flxed_effects = T
    )
    fold_result["RKHS"] <- cor(
      f_hat_val_gaussian_krmm,
      pheno_df[idx_val, trait_]
    )

    # train and predict with LASSO
    cv_fit_lasso_model <- cv.glmnet(
      intercept = TRUE, y = pheno_df[idx_train, trait_],
      x = as.matrix(geno_df[idx_train, ]),
      type.measure = "mse", alpha = 1.0, nfold = 10,
      parallel = TRUE
    )
    f_hat_val_lasso <- predict(cv_fit_lasso_model,
      newx = as.matrix(geno_df[idx_val, ]),
      s = "lambda.min"
    )
    fold_result["LASSO"] <- cor(
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

# convert to data frame format
df_result_ <- as.data.frame(df_result_)
df_ <- df_result_[, -match(
  "SVR_support_vectors",
  colnames(df_result_)
)]

# get methods names
df_ <- as.data.frame(apply(df_, 2, as.numeric))
method_names <- colnames(df_)

# initialize plot_ly boxplot graphic
boxplots_rpa_ <- plot_ly()

# add boxplots
for (method_ in method_names) {
  boxplots_rpa_ <- add_boxplot(
    boxplots_rpa_,
    y = df_[[method_]],
    name = method_,
    boxpoints = "all",
    jitter = 0.3,
    pointpos = -1.8
  )
}
# add layout
boxplots_rpa_ <- boxplots_rpa_ %>%
  layout(
    title = paste0(
      "Genomic prediction PA distributions of methods for ",
      trait_, ", based on ", snp_sample_size_, " SNPs across ",
      n_shuff_, " shuffling scenarios"
    ),
    yaxis = list(title = "Relative Prediction Accuracy (PA)", range = c(0, 1)),
    legend = list(title = list(text = "Prediction method"))
  )

# save boxplots_rpa_ graphics
saveWidget(boxplots_rpa_, file = paste0(
  output_pred_graphics_path, trait_, "/relative_prediction_accuracy_",
  trait_, "_", snp_sample_size_, "_SNP", ".html"
))

# save relative prediction accuracy results
rownames(df_) <- paste0("rpa_shuff_", 1:nrow(df_))
df_stat <- as.data.frame(rbind(apply(df_, 2, mean), apply(df_, 2, sd)))
rownames(df_stat) <- c("rpa_mean", "rpa_sd")
df_ <- rbind(df_, df_stat)
fwrite(df_,
  file = paste0(
    output_pred_results_path,
    "genomic_pred_results_", snp_sample_size_, "_SNP_",
    trait_, ".csv"
  ), row.names = T
)
