# script meant to detect phenotype outliers in the raw dataset
# note: text is formatted from Addins using Style active file from styler package
library(data.table)
library(stringr)
library(FactoMineR)
library(rstudioapi)
library(doParallel)
library(doRNG)
library(robustbase)
library(foreach)
library(parallel)
library(missForest)
library(Matrix)
setwd(dirname(getActiveDocumentContext()$path))
source("../../functions.R")

# set paths
pheno_dir_path <- "../../data/phenotype_data/"
pheno_out_dir_path <- "../../data/phenotype_data/outliers_per_env_phenotypes/"
raw_pheno_file_path <- paste0(pheno_dir_path, "phenotype_raw_data.csv")
clean_pheno_file_path <- paste0(pheno_dir_path, "phenotype_data.csv")

# define selected_traits_ and vars_to_keep_ for output
selected_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Scab_fruits", "Weight_sample",
  "Sample_size"
)
n_traits_ <- length(selected_traits_)

# define status for imputation using miss forest
use_miss_forest_imput_for_outlier_detection_ <- TRUE

# use pca for dimension reduction
use_pca_dim_reduc <- FALSE

# knowledge based outlier detection rules
test_sample_size_sup_to_fruit_number_ <- TRUE
test_sample_size_sup_to_value_ <- TRUE
size_value_ <- 20

# define level of risk alpha_ for outlier detection
alpha_ <- 0.01

# read raw pheno data
df_raw_ <- as.data.frame(fread(raw_pheno_file_path))
df_proxy_ <- df_raw_

# perform imputation or read already imputed data for df_proxy_
if (use_miss_forest_imput_for_outlier_detection_ &&
  !file.exists(paste0(
    pheno_dir_path,
    "phenotype_miss_forest_imputed_data.csv"
  ))) {
  # set cluster for parallelized computations
  cl <- makeCluster(n_traits_)
  registerDoParallel(cl)
  doRNG::registerDoRNG(seed = 123)

  df_proxy_[, selected_traits_] <- missForest(df_proxy_[, selected_traits_],
    parallelize = "forests",
    maxiter = 2,
    ntree = 50,
    verbose = T
  )$ximp

  # stop cluster
  stopCluster(cl)

  fwrite(df_proxy_,
    file = paste0(
      pheno_dir_path,
      "phenotype_miss_forest_imputed_data.csv"
    ),
    sep = ","
  )
} else if (use_miss_forest_imput_for_outlier_detection_) {
  df_proxy_ <- as.data.frame(fread(paste0(
    pheno_dir_path,
    "phenotype_miss_forest_imputed_data.csv"
  ), sep = ","))
}

# get unique environments (i.e.location-year) for trait_
env_years_ <- unique(str_extract(unique(df_proxy_$Envir), "\\d{4}"))
env_list_ <- reordered_cols(unique(df_proxy_$Envir),
  prefix_order_patterns = env_years_
)

# define list of rownames which will be used to delete outliers from df_raw_
list_rownames_out_env_ <- vector("list", length(env_list_))
names(list_rownames_out_env_) <- env_list_

for (env_ in env_list_) {
  print(env_)
  # get df for raw and imputed data by env_
  df_raw_env_ <- df_raw_[df_raw_$Envir == env_, ]
  df_proxy_env_ <- df_proxy_[df_proxy_$Envir == env_, ]

  # 1st outlier detection test based on data knowledge
  idx_rule_one_out_ <- NULL
  if (test_sample_size_sup_to_fruit_number_) {
    idx_rule_one_out_ <- which(df_raw_env_$Sample_size > df_raw_env_$Fruit_number)
  }

  # 2nd outlier detection test based on data knowledge
  idx_rule_two_out_ <- NULL
  if (test_sample_size_sup_to_value_) {
    idx_rule_two_out_ <- which(df_raw_env_$Sample_size > size_value_)
  }

  # 3rd outlier detection test based on PCA and Mahalanobis distance
  idx_three_out_ <- NULL
  tryCatch(
    {
      if (use_pca_dim_reduc) {
        # compute pca to reduce the number of dimensions for Mahalanobis distance
        # in an attempt to avoid the curse of dimensionality
        pca_obj <- PCA(df_proxy_env_[, selected_traits_],
          ncp = length(selected_traits_),
          graph = T
        )
        n_comp_required <- n_comp_required_for_percent_explained_var(pca_obj,
          percent_explained_var = 100
        )
        pca_comp_mat <- pca_obj$ind$coord[, 1:n_comp_required]

        # compute Mahalanobis distance based on minimum covariance determinant (MCD)
        mcd_obj <- covMcd(pca_comp_mat)
        maha_dist_ <- mahalanobis(pca_comp_mat,
          center = mcd_obj$center,
          cov = mcd_obj$cov
        )
      } else {
        # compute Mahalanobis distance based on minimum covariance determinant (MCD)
        mcd_obj <- covMcd(df_proxy_env_[, selected_traits_])
        maha_dist_ <- mahalanobis(df_proxy_env_[, selected_traits_],
          center = mcd_obj$center,
          cov = mcd_obj$cov
        )
      }
      # for pca case retrieve location and scale of original variables
      mcd_obj <- covMcd(df_proxy_env_[, selected_traits_])
      location_ <- signif(mcd_obj$center, 2)
      scale_ <- signif(sqrt(diag(mcd_obj$cov)), 2)

      maha_dist_threshold_ <- quantile(maha_dist_, probs = 1 - alpha_)
      idx_three_out_ <- as.numeric(which(maha_dist_ > maha_dist_threshold_))
    },

    # if covariance matrix is singular:
    error = function(e) {
      # find the nearest positive definite covariance matrix if it is singular
      cov_mat <- nearPD(cov(df_proxy_env_[, selected_traits_]))$mat
      location_ <- signif(colMeans(df_proxy_env_[, selected_traits_]), 2)
      scale_ <- signif(sqrt(diag(cov_mat)), 2)
      maha_dist_ <- mahalanobis(df_proxy_env_[, selected_traits_],
        center = location_,
        cov = cov_mat
      )

      maha_dist_threshold_ <- quantile(maha_dist_, probs = 1 - alpha_)
      idx_three_out_ <<- as.numeric(which(maha_dist_ > maha_dist_threshold_))
    }
  )

  # get all outliers indices and df_raw_env_ outliers
  idx_out_ <- sort(union(
    union(idx_rule_one_out_, idx_rule_two_out_),
    idx_three_out_
  ))
  df_raw_env_out_ <- df_raw_env_[idx_out_, ]
  list_rownames_out_env_[[env_]] <- rownames(df_raw_env_out_)

  # get location and scale parameters used for outlier detection
  df_loc_scale_ <- data.frame(
    "mcd_mean_used_for_out_detection" = location_,
    "mcd_standard_deviation_used_for_out_detection" = scale_
  )

  # write results
  fwrite(df_raw_env_out_, file = paste0(
    pheno_out_dir_path, env_,
    "_phenotype_outliers.csv"
  ))
  fwrite(df_loc_scale_, file = paste0(
    pheno_out_dir_path, env_,
    "_location_scale_parameters.csv"
  ))
}

# get all indices for outliers
idx_out_df_raw_ <- as.numeric(unlist(list_rownames_out_env_))

# write df_raw_ without outliers
fwrite(df_raw_[-idx_out_df_raw_, ],
  file = clean_pheno_file_path
)
