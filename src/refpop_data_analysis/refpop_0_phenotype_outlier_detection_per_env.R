# script meant to detect phenotype outliers in the raw dataset
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
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
library(rgl)
library(umap)
library(ggplot2)
library(plotly)
# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../functions.R")
# set options to increase memory and suppress warnings
options(expressions = 5e5)
options(warn = -1)

# set paths
# input and output data paths
pheno_dir_path <- "../../data/phenotype_data/"
raw_pheno_file_path <- paste0(pheno_dir_path, "phenotype_raw_data.csv")
pheno_no_outliers_file_path <- paste0(pheno_dir_path, "phenotype_raw_data_no_outliers.csv")
# outlier results data paths
pheno_outliers_results_path <- "../../results/phenotype_outlier_detection/"
pheno_outliers_dir_path <- paste0(pheno_outliers_results_path, "outliers_per_env_phenotypes/")
pheno_outliers_file_path <- paste0(pheno_outliers_results_path, "phenotype_raw_data_outliers.csv")
# output result path for phenotype graphics
output_pheno_graphics_path <- "../../results/graphics/phenotype_graphics/"


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
use_pca_dim_reduc <- TRUE

# knowledge based outlier detection rules
test_sample_size_sup_to_fruit_number_ <- TRUE
test_sample_size_sup_to_value_ <- TRUE
size_value_ <- 20

# define level of risk alpha_ for outlier detection
alpha_ <- 0.01

# plot confidence ellipse for outliers
plot_conf_ellipse_outliers_ <- TRUE

# umap parameters, most sensitive ones
use_umap_ <- F
random_state_umap_ <- 15 # 15, 30 and 50
n_neighbors_umap_ <- 15 # 15, 30, 50
min_dist_ <- 0.1

# read raw pheno data and define proxy for outlier detection
df_raw_ <- as.data.frame(fread(raw_pheno_file_path))
df_proxy_ <- df_raw_

# perform imputation or read already imputed data for df_proxy_
if (use_miss_forest_imput_for_outlier_detection_ &&
  !file.exists(paste0(
    pheno_dir_path,
    "phenotype_miss_forest_imputed_data_proxy.csv"
  ))) {
  # set cluster for parallelized computations
  cl <- makeCluster(n_traits_)
  registerDoParallel(cl)
  doRNG::registerDoRNG(seed = 123)

  df_proxy_[, selected_traits_] <- missForest(df_proxy_[, selected_traits_],
    parallelize = "forests",
    maxiter = 10,
    ntree = 100,
    verbose = T
  )$ximp

  # stop cluster
  stopCluster(cl)

  fwrite(df_proxy_,
    file = paste0(
      pheno_dir_path,
      "phenotype_miss_forest_imputed_data_proxy.csv"
    ),
    sep = ","
  )
} else if (use_miss_forest_imput_for_outlier_detection_) {
  df_proxy_ <- as.data.frame(fread(paste0(
    pheno_dir_path,
    "phenotype_miss_forest_imputed_data_proxy.csv"
  ), sep = ","))
}

# get unique environments (i.e.location-year) for trait_
env_years_ <- unique(str_extract(unique(df_proxy_$Envir), "\\d{4}"))
env_list_ <- reordered_cols(unique(df_proxy_$Envir),
  prefix_order_patterns = env_years_
)

# define list of rownames which will be used to delete outliers from df_raw_
list_rownames_outliers_env_ <- vector("list", length(env_list_))
names(list_rownames_outliers_env_) <- env_list_

for (env_ in env_list_) {
  print(env_)
  # get df for raw and imputed data by env_
  df_raw_env_ <- df_raw_[df_raw_$Envir == env_, ]
  df_proxy_env_ <- df_proxy_[df_proxy_$Envir == env_, ]

  # 1st outlier detection test based on data knowledge
  idx_rule_one_outliers_ <- NULL
  if (test_sample_size_sup_to_fruit_number_) {
    idx_rule_one_outliers_ <- which(df_raw_env_$Sample_size > df_raw_env_$Fruit_number)
  }

  # 2nd outlier detection test based on data knowledge
  idx_rule_two_outliers_ <- NULL
  if (test_sample_size_sup_to_value_) {
    idx_rule_two_outliers_ <- which(df_raw_env_$Sample_size > size_value_)
  }

  # 3rd outlier detection test based on PCA and Mahalanobis distance
  idx_three_outliers_ <- NULL
  tryCatch(
    {
      if (use_pca_dim_reduc) {
        # compute pca to reduce the number of dimensions for Mahalanobis distance
        # in an attempt to avoid the curse of dimensionality
        pca_obj <- PCA(df_proxy_env_[, selected_traits_],
          ncp = length(selected_traits_),
          graph = F
        )
        n_comp_required <- n_comp_required_for_percent_explained_var(pca_obj,
          percent_explained_var = 99.999999
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
      idx_three_outliers_ <- as.numeric(which(maha_dist_ > maha_dist_threshold_))
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
      idx_three_outliers_ <<- as.numeric(which(maha_dist_ > maha_dist_threshold_))
    }
  )

  # get all outliers indices and df_raw_env_ outliers
  idx_outliers_ <- sort(union(
    union(idx_rule_one_outliers_, idx_rule_two_outliers_),
    idx_three_outliers_
  ))
  df_raw_env_outliers_ <- df_raw_env_[idx_outliers_, ]
  list_rownames_outliers_env_[[env_]] <- rownames(df_raw_env_outliers_)

  # get location and scale parameters used for outlier detection
  df_loc_scale_ <- data.frame(
    "trait" = names(location_),
    "mean_used_for_outliers_detection" = location_,
    "standard_deviation_used_for_outliers_detection" = scale_
  )

  # write results
  fwrite(df_raw_env_outliers_, file = paste0(
    pheno_outliers_dir_path, env_,
    "_phenotype_outliers.csv"
  ))
  fwrite(df_loc_scale_,
    file = paste0(
      pheno_outliers_dir_path, env_,
      "_location_scale_parameters.csv"
    )
  )
}

# get all indices for outliers
idx_outliers_df_raw_ <- as.numeric(unlist(list_rownames_outliers_env_))

# write df_raw_ slice with outliers
fwrite(df_raw_[idx_outliers_df_raw_, ],
  file = pheno_outliers_file_path
)

# write df_raw_ slice without outliers
fwrite(df_raw_[-idx_outliers_df_raw_, ],
  file = pheno_no_outliers_file_path
)

# plots

# create a confidence ellipse for first env outlier detection
if (plot_conf_ellipse_outliers_) {
  # perform pca again but with graph = TRUE

  # get data for env_
  for (env_ in env_list_) {
    print(env_)

    tryCatch({
      df_raw_env_ <- df_raw_[df_raw_$Envir == env_, ]
      df_proxy_env_ <- df_proxy_[df_proxy_$Envir == env_, ]

      # 1st outlier detection test based on data knowledge
      idx_rule_one_outliers_ <- NULL
      if (test_sample_size_sup_to_fruit_number_) {
        idx_rule_one_outliers_ <- which(df_raw_env_$Sample_size > df_raw_env_$Fruit_number)
      }

      # 2nd outlier detection test based on data knowledge
      idx_rule_two_outliers_ <- NULL
      if (test_sample_size_sup_to_value_) {
        idx_rule_two_outliers_ <- which(df_raw_env_$Sample_size > size_value_)
      }

      # 3rd outlier detection test based on PCA and Mahalanobis distance
      idx_three_outliers_ <- NULL
      tryCatch(
        {
          # compute pca to reduce the number of dimensions for Mahalanobis distance
          # in an attempt to avoid the curse of dimensionality
          pca_obj <- PCA(df_proxy_env_[, selected_traits_],
            ncp = length(selected_traits_),
            graph = TRUE
          )
          n_comp_required <- n_comp_required_for_percent_explained_var(pca_obj,
            percent_explained_var = 99.999999
          )
          pca_comp_mat <- pca_obj$ind$coord[, 1:n_comp_required]

          # compute Mahalanobis distance based on minimum covariance determinant (MCD)
          mcd_obj <- covMcd(pca_comp_mat)
          maha_dist_ <- mahalanobis(pca_comp_mat,
            center = mcd_obj$center,
            cov = mcd_obj$cov
          )

          # for pca case retrieve location and scale of original variables
          mcd_obj <- covMcd(df_proxy_env_[, selected_traits_])
          location_ <- signif(mcd_obj$center, 2)
          scale_ <- signif(sqrt(diag(mcd_obj$cov)), 2)

          maha_dist_threshold_ <- quantile(maha_dist_, probs = 1 - alpha_)
          idx_three_outliers_ <- as.numeric(which(maha_dist_ > maha_dist_threshold_))
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
          idx_three_outliers_ <<- as.numeric(which(maha_dist_ > maha_dist_threshold_))
        }
      )

      # get all outliers indices and df_raw_env_ outliers
      idx_outliers_ <- sort(union(
        union(idx_rule_one_outliers_, idx_rule_two_outliers_),
        idx_three_outliers_
      ))

      # plot Mahalanobis dist against chi square with p = n_comp_required degrees of freedom to
      # test for normality hypothesis of data, the latter is visibly not multivariate normal
      # set.seed(123)
      # maha_dist_vect <- sort(maha_dist_)
      # chisq_vect <- sort(rchisq(
      #   n = length(maha_dist_),
      #   df = n_comp_required, ncp = 0
      # ))
      # qqplot(maha_dist_vect, chisq_vect)

      # plot cumulative percentage of explained variance as a function of the number of components
      n_comp_required <- n_comp_required_for_percent_explained_var(pca_obj,
        percent_explained_var = 99.999999
      )
      df <- data.frame(
        comp = 1:nrow(pca_obj$eig),
        cumulative_percentage = pca_obj$eig[, 3]
      )
      ggplot(df, aes(x = comp, y = cumulative_percentage)) +
        geom_line() +
        labs(
          x = "Number of components",
          y = "Cumulative percentage of explained variance",
          title = paste0("Cumulative percentage of explained variance by \n number of components for ", env_)
        ) +
        theme_minimal() +
        geom_vline(xintercept = n_comp_required, linetype = "dashed", color = "red") +
        annotate("text",
          x = n_comp_required - 1.5, y = max(df$cumulative_percentage) - 5,
          label = "Nb. of components \n required to reach \n at least 99.99% of \n explained variance",
          color = "red",
          size = 2
        ) +
        theme(
          # title
          plot.title = element_text(hjust = 0.5, size = 8),
          # axis titles
          axis.title = element_text(size = 8),
          # axis text
          axis.text = element_text(size = 6),
          # legend title
          legend.title = element_text(size = 8),
          # legend text
          legend.text = element_text(size = 6),
        )
      ggsave(paste0(
        output_pheno_graphics_path,
        "_outlier_detection/nb_comp_99_variance_",
        env_, ".pdf"
      ))

      # 2D plot of outliers using pca, note that two first components may not explain all the variance,
      # hence the outliers may not all lie outside the 99%-confidence ellipse
      n_comp <- 2
      pdf(paste0(
        output_pheno_graphics_path,
        "_outlier_detection/pca_cor_circle_",
        env_, ".pdf"
      ))
      plot(pca_obj,
        choix = "varcor",
        title = paste0("PCA graph of variables for ", env_)
      )
      dev.off()

      pca_comp_mat <- as.data.frame(pca_obj$ind$coord[, 1:n_comp_required])
      colnames_pca_ <- colnames(pca_comp_mat)

      for (i in 1:length(colnames_pca_)) {
        colnames_pca_[i] <- paste0(
          colnames_pca_[i], " (",
          signif(pca_obj$eig[i, 2][1], 4), "%)"
        )
      }
      colnames(pca_comp_mat) <- colnames_pca_

      data_outliers <- pca_comp_mat[idx_outliers_, 1:n_comp]
      data_non_outliers <- pca_comp_mat[-idx_outliers_, 1:n_comp]

      ggplot() +
        geom_point(
          data = data_non_outliers,
          aes(
            x = !!rlang::sym(colnames_pca_[1]),
            y = !!rlang::sym(colnames_pca_[2]),
            color = "Non-Outlier"
          )
        ) +
        geom_point(
          data = data_outliers,
          aes(
            x = !!rlang::sym(colnames_pca_[1]),
            y = !!rlang::sym(colnames_pca_[2]),
            color = "Outlier"
          )
        ) +
        stat_ellipse(
          data = pca_comp_mat,
          aes(
            x = !!rlang::sym(colnames_pca_[1]),
            y = !!rlang::sym(colnames_pca_[2])
          ),
          level = 0.99, geom = "polygon",
          fill = "red", alpha = 0.2
        ) +
        scale_color_manual(
          values = c("Non-Outlier" = "blue", "Outlier" = "red"),
          labels = c("Non-Outlier", "Outlier")
        ) +
        labs(color = "Outlier status") +
        xlab(colnames_pca_[1]) +
        ylab(colnames_pca_[2]) +
        theme_minimal() +
        ggtitle(paste0(
          "Two dimensional 99%-confidence ellipse associated with ", env_,
          ", for outlier \n detection in phenotypes using PCA and Mahalanobis distance across ",
          ncol(pca_comp_mat), " dimensions."
        )) +
        theme(
          # title
          plot.title = element_text(hjust = 0.5, size = 8),
          # axis titles
          axis.title = element_text(size = 8),
          # axis text
          axis.text = element_text(size = 6),
          # legend title
          legend.title = element_text(size = 8),
          # legend text
          legend.text = element_text(size = 6)
        )
      ggsave(
        paste0(
          output_pheno_graphics_path,
          "_outlier_detection/2d_plot_outlier_detection_",
          env_, ".pdf"
        )
      )

      # 3D plot of outliers using pca, note that two first components may not explain all the variance,
      # hence the outliers may not all lie outside the 99%-confidence ellipse
      pca_df_ <- pca_comp_mat
      pca_df_ <- pca_df_[, 1:3]

      # add a column to indicate outliers
      pca_df_$Outlier_status <- NA
      pca_df_$Outlier_status[idx_outliers_] <- "Outlier"
      pca_df_$Outlier_status[-idx_outliers_] <- "Non outlier"

      rgl.open()
      view3d(
        theta = 60, phi = 30, fov = 60, zoom = 1.5,
        type = "userviewpoint"
      )
      # create a new 3D plot
      plot3d(
        pca_df_[, -match(
          "Outlier_status",
          colnames(pca_df_)
        )],
        col = ifelse(pca_df_$Outlier_status == "Outlier", "red", "blue"),
        type = "s",
        size = 4,
        main = paste0(
          "Three-dimensional 99%-confidence ellipse associated with ",
          env_, " (", ncol(pca_comp_mat), " dimensions used for OD)"
        ),
        cex = 0.6
      )
      # draw the confidence ellipse
      ellipse <- ellipse3d(
        cov(pca_df_[, -match(
          "Outlier_status",
          colnames(pca_df_)
        )]),
        centre = colMeans(pca_df_[, -match(
          "Outlier_status",
          colnames(pca_df_)
        )]),
        level = 0.99
      )
      shade3d(ellipse, col = "pink", alpha = 0.5)
      # add legend
      legend3d(
        "topright",
        legend = c("Outlier", "Non outlier"),
        col = c("red", "blue"),
        pch = 16,
        magnify = 2,
        cex = 1.2
      )
      rgl.snapshot(paste0(
        output_pheno_graphics_path,
        "_outlier_detection/3d_plot_outlier_detection_",
        env_, ".png"
      ))
      rgl.clear()
    })
  }
}
