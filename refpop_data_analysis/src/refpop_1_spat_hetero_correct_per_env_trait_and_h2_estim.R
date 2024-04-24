# script meant to correct for spatial heterogenity for a specific trait_
# note: text is formatted from Addins using Style active file from styler package
library(tidyverse)
library(tidyr)
library(data.table)
library(lubridate)
library(ggplot2)
library(emmeans)
library(SpATS)
library(stringr)
library(rstudioapi)
library(lme4)
library(anytime)
library(foreach)
library(parallel)
library(doParallel)
options(warn = -1)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
setwd(dirname(getActiveDocumentContext()$path))
source("../../functions.R")

# set paths
pheno_file_path_ <- "../../data/phenotype_data/phenotype_raw_data_no_outliers.csv"
output_file_path <- "../../data/phenotype_data/spats_per_env_adjusted_phenotypes/"
output_pheno_graphics_path <- "../../data/graphics/pheno_graphics/"

# define function(s) and package(s) to export for parallelization
func_to_export_ <- c("fread")
pkgs_to_export_ <- c(
  "data.table", "stringr", "SpATS", "lme4",
  "lubridate", "emmeans", "ggplot2", "tidyr"
)

# define selected_traits_ and vars_to_keep_ for output
selected_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Scab_fruits", "Weight_sample",
  "Sample_size"
)

vars_to_keep_ <- c(
  "Envir", "Management", "Row",
  "Position", "Genotype"
)

# define parameters for computations
convert_date_to_days_ <- FALSE # true only if Flowering_begin has not already been converted to days
plot_h2_per_trait_ <- TRUE
h2_mad_value_factor <- 2.5
min_obs_lmer_ <- 5 # cannot fit lmer if less than that.. Note  5 is pretty small
# and doesn't necessarily make sense either, its somewhat arbitrary

# get pheno_df and detect attributes, e.g. number of modalities, for specific variables
pheno_df_ <- as.data.frame(fread(pheno_file_path_))
management_types <- unique(df_$Management)
n_management <- length(management_types)

# parallelize treatments for each trait_, for sequential treatment replace %dopar% by %do%
foreach(
  trait_ = selected_traits_,
  .export = func_to_export_,
  .packages = pkgs_to_export_
) %dopar% {
  print(paste0("performing computation for ", trait_))

  # keep variables of interest
  df_ <- pheno_df_[, c(vars_to_keep_, trait_)]

  # if trait_ is flowering start, convert date to days if necessary
  if (identical(trait_, "Flowering_begin") && convert_date_to_days_) {
    df_[, trait_] <- as.numeric(strftime(
      as.data.frame(df_)[, trait_],
      format = "%j"
    ))
  }

  # get unique environments (i.e.location-year) for trait_
  env_years_ <- unique(str_extract(unique(df_$Envir), "\\d{4}"))
  env_list_ <- reordered_cols(unique(df_$Envir),
    prefix_order_patterns = env_years_
  )

  # compute individual-location clonal mean h2 for raw phenotypes across
  # all managements and by management type
  env_h2_pheno_vect_ <- rep(0, length(env_list_))
  names(env_h2_pheno_vect_) <- env_list_

  env_h2_manage_list_ <- lapply(1:n_management, function(x) env_h2_pheno_vect_)
  names(env_h2_manage_list_) <- management_types

  for (env_ in env_list_)
  {
    df_trait_env_ <- df_[df_$Envir == env_, c(
      "Genotype", "Management",
      "Envir", trait_
    )]
    # drop na for trait_ for regression
    df_trait_env_ <- df_trait_env_ %>% drop_na(all_of(trait_))

    if (length(unique(df_trait_env_$Genotype)) > min_obs_lmer_) {
      # get average number of replications for genotypes
      nr_bar_ <- mean(table(df_trait_env_$Genotype))
      try(
        {
          # compute lmer mixed model  y = mu + g + eps with g and eps as random factors
          trait_lmer_mod_ <- lmer(as.formula(paste0(trait_, " ~ 1 + (1 | Genotype)")),
            data = df_trait_env_
          )
          # compute individual clonal mean h2 for each environment
          env_h2_pheno_vect_[[env_]] <- compute_indiv_location_clonal_mean_h2(
            trait_lmer_mod_, nr_bar_
          )
        },
        silent = TRUE
      )
      # compute individual clonal mean h2 for each management per environment
      for (manage_type_ in management_types) {
        print(env_)
        try(
          {
            df_trait_manage_env_ <- df_trait_env_[
              df_trait_env_$Management == manage_type_,
            ]
            if (nrow(df_trait_manage_env_) > 1 &&
              length(unique(df_trait_manage_env_$Genotype)) > min_obs_lmer_) {
              nr_bar_manage_env_ <- mean(table(df_trait_manage_env_$Genotype))

              trait_lmer_mod_ <- lmer(
                as.formula(paste0(
                  trait_,
                  " ~ 1 + (1 | Genotype)"
                )),
                data = df_trait_manage_env_
              )
              env_h2_manage_list_[[manage_type_]][[env_]] <-
                compute_indiv_location_clonal_mean_h2(
                  trait_lmer_mod_, nr_bar_manage_env_
                )
            }
          },
          silent = TRUE
        )
      }
    }
  }

  # perform a spatial heterogeneity correction for each environment
  list_spats_envir_ <- vector("list", length(env_list_))
  names(list_spats_envir_) <- env_list_
  for (env_ in env_list_) {
    print(env_)
    list_spats_envir_[[env_]] <- spat_hetero_env_correct_trait(
      trait_,
      env_,
      df_
    )
  }
  df_ <- do.call(rbind, list_spats_envir_)
  df_ <- drop_na(df_)

  # compute individual-location clonal mean h2 for adjusted phenotypes
  env_h2_adj_pheno_vect_ <- rep(0, length(env_list_))
  names(env_h2_adj_pheno_vect_) <- env_list_
  for (env_ in env_list_)
  {
    df_trait_env_ <- df_[df_$Envir == env_, c("Genotype", "Envir", "spats_adj_pheno")]

    if (length(unique(df_trait_env_$Genotype)) > min_obs_lmer_) {
      try(
        {
          # get average number of replications for genotypes
          nr_bar_ <- mean(table(df_trait_env_$Genotype))

          # compute lmer mixed model  y = mu + g + eps with g and eps as random factors
          trait_lmer_mod_ <- lmer(spats_adj_pheno ~ 1 + (1 | Genotype),
            data = df_trait_env_
          )
          # compute individual clonal mean h2 for each environment
          env_h2_adj_pheno_vect_[[env_]] <- compute_indiv_location_clonal_mean_h2(
            trait_lmer_mod_, nr_bar_
          )
        },
        silent = TRUE
      )
    }
  }

  # keep envir with h2 different from zero (i.e. envir with no data)
  env_h2_adj_pheno_vect_ <- env_h2_adj_pheno_vect_[env_h2_adj_pheno_vect_ > 0]

  # apply mad to detect outlier(s)
  h2_mad_value <- mad(env_h2_adj_pheno_vect_, constant = 1)
  h2_median <- median(env_h2_adj_pheno_vect_)
  h2_threshold <- h2_mad_value_factor * h2_mad_value

  # detect environment(s) to be excluded from adjusted pheno computed h2
  idx_outliers_h2_adj_pheno_ <- which(abs(env_h2_adj_pheno_vect_ -
    h2_median) > h2_threshold &
    env_h2_adj_pheno_vect_ < h2_median)

  outliers_h2_adj_pheno_ <- names(env_h2_adj_pheno_vect_)[
    idx_outliers_h2_adj_pheno_
  ]

  # remove the environments to be excluded from the lists of heritabilities
  if (length(idx_outliers_h2_adj_pheno_) > 0) {
    env_h2_adj_pheno_vect_ <- env_h2_adj_pheno_vect_[-match(
      outliers_h2_adj_pheno_,
      names(env_h2_adj_pheno_vect_)
    )]
  }

  # exclude identified environments by keeping only those with an h2 over
  # the defined threshold
  env_to_keep_ <- unique(names(env_h2_adj_pheno_vect_))
  df_ <- df_[which(df_$Envir %in% env_to_keep_), ]

  # apply the same treatment with mad, as above for h2 computed from adjusted phenotypes,
  # for h2 computed from raw phenotypes
  # keep envir with h2 different from zero (i.e. envir with no data)
  env_h2_pheno_vect_ <- env_h2_pheno_vect_[env_h2_pheno_vect_ > 0]

  # apply mad to detect outlier(s)
  h2_mad_value <- mad(env_h2_pheno_vect_, constant = 1)
  h2_median <- median(env_h2_pheno_vect_)
  h2_threshold <- h2_mad_value_factor * h2_mad_value

  # detect environment(s) to be excluded from adjusted pheno computed h2
  idx_outliers_h2_pheno_ <- which(abs(env_h2_pheno_vect_ -
    h2_median) > h2_threshold &
    env_h2_pheno_vect_ < h2_median)

  outliers_h2_pheno_ <- names(env_h2_pheno_vect_)[
    idx_outliers_h2_pheno_
  ]

  # remove the environments to be excluded from the lists of heritabilities
  if (length(idx_outliers_h2_pheno_) > 0) {
    env_h2_pheno_vect_ <- env_h2_pheno_vect_[-match(
      outliers_h2_pheno_,
      names(env_h2_pheno_vect_)
    )]
  }

  # define rename exceptions
  exception_cols <- c(
    "Genotype", "Envir",
    "Management", "Row", "Position",
    "R", "P", trait_
  )

  # rename columns excluding the exception columns
  new_names <- colnames(df_)
  new_names[!(new_names %in% exception_cols)] <- paste0(trait_, "_", new_names[
    !(new_names %in% exception_cols)
  ])

  # replace the existing column names with the new names
  colnames(df_) <- new_names

  # write adjusted phenotype, from spatial heterogeneity correction, to long format
  fwrite(df_, paste0(
    output_file_path, trait_,
    "_spats_adjusted_phenotypes_long_format.csv"
  ))

  # get and write retained environments for traits according to h2 outlier eviction
  df_trait_retain_env_ <- data.frame(
    "trait" = trait_,
    "retained_env_based_on_h2" =
      paste0(env_to_keep_, collapse = ", ")
  )
  fwrite(df_trait_retain_env_, paste0(
    output_file_path, trait_,
    "_retained_env_based_on_h2.csv"
  ))

  # boxplots of heritabilities for adjusted and raw phenotypes accross environments
  if (plot_h2_per_trait_) {
    # convert to df
    df_h2_adj <- data.frame(
      Environment = names(env_h2_adj_pheno_vect_),
      h2_adj = env_h2_adj_pheno_vect_
    )
    df_h2 <- data.frame(
      Environment = names(env_h2_pheno_vect_),
      h2 = env_h2_pheno_vect_
    )

    # merge df
    df_combined <- merge(df_h2_adj, df_h2, by = "Environment")

    # create boxplots
    ggplot(data = df_combined) +
      geom_boxplot(
        aes(
          x = "h2 computed from adjusted phenotypes corrected with SpATS",
          y = h2_adj
        ),
        width = 0.25, fill = "blue",
        alpha = 0.7, position = position_dodge(width = 0.75)
      ) +
      geom_boxplot(aes(x = "h2 computed from raw phenotypes", y = h2),
        width = 0.25, fill = "red",
        alpha = 0.7, position = position_dodge(width = 0.75)
      ) +
      scale_x_discrete(
        labels = c(
          "h2 computed from adjusted phenotypes \n corrected with SpATS",
          "h2 computed from raw phenotypes"
        ),
        name = NULL
      ) + # Supprimer le nom de l'axe x
      labs(
        y = paste0("Individual-location clonal mean heritability (h2)"),
        title = paste0("Comparison of individual-location clonal mean heritability (h2)
                        \n computed from raw and adjusted phenotypes for ", trait_)
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 1)

    ggsave(paste0(
      output_pheno_graphics_path, "individual_location_clonal_mean_h2_",
      trait_, ".jpg"
    ), dpi = 600)
  }
}
# stop cluster
stopCluster(cl)

# reformat data of retained environments for traits according to h2 outlier eviction
file_list <- list.files(
  path = output_file_path,
  pattern = "_retained_env_based_on_h2.csv",
  full.names = TRUE
)

# concatenate all data for these files into a single source and write the output
df_all_trait_retain_env_ <- map_df(file_list, fread)
fwrite(df_all_trait_retain_env_, paste0(
  output_file_path,
  "../all_traits_retained_env_based_on_non_missing_data_h2_computed_from_spats_adj_phenotypes.csv"
))

# remove single files for df_trait_retain_env_
file.remove(file_list)
