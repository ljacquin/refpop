# script meant to correct for spatial heterogenity within environments for all traits
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
library(devtools)
if ("refpop_env" %in% conda_list()$name) {
  use_condaenv("refpop_env")
}
library(tidyverse)
library(tidyr)
library(data.table)
library(lubridate)
library(plotly)
library(htmlwidgets)
library(emmeans)
library(SpATS)
library(stringr)
library(lme4)
library(anytime)
library(foreach)
library(parallel)
library(doParallel)

# define computation mode, i.e. "local" or "cluster"
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

# set paths
pheno_dir_path_ <- "../../data/phenotype_data/"

pheno_file_path_ <- paste0(
  pheno_dir_path_,
  "raw_phenotype_data_correc_manage_type_no_md_outliers.csv"
)

output_spats_file_path <- paste0(
  pheno_dir_path_,
  "spats_per_env_adjusted_phenotypes/"
)

# outlier results data paths
pheno_outliers_results_path <- "../../results/phenotype_outlier_detection/"

# output result path for phenotype graphics
output_pheno_graphics_path <- "../../results/graphics/phenotype_graphics/"

# define function(s) and package(s) to export for parallelization
func_to_export_ <- c("fread")
pkgs_to_export_ <- c(
  "data.table", "stringr", "SpATS", "lme4",
  "lubridate", "emmeans", "plotly", "tidyr", "htmlwidgets"
)

# define selected_traits_ and vars_to_keep_ for output, note Sample_size cannot
# be considered as a real trait
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

excluded_pseudo_trait_for_save_ <- "Sample_size"

# colors for boxplots
blue_gradient <- c("#3D9BC5", "#005AB5", "#00407A", "#002A66")
yellow_orange_gradient <- colorRampPalette(c("#149414", "#3b5534"))(4)
bp_colors_ <- c(blue_gradient, yellow_orange_gradient)

# define parameters for computations
convert_date_to_days_ <- F # true only if Flowering_begin has not already been converted to days
plot_h2_per_trait_ <- T
h2_mad_value_factor <- 2.5
min_obs_lmer_ <- 5 # cannot fit lmer if less than that.. Note  5 is pretty small
# and doesn't necessarily make sense either, its somewhat arbitrary

# get pheno_df and detect attributes, e.g. number of modalities or levels for specific variables
pheno_df_ <- as.data.frame(fread(pheno_file_path_))
management_types <- unique(pheno_df_$Management)
n_management <- length(management_types)

# define a list for missing data, model singularity and h2 outliers
# to save these infos
miss_data_singular_model_h2_out_list_ <<- vector("list", length(selected_traits_))
names(miss_data_singular_model_h2_out_list_) <- selected_traits_

# parallelize treatments for each trait_, for sequential treatment replace %dopar% by %do%
# save and return errors in miss_data_singular_model_h2_out_vect_
cl <- makeCluster(detectCores())
registerDoParallel(cl)

miss_data_singular_model_h2_out_vect_ <-
  foreach(
    trait_ = selected_traits_,
    .export = func_to_export_,
    .packages = pkgs_to_export_,
    .combine = c
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

    # get unique environments (i.e.location-year-management) for trait_
    env_years_ <- unique(str_extract(unique(df_$Envir), "\\d{4}"))
    env_list_ <- reordered_cols(unique(df_$Envir),
      prefix_order_patterns = env_years_
    )

    # initialize lists for individual-location clonal mean h2 for raw and spats
    # adjusted phenotypes across all managements, and by management type

    # list of h2 by env for raw phenotypes
    env_h2_raw_pheno_vect_ <- rep(0, length(env_list_))
    names(env_h2_raw_pheno_vect_) <- env_list_

    # list of h2 by env for spats adjusted phenotypes
    env_h2_adj_pheno_vect_ <- rep(0, length(env_list_))
    names(env_h2_adj_pheno_vect_) <- env_list_

    # list of h2 by management per env for raw phenotypes
    env_h2_raw_pheno_manage_list_ <- lapply(
      1:n_management,
      function(x) env_h2_raw_pheno_vect_
    )
    names(env_h2_raw_pheno_manage_list_) <- management_types

    # list of h2 by management per env for adj phenotypes
    env_h2_adj_pheno_manage_list_ <- lapply(
      1:n_management,
      function(x) env_h2_adj_pheno_vect_
    )
    names(env_h2_adj_pheno_manage_list_) <- management_types

    # initialize list for spatial heterogeneity correction for each environment
    list_spats_envir_ <- vector("list", length(env_list_))
    names(list_spats_envir_) <- env_list_

    # perform h2 computations by environment management type for raw and
    # adjusted phenotypes for analyzed trait
    for (env_ in env_list_)
    {
      df_trait_env_ <- df_[df_$Envir == env_, c(
        "Genotype", "Management",
        "Envir", trait_
      )]

      # drop na for trait_ to fit regression
      df_trait_env_ <- df_trait_env_ %>% drop_na(all_of(trait_))
      df_trait_env_ <- droplevels(df_trait_env_)

      if (!is.null(df_trait_env_) && nrow(df_trait_env_) > 1) {
        if (length(unique(df_trait_env_$Genotype)) > min_obs_lmer_) {
          # get average number of replications for genotypes
          nr_bar_ <- mean(table(df_trait_env_$Genotype))
          tryCatch(
            {
              # compute lmer mixed model  y = mu + g + eps with g and eps as random factors
              trait_lmer_mod_ <- lmer(as.formula(paste0(trait_, " ~ 1 + (1 | Genotype)")),
                data = df_trait_env_
              )
              # compute individual clonal mean h2 for each environment
              env_h2_raw_pheno_vect_[[env_]] <- compute_indiv_location_clonal_mean_h2(
                trait_lmer_mod_, nr_bar_
              )
            },
            error = function(e) {
              condition_message_ <<- as.character(capture.output(
                conditionMessage(e)
              ))
              miss_data_singular_model_h2_out_list_[[trait_]] <<- c(
                miss_data_singular_model_h2_out_list_[[trait_]],
                paste0(
                  "Error for ", env_,
                  " during individual clonal mean h2 raw computation: ",
                  condition_message_
                )
              )
              condition_message_ <<- NULL
            }
          )
          # compute individual clonal mean h2 for each management per environment
          for (manage_type_ in management_types) {
            tryCatch(
              {
                df_trait_manage_env_ <- df_trait_env_[
                  df_trait_env_$Management == manage_type_,
                ]
                df_trait_manage_env_ <- droplevels(df_trait_manage_env_)

                if (nrow(df_trait_manage_env_) > 1 &&
                  length(unique(df_trait_manage_env_$Genotype)) > min_obs_lmer_) {
                  nr_bar_manage_ <- mean(table(df_trait_manage_env_$Genotype))

                  trait_manage_lmer_mod_ <- lmer(
                    as.formula(paste0(
                      trait_,
                      " ~ 1 + (1 | Genotype)"
                    )),
                    data = df_trait_manage_env_
                  )
                  env_h2_raw_pheno_manage_list_[[manage_type_]][[env_]] <-
                    compute_indiv_location_clonal_mean_h2(
                      trait_manage_lmer_mod_, nr_bar_manage_
                    )
                }
              },
              error = function(e) {
                condition_message_ <<- as.character(capture.output(
                  conditionMessage(e)
                ))
                miss_data_singular_model_h2_out_list_[[trait_]] <<- c(
                  miss_data_singular_model_h2_out_list_[[trait_]],
                  paste0(
                    "Error for ", env_, ", for management type ", manage_type_,
                    ", during individual clonal mean h2 raw computation: ",
                    condition_message_
                  )
                )
                condition_message_ <<- NULL
              }
            )
          }
        } else {
          miss_data_singular_model_h2_out_list_[[trait_]] <- c(
            miss_data_singular_model_h2_out_list_[[trait_]],
            paste0(
              "Error for ", env_, ": too few genotypes in environment"
            )
          )
        }
      } else {
        miss_data_singular_model_h2_out_list_[[trait_]] <- c(
          miss_data_singular_model_h2_out_list_[[trait_]],
          paste0(
            "Error for ", env_, ": no data for trait in environment "
          )
        )
      }

      # perform a spatial heterogeneity correction for each environment
      list_spats_env_ <- spat_hetero_env_correct_trait(
        trait_,
        env_,
        df_
      )
      list_spats_envir_[[env_]] <- list_spats_env_$df_envir_
      if (list_spats_env_$message_ != "no error") {
        miss_data_singular_model_h2_out_list_[[trait_]] <- c(
          miss_data_singular_model_h2_out_list_[[trait_]],
          paste0(
            "Error for ", env_, " during SpATs: ",
            list_spats_env_$message_
          )
        )
      }

      # test if any data is available after spatial heterogeneity correction
      if (!is.null(list_spats_envir_[[env_]]) &&
        nrow(list_spats_envir_[[env_]]) > 1) {
        # compute individual-location clonal mean h2 for adjusted phenotypes by env
        df_trait_spats_env_ <- list_spats_envir_[[env_]][, c(
          "Genotype", "Management",
          "Envir", "spats_adj_pheno"
        )]
        # drop na for trait_ to fit regression
        df_trait_spats_env_ <- df_trait_spats_env_ %>% drop_na(all_of("spats_adj_pheno"))
        df_trait_spats_env_ <- droplevels(df_trait_spats_env_)

        if (length(unique(df_trait_spats_env_$Genotype)) > min_obs_lmer_) {
          tryCatch(
            {
              # get average number of replications for genotypes
              nr_bar_spats_ <- mean(table(df_trait_spats_env_$Genotype))

              # compute lmer mixed model  y = mu + g + eps with g and eps as random factors
              trait_spats_lmer_mod_ <- lmer(spats_adj_pheno ~ 1 + (1 | Genotype),
                data = df_trait_spats_env_
              )
              # compute individual clonal mean h2 for each environment
              env_h2_adj_pheno_vect_[[env_]] <- compute_indiv_location_clonal_mean_h2(
                trait_spats_lmer_mod_, nr_bar_spats_
              )
            },
            error = function(e) {
              condition_message_ <<- as.character(capture.output(
                conditionMessage(e)
              ))
              miss_data_singular_model_h2_out_list_[[trait_]] <<- c(
                miss_data_singular_model_h2_out_list_[[trait_]],
                paste0(
                  "Error for ", env_,
                  " after SpATs during individual clonal mean h2 adj computation: ",
                  condition_message_
                )
              )
              condition_message_ <<- NULL
            }
          )
          # compute individual clonal mean h2 for each management per environment
          for (manage_type_ in management_types) {
            tryCatch(
              {
                df_trait_spats_manage_env_ <- df_trait_spats_env_[
                  df_trait_spats_env_$Management == manage_type_,
                ]
                df_trait_spats_manage_env_ <- droplevels(df_trait_spats_manage_env_)

                if (nrow(df_trait_spats_manage_env_) > 1 &&
                  length(unique(df_trait_spats_manage_env_$Genotype)) > min_obs_lmer_) {
                  nr_bar_spats_manage_ <- mean(table(df_trait_spats_manage_env_$Genotype))

                  trait_spats_manage_lmer_mod_ <- lmer(
                    spats_adj_pheno ~ 1 + (1 | Genotype),
                    data = df_trait_spats_manage_env_
                  )
                  env_h2_adj_pheno_manage_list_[[manage_type_]][[env_]] <-
                    compute_indiv_location_clonal_mean_h2(
                      trait_spats_manage_lmer_mod_, nr_bar_spats_manage_
                    )
                }
              },
              error = function(e) {
                condition_message_ <<- as.character(capture.output(
                  conditionMessage(e)
                ))
                miss_data_singular_model_h2_out_list_[[trait_]] <<- c(
                  miss_data_singular_model_h2_out_list_[[trait_]],
                  paste0(
                    "Error for ", env_, " after SpATs, for management type ", manage_type_,
                    ", during individual clonal mean h2 adj computation: ",
                    condition_message_
                  )
                )
                condition_message_ <<- NULL
              }
            )
          }
        } else {
          miss_data_singular_model_h2_out_list_[[trait_]] <- c(
            miss_data_singular_model_h2_out_list_[[trait_]],
            paste0(
              "Error for ", env_, " after SpATS: too few genotypes in environment"
            )
          )
        }
      }
    }

    # concatenate list elements for spatial heterogeneity correction into a single df_
    df_ <- do.call(rbind, list_spats_envir_)
    df_ <- drop_na(df_)

    # get h2 outliers associated to raw phenotypes
    env_h2_raw_pheno_vect_ <- env_h2_raw_pheno_vect_[env_h2_raw_pheno_vect_ > 0]
    h2_raw_pheno_mad_out <- get_outliers_vect_mad(
      env_h2_raw_pheno_vect_,
      h2_mad_value_factor
    )
    if (length(h2_raw_pheno_mad_out$idx_outliers) > 0) {
      env_h2_raw_out_ <- env_h2_raw_pheno_vect_[h2_raw_pheno_mad_out$idx_outliers]
      miss_data_singular_model_h2_out_list_[[trait_]] <- c(
        miss_data_singular_model_h2_out_list_[[trait_]],
        paste0(
          "h2_raw outlier detected for ",
          names(env_h2_raw_out_), ": ",
          as.numeric(env_h2_raw_out_)
        )
      )
    }
    env_h2_raw_pheno_vect_ <- h2_raw_pheno_mad_out$data_no_outliers

    # get h2 outliers associated to spats adjusted phenotypes
    env_h2_adj_pheno_vect_ <- env_h2_adj_pheno_vect_[env_h2_adj_pheno_vect_ > 0]
    h2_adj_pheno_mad_out <- get_outliers_vect_mad(
      env_h2_adj_pheno_vect_,
      h2_mad_value_factor
    )
    if (length(h2_adj_pheno_mad_out$idx_outliers) > 0) {
      env_h2_adj_out_ <- env_h2_adj_pheno_vect_[h2_adj_pheno_mad_out$idx_outliers]
      miss_data_singular_model_h2_out_list_[[trait_]] <- c(
        miss_data_singular_model_h2_out_list_[[trait_]],
        paste0(
          "h2_adj outlier detected for ",
          names(env_h2_adj_out_), ": ",
          as.numeric(env_h2_adj_out_)
        )
      )
    }
    env_h2_adj_pheno_vect_ <- h2_adj_pheno_mad_out$data_no_outliers

    # get h2 outliers associated to raw phenotypes for type 1 & 2 managements
    env_h2_raw_pheno_manage_list_ <- lapply(
      env_h2_raw_pheno_manage_list_,
      function(x) x[x > 0]
    )
    h2_raw_pheno_manage_types_mad_out <- get_outliers_list_mad(
      env_h2_raw_pheno_manage_list_,
      h2_mad_value_factor
    )
    # get outliers for h2_raw type 1 management
    if (length(h2_raw_pheno_manage_types_mad_out$idx_outliers[[2]]) > 0) {
      env_h2_raw_manage_type_1_out_ <- env_h2_raw_pheno_manage_list_[[2]][
        h2_raw_pheno_manage_types_mad_out$idx_outliers[[2]]
      ]
      miss_data_singular_model_h2_out_list_[[trait_]] <- c(
        miss_data_singular_model_h2_out_list_[[trait_]],
        paste0(
          "h2_raw management type 1 outlier detected for ",
          names(env_h2_raw_manage_type_1_out_), ": ",
          as.numeric(env_h2_raw_manage_type_1_out_)
        )
      )
    }
    # get outliers for h2_raw type 2 management
    if (length(h2_raw_pheno_manage_types_mad_out$idx_outliers[[4]]) > 0) {
      env_h2_raw_manage_type_2_out_ <- env_h2_raw_pheno_manage_list_[[4]][
        h2_raw_pheno_manage_types_mad_out$idx_outliers[[4]]
      ]
      miss_data_singular_model_h2_out_list_[[trait_]] <- c(
        miss_data_singular_model_h2_out_list_[[trait_]],
        paste0(
          "h2_raw management type 2 outlier detected for ",
          names(env_h2_raw_manage_type_2_out_), ": ",
          as.numeric(env_h2_raw_manage_type_2_out_)
        )
      )
    }
    # get outliers for h2_raw type 3 management
    if (length(h2_raw_pheno_manage_types_mad_out$idx_outliers[[3]]) > 0) {
      env_h2_raw_manage_type_3_out_ <- env_h2_raw_pheno_manage_list_[[3]][
        h2_raw_pheno_manage_types_mad_out$idx_outliers[[3]]
      ]
      miss_data_singular_model_h2_out_list_[[trait_]] <- c(
        miss_data_singular_model_h2_out_list_[[trait_]],
        paste0(
          "h2_raw management type 3 outlier detected for ",
          names(env_h2_raw_manage_type_3_out_), ": ",
          as.numeric(env_h2_raw_manage_type_3_out_)
        )
      )
    }

    env_h2_raw_pheno_manage_list_ <- h2_raw_pheno_manage_types_mad_out$data_no_outliers

    # get h2 outliers associated to adjusted phenotypes for type 1 & 2 managements
    env_h2_adj_pheno_manage_list_ <- lapply(
      env_h2_adj_pheno_manage_list_,
      function(x) x[x > 0]
    )
    h2_adj_pheno_manage_types_mad_out <- get_outliers_list_mad(
      env_h2_adj_pheno_manage_list_,
      h2_mad_value_factor
    )
    # get outliers for h2_adj type 1 management
    if (length(h2_adj_pheno_manage_types_mad_out$idx_outliers[[2]]) > 0) {
      env_h2_adj_manage_type_1_out_ <- env_h2_adj_pheno_manage_list_[[2]][
        h2_adj_pheno_manage_types_mad_out$idx_outliers[[2]]
      ]
      miss_data_singular_model_h2_out_list_[[trait_]] <- c(
        miss_data_singular_model_h2_out_list_[[trait_]],
        paste0(
          "h2_adj management type 1 outlier detected for ",
          names(env_h2_adj_manage_type_1_out_), ": ",
          as.numeric(env_h2_adj_manage_type_1_out_)
        )
      )
    }
    # get outliers for h2_adj type 2 management
    if (length(h2_adj_pheno_manage_types_mad_out$idx_outliers[[4]]) > 0) {
      env_h2_adj_manage_type_2_out_ <- env_h2_adj_pheno_manage_list_[[4]][
        h2_adj_pheno_manage_types_mad_out$idx_outliers[[4]]
      ]
      miss_data_singular_model_h2_out_list_[[trait_]] <- c(
        miss_data_singular_model_h2_out_list_[[trait_]],
        paste0(
          "h2_adj management type 2 outlier detected for ",
          names(env_h2_adj_manage_type_2_out_), ": ",
          as.numeric(env_h2_adj_manage_type_2_out_)
        )
      )
    }
    # get outliers for h2_adj type 3 management
    if (length(h2_adj_pheno_manage_types_mad_out$idx_outliers[[3]]) > 0) {
      env_h2_adj_manage_type_3_out_ <- env_h2_adj_pheno_manage_list_[[3]][
        h2_adj_pheno_manage_types_mad_out$idx_outliers[[3]]
      ]
      miss_data_singular_model_h2_out_list_[[trait_]] <- c(
        miss_data_singular_model_h2_out_list_[[trait_]],
        paste0(
          "h2_adj management type 3 outlier detected for ",
          names(env_h2_adj_manage_type_3_out_), ": ",
          as.numeric(env_h2_adj_manage_type_3_out_)
        )
      )
    }
    env_h2_adj_pheno_manage_list_ <- h2_adj_pheno_manage_types_mad_out$data_no_outliers

    # detect BUFFER management, if any, and remove them
    idx_manage_type_buffer <- which(str_detect(miss_data_singular_model_h2_out_list_[[trait_]],
      pattern = "BUFFER"
    ))
    if (length(idx_manage_type_buffer) > 0) {
      miss_data_singular_model_h2_out_list_[[trait_]] <-
        miss_data_singular_model_h2_out_list_[[trait_]][-idx_manage_type_buffer]
    }

    # sort list
    miss_data_singular_model_h2_out_list_[[trait_]] <-
      sort(miss_data_singular_model_h2_out_list_[[trait_]], decreasing = T)

    # keep environments considered as non outliers for adjusted phenotypes since
    # this data will be used for further analysis and prediction
    env_to_keep_ <- unique(names(env_h2_adj_pheno_vect_))
    df_ <- df_[which(df_$Envir %in% env_to_keep_), ]

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
      output_spats_file_path, trait_,
      "_spats_adjusted_phenotypes_long_format.csv"
    ))

    # get and write retained environments for traits according to h2 inliers
    df_trait_env_retain_ <- data.frame(
      "trait" = trait_,
      "env_retained_based_on_inlier_h2" =
        paste0(env_to_keep_, collapse = ", ")
    )
    fwrite(df_trait_env_retain_, paste0(
      output_spats_file_path, trait_,
      "_env_retained_based_on_inlier_h2.csv"
    ))

    # boxplots of heritabilities for adjusted and raw phenotypes accross environments
    if (plot_h2_per_trait_) {
      # convert to df
      list_raw_adj_data <- vector("list", 8)

      names(list_raw_adj_data) <- c(
        "h2_raw",
        "h2_raw_manage_type_1",
        "h2_raw_manage_type_2",
        "h2_adj_manage_type_3",
        "h2_adj",
        "h2_adj_manage_type_1",
        "h2_adj_manage_type_2",
        "h2_adj_manage_type_3"
      )

      list_raw_adj_data[["h2_raw"]] <- na.omit(env_h2_raw_pheno_vect_)
      list_raw_adj_data[["h2_raw_manage_type_1"]] <- na.omit(env_h2_raw_pheno_manage_list_[[2]])
      list_raw_adj_data[["h2_raw_manage_type_2"]] <- na.omit(env_h2_raw_pheno_manage_list_[[4]])
      list_raw_adj_data[["h2_raw_manage_type_3"]] <- na.omit(env_h2_raw_pheno_manage_list_[[3]])

      list_raw_adj_data[["h2_adj"]] <- na.omit(env_h2_adj_pheno_vect_)
      list_raw_adj_data[["h2_adj_manage_type_1"]] <- na.omit(env_h2_adj_pheno_manage_list_[[2]])
      list_raw_adj_data[["h2_adj_manage_type_2"]] <- na.omit(env_h2_adj_pheno_manage_list_[[4]])
      list_raw_adj_data[["h2_adj_manage_type_3"]] <- na.omit(env_h2_raw_pheno_manage_list_[[3]])

      # get all registered environments for the computed h2
      all_env_names <- unique(unlist(lapply(list_raw_adj_data, names)))

      # initialize a data frame for the computed h2
      df_h2_raw_adj_and_manage_types_ <- data.frame(
        Environment = all_env_names,
        h2_raw = rep(NA, length(all_env_names)),
        h2_raw_manage_type_1 = rep(NA, length(all_env_names)),
        h2_raw_manage_type_2 = rep(NA, length(all_env_names)),
        h2_raw_manage_type_3 = rep(NA, length(all_env_names)),
        h2_adj = rep(NA, length(all_env_names)),
        h2_adj_manage_type_1 = rep(NA, length(all_env_names)),
        h2_adj_manage_type_2 = rep(NA, length(all_env_names)),
        h2_adj_manage_type_3 = rep(NA, length(all_env_names))
      )

      # assign the values in the corresponding columns
      for (env_name in all_env_names) {
        ##
        if (env_name %in% names(list_raw_adj_data[["h2_raw"]])) {
          df_h2_raw_adj_and_manage_types_[
            df_h2_raw_adj_and_manage_types_$Environment == env_name,
            "h2_raw"
          ] <- list_raw_adj_data[["h2_raw"]][[env_name]]
        }
        if (env_name %in% names(list_raw_adj_data[["h2_raw_manage_type_1"]])) {
          df_h2_raw_adj_and_manage_types_[
            df_h2_raw_adj_and_manage_types_$Environment == env_name,
            "h2_raw_manage_type_1"
          ] <- list_raw_adj_data[["h2_raw_manage_type_1"]][[env_name]]
        }
        if (env_name %in% names(list_raw_adj_data[["h2_raw_manage_type_2"]])) {
          df_h2_raw_adj_and_manage_types_[
            df_h2_raw_adj_and_manage_types_$Environment == env_name,
            "h2_raw_manage_type_2"
          ] <- list_raw_adj_data[["h2_raw_manage_type_2"]][[env_name]]
        }
        if (env_name %in% names(list_raw_adj_data[["h2_raw_manage_type_3"]])) {
          df_h2_raw_adj_and_manage_types_[
            df_h2_raw_adj_and_manage_types_$Environment == env_name,
            "h2_raw_manage_type_3"
          ] <- list_raw_adj_data[["h2_raw_manage_type_3"]][[env_name]]
        }
        ##
        if (env_name %in% names(list_raw_adj_data[["h2_adj"]])) {
          df_h2_raw_adj_and_manage_types_[
            df_h2_raw_adj_and_manage_types_$Environment == env_name,
            "h2_adj"
          ] <- list_raw_adj_data[["h2_adj"]][[env_name]]
        }
        if (env_name %in% names(list_raw_adj_data[["h2_adj_manage_type_1"]])) {
          df_h2_raw_adj_and_manage_types_[
            df_h2_raw_adj_and_manage_types_$Environment == env_name,
            "h2_adj_manage_type_1"
          ] <- list_raw_adj_data[["h2_adj_manage_type_1"]][[env_name]]
        }
        if (env_name %in% names(list_raw_adj_data[["h2_adj_manage_type_2"]])) {
          df_h2_raw_adj_and_manage_types_[
            df_h2_raw_adj_and_manage_types_$Environment == env_name,
            "h2_adj_manage_type_2"
          ] <- list_raw_adj_data[["h2_adj_manage_type_2"]][[env_name]]
        }
        if (env_name %in% names(list_raw_adj_data[["h2_adj_manage_type_3"]])) {
          df_h2_raw_adj_and_manage_types_[
            df_h2_raw_adj_and_manage_types_$Environment == env_name,
            "h2_adj_manage_type_3"
          ] <- list_raw_adj_data[["h2_adj_manage_type_3"]][[env_name]]
        }
      }

      # create boxplots
      boxplots_ <- plot_ly(df_h2_raw_adj_and_manage_types_, type = "box")

      # initialize a counter for the colors
      color_counter <- 1

      # add each column as a trace on the graphic
      for (col in names(df_h2_raw_adj_and_manage_types_)[-1]) {
        # add legend
        legend_label <- paste0(col, " (", sum(!is.na(df_h2_raw_adj_and_manage_types_[[col]])), " envir.)")

        # assign color based on the counter
        color_to_use <- bp_colors_[color_counter]

        # increment the color counter for the next loop
        color_counter <- color_counter + 1

        # reset color counter if it exceeds the number of available colors
        if (color_counter > length(bp_colors_)) {
          color_counter <- 1
        }

        # add trace with legend name
        boxplots_ <- add_trace(
          boxplots_,
          y = df_h2_raw_adj_and_manage_types_[[col]],
          name = legend_label,
          boxpoints = "all",
          jitter = 0.3,
          pointpos = -1.8,
          marker = list(color = color_to_use),
          fillcolor = color_to_use,
          line = list(color = color_to_use),
          text = paste("Environment: ", df_h2_raw_adj_and_manage_types_$Environment)
        )
      }

      # add axes and title
      boxplots_ <- layout(boxplots_,
        yaxis = list(title = "Individual-location clonal mean heritability (h2)", range = c(0, 1.0)),
        xaxis = list(title = ""),
        title = list(
          text = paste0(
            "Individual-location clonal mean heritability (h2) for ", trait_,
            ", computed from raw and adjusted \n phenotypes obtained with SpATS, for the considered environments and management types"
          ),
          font = list(size = 15)
        ),
        margin = list(
          l = 100, # adjust margin to create space for the title
          r = 50, # adjust margin to the right
          b = 100, # adjust inferior margin
          t = 100 # adjust superior margin
        )
      )

      # add hovermode
      boxplots_ <- boxplots_ %>%
        layout(hovermode = "closest")

      # add annotations
      annotations <- list(
        list(
          text = "&#8226; <b>h2_raw :</b> h2 computed from raw phenotypes for environments",
          x = 1.2, y = 0.2,
          xref = "paper", yref = "paper",
          xanchor = "right", yanchor = "bottom",
          showarrow = F,
          font = list(size = 11)
        ),
        list(
          text = "&#8226; <b>h2_adj :</b> h2 computed from adjusted phenotypes for environments",
          x = 1.2, y = 0.15,
          xref = "paper", yref = "paper",
          xanchor = "right", yanchor = "bottom",
          showarrow = F,
          font = list(size = 11)
        ),
        list(
          text = "&#8226; <b>h2_raw_manage_x :</b> h2_raw computed for environments associated to management x",
          x = 1.2, y = 0.1,
          xref = "paper", yref = "paper",
          xanchor = "right", yanchor = "bottom",
          showarrow = F,
          font = list(size = 11)
        ),
        list(
          text = "&#8226; <b>h2_adj_manage_x :</b> h2_adj computed for environments associated to management x",
          x = 1.2, y = 0.05,
          xref = "paper", yref = "paper",
          xanchor = "right", yanchor = "bottom",
          showarrow = F,
          font = list(size = 11)
        )
      )

      # add annotations to layout
      boxplots_ <- boxplots_ %>%
        layout(annotations = annotations)

      if (!identical(trait_, excluded_pseudo_trait_for_save_)) {
        # save boxplots_ graphics
        saveWidget(boxplots_, file = paste0(
          output_pheno_graphics_path, "individual_location_clonal_mean_h2_",
          trait_, ".html"
        ))
      }
    }
    return(paste0(
      trait_, ": ",
      miss_data_singular_model_h2_out_list_[[trait_]]
    ))
  }

# stop cluster
stopCluster(cl)

# write errors saved in miss_data_singular_model_h2_out_vect_
writeLines(miss_data_singular_model_h2_out_vect_, paste0(
  pheno_outliers_results_path,
  "envir_per_trait_with_miss_data_or_singular_model_or_h2_outliers.csv"
))

# reformat data of retained environments for traits according to h2 outlier eviction
file_list <- list.files(
  path = output_spats_file_path,
  pattern = "_env_retained_based_on_inlier_h2.csv",
  full.names = T
)

# concatenate all data for these files into a single source and write the output
df_all_trait_env_retain_ <- map_df(file_list, fread)
fwrite(df_all_trait_env_retain_, paste0(
  pheno_outliers_results_path,
  "envir_per_trait_retained_based_on_inlier_h2_for_adj_phenotypes.csv"
))

# write new phenotype data source with retained environments
sel_env_ <- unique(str_split(df_all_trait_env_retain_[[2]], pattern = ", ")[[2]])
pheno_df_ <- pheno_df_[pheno_df_$Envir %in% sel_env_, ]

fwrite(pheno_df_, paste0(
  pheno_dir_path_,
  "phenotype_data.csv"
))

# remove single files for df_trait_env_retain_
file.remove(file_list)
