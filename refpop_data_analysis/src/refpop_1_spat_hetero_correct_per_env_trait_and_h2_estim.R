# script meant to correct for spatial heterogenity for a specific trait_
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(tidyverse)
library(tidyr)
library(data.table)
library(lubridate)
library(plotly)
library(htmlwidgets)
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
pheno_dir_path_ <- "../../data/phenotype_data/"
pheno_file_path_ <- paste0(pheno_dir_path_, "phenotype_raw_data_no_outliers.csv")
spats_out_file_path <- paste0(pheno_dir_path_, "spats_per_env_adjusted_phenotypes/")
output_pheno_graphics_path <- "../../data/graphics/pheno_graphics/"

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

# define parameters for computations
convert_date_to_days_ <- FALSE # true only if Flowering_begin has not already been converted to days
plot_h2_per_trait_ <- TRUE
h2_mad_value_factor <- 2.5
min_obs_lmer_ <- 5 # cannot fit lmer if less than that.. Note  5 is pretty small
# and doesn't necessarily make sense either, its somewhat arbitrary

# get pheno_df and detect attributes, e.g. number of modalities, for specific variables
pheno_df_ <- as.data.frame(fread(pheno_file_path_))
management_types <- unique(pheno_df_$Management)
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

  # initialize lists for individual-location clonal mean h2 for raw and spats
  # adjusted phenotypes across all managements, and by management type

  # list of h2 by env for raw phenotypes
  env_h2_pheno_vect_ <- rep(0, length(env_list_))
  names(env_h2_pheno_vect_) <- env_list_

  # list of h2 by env for spats adjusted phenotypes
  env_h2_adj_pheno_vect_ <- rep(0, length(env_list_))
  names(env_h2_adj_pheno_vect_) <- env_list_

  # list of h2 by management per env for raw phenotypes
  env_h2_pheno_manage_list_ <- lapply(
    1:n_management,
    function(x) env_h2_pheno_vect_
  )
  names(env_h2_pheno_manage_list_) <- management_types

  # list of h2 by management per env for raw phenotypes
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
        try(
          {
            df_trait_manage_env_ <- df_trait_env_[
              df_trait_env_$Management == manage_type_,
            ]
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
              env_h2_pheno_manage_list_[[manage_type_]][[env_]] <-
                compute_indiv_location_clonal_mean_h2(
                  trait_manage_lmer_mod_, nr_bar_manage_
                )
            }
          },
          silent = TRUE
        )
      }
    }

    # perform a spatial heterogeneity correction for each environment
    list_spats_envir_[[env_]] <- spat_hetero_env_correct_trait(
      trait_,
      env_,
      df_
    )

    # test if any data is available after spatial heterogeneity correction
    if (nrow(list_spats_envir_[[env_]]) > 1) {
      # compute individual-location clonal mean h2 for adjusted phenotypes by env
      df_trait_spats_env_ <- list_spats_envir_[[env_]][, c(
        "Genotype", "Management",
        "Envir", "spats_adj_pheno"
      )]
      # drop na for trait_ to fit regression
      df_trait_spats_env_ <- df_trait_spats_env_ %>% drop_na(all_of("spats_adj_pheno"))

      if (length(unique(df_trait_spats_env_$Genotype)) > min_obs_lmer_) {
        try(
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
          silent = TRUE
        )
        # compute individual clonal mean h2 for each management per environment
        for (manage_type_ in management_types) {
          try(
            {
              df_trait_spats_manage_env_ <- df_trait_spats_env_[
                df_trait_spats_env_$Management == manage_type_,
              ]
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
            silent = TRUE
          )
        }
      }
    }
  }

  # concatenate list elements for spatial heterogeneity correction into a single df_
  df_ <- do.call(rbind, list_spats_envir_)
  df_ <- drop_na(df_)

  # h2 outlier treatment for spats adjusted phenotypes
  env_h2_adj_pheno_vect_ <- env_h2_adj_pheno_vect_[env_h2_adj_pheno_vect_ > 0]
  h2_adj_pheno_mad_out <- get_outliers_vect_mad(env_h2_adj_pheno_vect_, h2_mad_value_factor)
  env_h2_adj_pheno_vect_ <- h2_adj_pheno_mad_out$data_no_outliers

  # exclude identified environments by keeping only those with an h2 over
  # the defined threshold
  env_to_keep_ <- unique(names(env_h2_adj_pheno_vect_))
  df_ <- df_[which(df_$Envir %in% env_to_keep_), ]

  # h2 outlier treatment for raw phenotypes
  env_h2_pheno_vect_ <- env_h2_pheno_vect_[env_h2_pheno_vect_ > 0]
  h2_pheno_mad_out <- get_outliers_vect_mad(env_h2_pheno_vect_, h2_mad_value_factor)
  env_h2_pheno_vect_ <- h2_pheno_mad_out$data_no_outliers

  # h2 outlier treatment for raw phenotypes with specific management
  env_h2_pheno_manage_list_ <- lapply(env_h2_pheno_manage_list_, function(x) x[x > 0])
  h2_pheno_mad_out <- get_outliers_list_mad(env_h2_pheno_manage_list_, h2_mad_value_factor)
  env_h2_pheno_manage_list_ <- h2_pheno_mad_out$data_no_outliers

  # h2 outlier treatment for raw phenotypes with specific management
  env_h2_adj_pheno_manage_list_ <- lapply(env_h2_adj_pheno_manage_list_, function(x) x[x > 0])
  h2_pheno_mad_out <- get_outliers_list_mad(env_h2_adj_pheno_manage_list_, h2_mad_value_factor)
  env_h2_adj_pheno_manage_list_ <- h2_pheno_mad_out$data_no_outliers

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
    spats_out_file_path, trait_,
    "_spats_adjusted_phenotypes_long_format.csv"
  ))

  # get and write retained environments for traits according to h2 outliers removed
  df_trait_env_retain_ <- data.frame(
    "trait" = trait_,
    "env_retained_based_on_h2" =
      paste0(env_to_keep_, collapse = ", ")
  )
  fwrite(df_trait_env_retain_, paste0(
    spats_out_file_path, trait_,
    "_env_retained_based_on_h2.csv"
  ))

  # boxplots of heritabilities for adjusted and raw phenotypes accross environments
  if (plot_h2_per_trait_) {
    # convert to df
    list_raw_adj_data <- vector("list", 6)

    names(list_raw_adj_data) <- c(
      "h2_raw",
      "h2_manage_type_1",
      "h2_manage_type_2",
      "h2_adj",
      "h2_adj_manage_type_1",
      "h2_adj_manage_type_2"
    )

    list_raw_adj_data[["h2_raw"]] <- na.omit(env_h2_pheno_vect_)
    list_raw_adj_data[["h2_manage_type_1"]] <- na.omit(env_h2_pheno_manage_list_[[3]])
    list_raw_adj_data[["h2_manage_type_2"]] <- na.omit(env_h2_pheno_manage_list_[[2]])

    list_raw_adj_data[["h2_adj"]] <- na.omit(env_h2_adj_pheno_vect_)
    list_raw_adj_data[["h2_adj_manage_type_1"]] <- na.omit(env_h2_adj_pheno_manage_list_[[3]])
    list_raw_adj_data[["h2_adj_manage_type_2"]] <- na.omit(env_h2_adj_pheno_manage_list_[[2]])

    # get all registered environments for the computed h2
    all_env_names <- unique(unlist(lapply(list_raw_adj_data, names)))

    # initialize a data frame for the computed h2
    df_h2_raw_adj_and_manage_types_ <- data.frame(
      Environment = all_env_names,
      h2_raw = rep(NA, length(all_env_names)),
      h2_manage_type_1 = rep(NA, length(all_env_names)),
      h2_manage_type_2 = rep(NA, length(all_env_names)),
      h2_adj = rep(NA, length(all_env_names)),
      h2_adj_manage_type_1 = rep(NA, length(all_env_names)),
      h2_adj_manage_type_2 = rep(NA, length(all_env_names))
    )

    # assign the values in the corresponding columns
    for (env_name in all_env_names) {
      if (env_name %in% names(list_raw_adj_data[["h2_raw"]])) {
        df_h2_raw_adj_and_manage_types_[
          df_h2_raw_adj_and_manage_types_$Environment == env_name,
          "h2_raw"
        ] <- list_raw_adj_data[["h2_raw"]][[env_name]]
      }

      if (env_name %in% names(list_raw_adj_data[["h2_manage_type_1"]])) {
        df_h2_raw_adj_and_manage_types_[
          df_h2_raw_adj_and_manage_types_$Environment == env_name,
          "h2_manage_type_1"
        ] <- list_raw_adj_data[["h2_manage_type_1"]][[env_name]]
      }

      if (env_name %in% names(list_raw_adj_data[["h2_manage_type_2"]])) {
        df_h2_raw_adj_and_manage_types_[
          df_h2_raw_adj_and_manage_types_$Environment == env_name,
          "h2_manage_type_2"
        ] <- list_raw_adj_data[["h2_manage_type_2"]][[env_name]]
      }

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
    }

    # create boxplots
    boxplots_ <- plot_ly(df_h2_raw_adj_and_manage_types_, type = "box")

    # add each column as a trace on the graphic
    for (col in names(df_h2_raw_adj_and_manage_types_)[-1]) {
      # add legend
      legend_label <- paste0(col, " (", sum(!is.na(df_h2_raw_adj_and_manage_types_[[col]])), " envir.)")
      # add trace with legend name
      boxplots_ <- add_trace(boxplots_,
        y = df_h2_raw_adj_and_manage_types_[[col]],
        name = legend_label,
        boxpoints = "all", jitter = 0.3, pointpos = -1.8,
        text = paste("Environment: ", df_h2_raw_adj_and_manage_types_$Environment)
      )
    }

    # add axe and title
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
        l = 100,     # adjust margin to create space for the title
        r = 50,      # adjust margin to the right
        b = 100,     # adjust inferior margin
        t = 100      # adjust superior margin
      )
    )

    # add hovermode
    boxplots_ <- boxplots_ %>%
      layout(hovermode = "closest")

    # add annotations next to the graphics
    annotations <- list(
      list(
        text = "&#8226; <b>h2_raw :</b> h2 computed from raw phenotypes for environments",
        x = 1.2, y = 0.2,
        xref = "paper", yref = "paper",
        xanchor = "right", yanchor = "bottom",
        showarrow = FALSE,
        font = list(size = 11)
      ),
      list(
        text = "&#8226; <b>h2_adj :</b> h2 computed from adjusted phenotypes for environments",
        x = 1.2, y = 0.15,
        xref = "paper", yref = "paper",
        xanchor = "right", yanchor = "bottom",
        showarrow = FALSE,
        font = list(size = 11)
      ),
      list(
        text = "&#8226; <b>h2_raw_manage_x :</b> h2 computed for environments associated to management x",
        x = 1.2, y = 0.1,
        xref = "paper", yref = "paper",
        xanchor = "right", yanchor = "bottom",
        showarrow = FALSE,
        font = list(size = 11)
      ),
      list(
        text = "&#8226; <b>h2_adj_manage_x :</b> h2_adj computed for environments associated to management x",
        x = 1.2, y = 0.05,
        xref = "paper", yref = "paper",
        xanchor = "right", yanchor = "bottom",
        showarrow = FALSE,
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
}
# stop cluster
stopCluster(cl)

# reformat data of retained environments for traits according to h2 outlier eviction
file_list <- list.files(
  path = spats_out_file_path,
  pattern = "_env_retained_based_on_h2.csv",
  full.names = TRUE
)

# concatenate all data for these files into a single source and write the output
df_all_trait_env_retain_ <- map_df(file_list, fread)
fwrite(df_all_trait_env_retain_, paste0(
  spats_out_file_path,
  "../all_traits_env_retained_based_on_non_missing_data_and_mad_h2_computed_from_spats_adj_phenotypes.csv"
))

# remove single files for df_trait_env_retain_
file.remove(file_list)
