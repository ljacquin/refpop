# script meant to correct for spatial heterogenity for a specific trait_
# note: text is formatted from Addins using Style active file from styler package
library(tidyverse)
library(tidyr)
library(lubridate)
library(ggplot2)
library(emmeans)
library(SpATS)
library(stringr)
library(rstudioapi)
library(data.table)
library(lme4)
library(anytime)
setwd(dirname(getActiveDocumentContext()$path))
source("../../functions.R")

# set paths
pheno_file_path_ <- "../../data/phenotype_data/phenotype_data.csv"
output_file_path <- "../../data/phenotype_data/spats_adjusted_phenotypes/"
output_pheno_graphics_path <- "../../data/graphics/pheno_graphics/"

# define selected_traits_, vars_to_keep_ for output and parameters for computations
selected_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Scab_fruits", "Weight_sample"
)
vars_to_keep_ <- c("Envir", "Management", "Row", "Position", "Genotype")

convert_date_to_days_ <- FALSE # true only if Flowering_begin has not already been converted to days 
use_ggplot2_ <- FALSE
h2_mad_value_factor <- 2.5
min_obs_lmer_ <- 5            # cannot fit lmer if less than that.. Note  5 is pretty small
                              # and doesn't necessarily make sense either, its somewhat arbitrary

for (trait_ in selected_traits_) {
  
  print(paste0("performing computation for ", trait_))

  # get pheno_df
  df_ <- as.data.frame(fread(pheno_file_path_))
  vars_to_keep <- c(vars_to_keep_, trait_)
  df_ <- df_[, vars_to_keep]

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

  # compute individual-location clonal mean h2 for raw phenotypes
  env_h2_pheno_list_ <- rep(0, length(env_list_))
  names(env_h2_pheno_list_) <- env_list_
  for (env_ in env_list_)
  {
    df_trait_env_ <- df_[df_$Envir == env_, c("Genotype", "Envir", trait_)]

    # drop na for trait_
    df_trait_env_ <- df_trait_env_ %>% drop_na(all_of(trait_))

    if (nrow(df_trait_env_) > min_obs_lmer_) {
      # get average number of replications for genotypes
      nr_bar_ <- mean(table(df_trait_env_$Genotype))
      try(
        {
          # compute lmer mixed model  y = mu + g + eps with g and eps as random factors
          trait_lmer_mod_ <- lmer(as.formula(paste0(trait_, " ~ 1 + (1 | Genotype)")),
            data = df_trait_env_
          )
          # compute individual clonal mean h2 for each environment
          env_h2_pheno_list_[[env_]] <- compute_indiv_location_clonal_mean_h2(
            trait_lmer_mod_, nr_bar_
          )
        },
        silent = TRUE
      )
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
  env_h2_adj_pheno_list_ <- rep(0, length(env_list_))
  names(env_h2_adj_pheno_list_) <- env_list_
  for (env_ in env_list_)
  {
    df_trait_env_ <- df_[df_$Envir == env_, c("Genotype", "Envir", "spats_adj_pheno")]

    if (nrow(df_trait_env_) > min_obs_lmer_) {
      try(
        {
          # get average number of replications for genotypes
          nr_bar_ <- mean(table(df_trait_env_$Genotype))

          # compute lmer mixed model  y = mu + g + eps with g and eps as random factors
          trait_lmer_mod_ <- lmer(spats_adj_pheno ~ 1 + (1 | Genotype),
            data = df_trait_env_
          )
          # compute individual clonal mean h2 for each environment
          env_h2_adj_pheno_list_[[env_]] <- compute_indiv_location_clonal_mean_h2(
            trait_lmer_mod_, nr_bar_
          )
        },
        silent = TRUE
      )
    }
  }
  
  # keep envir with h2 different from zero (i.e. envir with no data)
  env_h2_adj_pheno_list_ <- env_h2_adj_pheno_list_[env_h2_adj_pheno_list_>0]
  
  # apply mad to detect outlier(s)
  h2_mad_value <- mad(env_h2_adj_pheno_list_, constant = 1)
  h2_median <- median(env_h2_adj_pheno_list_)
  h2_threshold <- h2_mad_value_factor * h2_mad_value
  
  # detect environment(s) to be excluded from adjusted pheno computed h2
  idx_env_to_exclude_h2_adj_pheno_ <-  which(abs(env_h2_adj_pheno_list_ - 
                                        h2_median) > h2_threshold &
                                      env_h2_adj_pheno_list_ < h2_median)
  
  env_to_exclude_h2_adj_pheno_ <- names(env_h2_adj_pheno_list_)[
    idx_env_to_exclude_h2_adj_pheno_
  ]

  # remove the environments to be excluded from the lists of heritabilities
  if (length(idx_env_to_exclude_h2_adj_pheno_) > 0) {
    env_h2_adj_pheno_list_ <- env_h2_adj_pheno_list_[-match(
      env_to_exclude_h2_adj_pheno_,
      names(env_h2_adj_pheno_list_)
    )]
  }

  # exclude identified environments by keeping only those with an h2 over
  # the defined threshold
  env_to_keep_ <- unique(names(env_h2_adj_pheno_list_))
  df_ <- df_[which(df_$Envir %in% env_to_keep_), ]
  
  # boxplots of heritabilities
  if (use_ggplot2_) {
    # convert to df
    df_h2_adj <- data.frame(
      Environment = names(env_h2_adj_pheno_list_),
      h2_adj = env_h2_adj_pheno_list_
    )
    df_h2 <- data.frame(
      Environment = names(env_h2_pheno_list_),
      h2 = env_h2_pheno_list_
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

  # define rename exceptions
  exception_cols <- c("Genotype", "Envir",
                      "Management", "Row", "Position",
                      "R", "P", trait_)
  
  # rename columns excluding the exception columns
  new_names <- colnames(df_)
  new_names[!(new_names %in% exception_cols)] <- paste0(trait_, "_", new_names[
    !(new_names %in% exception_cols)])
  
  # replace the existing column names with the new names
  colnames(df_) <- new_names
  
  # write long format data
  fwrite(df_, paste0(
    output_file_path, trait_,
    "_spats_adjusted_phenotypes_long_format.csv"
  ))
}
