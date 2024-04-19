# functions which corrects for trait spatial heterogeneity for a site
spat_hetero_env_correct_trait <- function(trait_, envir_, df_,
                                                      min_obs_ = 5) {
  tryCatch(
    {
      df_envir_ <- df_ %>% subset(Envir == envir_)
      df_envir_ <- drop_na(df_envir_)
      if (nrow(df_envir_) > min_obs_) {
        df_envir_$R <- as.factor(df_envir_$Row)
        df_envir_$P <- as.factor(df_envir_$Position)
        df_envir_$Genotype <- as.factor(df_envir_$Genotype)
        spats_formula <- as.formula(~ PSANOVA(Position, Row))
        spats_model <- SpATS(
          response = trait_,
          spatial = spats_formula,
          genotype = "Genotype", genotype.as.random = TRUE, fixed = NULL,
          random = ~ R + P, data = df_envir_,
          control = controlSpATS()
        )
        spats_geno_pred <- predict(spats_model, which = "Genotype")[
          , c("Genotype", "predicted.values")
        ]
        colnames(spats_geno_pred) <- c("Genotype", "spats_geno_pred_values")
        df_envir_$spats_residuals <- spats_model$residuals
        df_envir_ <- merge(df_envir_, spats_geno_pred, by = "Genotype", all = TRUE)
        df_envir_$spats_adj_pheno <- df_envir_$spats_geno_pred_values +
          df_envir_$spats_residuals
        return(df_envir_)
      } else {
        return(data.frame())
      }
    },
    error = function(e) {
      cat(
        "Error with spat_hetero_env_correct_trait,
      here is the possible issue with data and/or computation : ",
        conditionMessage(e), "\n"
      )
    }
  )
}


# function to re-order column names according to a specified order
reordered_cols <- function(col_names, prefix_order_patterns) {
  tryCatch(
    {
      # extract prefixes
      prefixes_ <- sub("_[0-9]{4}$", "", col_names)

      # order according to specified order
      ordered_countries <- factor(prefixes_, levels = prefix_order_patterns)

      # order the names
      ordered_cols <- col_names[order(ordered_countries, col_names)]

      return(ordered_cols)
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function to detect dates (does not work for all dates format, be careful)
contains_date <- function(string) {
  tryCatch(
    {
      !is.na(anydate(string))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function which computes individual location clonal mean heritability
compute_indiv_location_clonal_mean_h2 <- function(lmer_mod_, nr_bar_) {
  tryCatch(
    {
      var_cov <- VarCorr(lmer_mod_)
      sigma2G <- as.numeric(attr(var_cov$Genotype, "stddev")^2)
      sigma2E <- as.numeric(sigma(lmer_mod_)^2)
      h2 <- sigma2G / (sigma2G + sigma2E / nr_bar_)
      return(h2)
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function which multi-location clonal mean heritability (WIP)
compute_multi_location_clonal_mean_h2 <- function(lmer_mod_, nr_bar_, nl) {
  tryCatch(
    {
      var_cov <- VarCorr(lmer_mod_)
      sigma2G <- as.numeric(attr(var_cov$Genotype, "stddev")^2)
      sigma2Gl <- as.numeric(attr(var_cov$`Genotype:Envir`, "stddev")^2)
      sigma2E <- as.numeric(sigma(lmer_mod_)^2)
      h2 <- sigma2G / (sigma2G + sigma2Gl/nl + sigma2E/(nr_bar_*nl))
      return(h2)
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# (DEPRECATED) functions which corrects for trait spatial heterogeneity in raw_data
correct_trait_spatial_heterogeneity_raw_data <- function(
    trait_, geno_df_, df_year_, raw_data_site_names_, year_) {
  tryCatch(
    {
      for (site_name_ in raw_data_site_names_) {
        df_site_ <- df_year_ %>% subset(Country == site_name_)
        df_site_ <- df_site_[, -match("Country", colnames(df_site_))]
        df_site_ <- drop_na(df_site_)
        df_site_$R <- as.factor(df_site_$Row)
        df_site_$P <- as.factor(df_site_$Position)
        df_site_$Genotype <- as.factor(df_site_$Genotype)
        spats_formula <- as.formula(~ PSANOVA(Position, Row))
        spats_model <- SpATS(
          response = trait_,
          spatial = spats_formula,
          genotype = "Genotype", genotype.as.random = TRUE, fixed = NULL,
          random = ~ R + P, data = df_site_,
          control = controlSpATS()
        )
        spats_pred <- predict(spats_model, which = "Genotype")
        spats_pred <- spats_pred[, c("Genotype", "predicted.values")]
        colnames(spats_pred)[2] <- paste0(site_name_, "_", year_)
        spats_pred$Genotype <- as.character(spats_pred$Genotype)
        geno_df_ <- merge(geno_df_, spats_pred, by = "Genotype", all = TRUE)
      }
      return(geno_df_)
    },
    error = function(e) {
      cat(
        "Error with correct_trait_spatial_heterogeneity_raw_data,
      here is the possible issue with data and/or computation : ",
        conditionMessage(e), "\n"
      )
    }
  )
}

# (DEPRECATED) functions which corrects for trait spatial heterogeneity in raw_data
correct_trait_spatial_heterogeneity_munq_data <- function(
    trait_, geno_df_, list_df_year_, data_site_names_year_) {
  tryCatch(
    {
      for (site_name_ in data_site_names_year_) {
        print(site_name_)
        df_site_ <- list_df_year_[[site_name_]]
        df_site_[is.na(df_site_$MUNQ), "MUNQ"] <- df_site_[is.na(df_site_$MUNQ), "Genotype"]
        df_site_ <- df_site_[, -match("Genotype", colnames(df_site_))]
        df_site_ <- drop_na(df_site_)
        df_site_$R <- as.factor(df_site_$Row)
        df_site_$P <- as.factor(df_site_$Position)
        df_site_$MUNQ <- as.factor(df_site_$MUNQ)
        spats_formula <- as.formula(~ PSANOVA(Position, Row))
        spats_model <- SpATS(
          response = trait_,
          spatial = spats_formula,
          genotype = "MUNQ", genotype.as.random = TRUE, fixed = NULL,
          random = ~ R + P, data = df_site_,
          control = controlSpATS()
        )
        spats_pred <- predict(spats_model, which = "MUNQ")
        spats_pred <- spats_pred[, c("MUNQ", "predicted.values")]
        spats_pred$MUNQ <- as.character(spats_pred$MUNQ)
        colnames(spats_pred) <- c("Genotype", site_name_)
        geno_df_ <- merge(geno_df_, spats_pred, by = "Genotype", all = TRUE)
      }
      return(geno_df_)
    },
    error = function(e) {
      cat(
        "Error with correct_trait_spatial_heterogeneity_raw_data,
      here is the possible issue with data and/or computation : ",
        conditionMessage(e), "\n"
      )
    }
  )
}

