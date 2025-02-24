# script meant to compute statistics and perform genotype x environment x management (GEM) analyses
# for the refpop. Note: text is formatted from Addins using Style active file from styler package

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
library(car)
library(stringr)
library(grid)
library(ggplot2)
library(ggcorrplot)
library(gridExtra)
library(reshape2)
library(lme4)
library(rstatix)
library(dplyr)
library(anytime)
library(foreach)
library(parallel)
library(doParallel)
library(future.apply)

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

h2_outlier_path <- "../../results/phenotype_outlier_detection/"

h2_file_path <- paste0(
  pheno_dir_path_,
  "h2_per_env_and_management/"
)

# output result path for phenotype graphics
output_gem_graphics_path <- "../../results/graphics/gem_graphics/"

# define function(s) and package(s) to export for parallelization
func_to_export_ <- c("fread")
pkgs_to_export_ <- c(
  "data.table", "stringr", "SpATS", "lme4",
  "lubridate", "emmeans", "plotly", "tidyr", "htmlwidgets"
)

# should lmer be refitted in case of singularity ?
refit_lmer_ <- F

# define vector of sites
vect_sites_ <- c("CHE", "BEL", "ITA", "FRA", "POL", "ESP")

# define selected_traits_ and vars_to_keep_ for output, note Weight_sample and
# Sample_size cannot be considered as a real trait
selected_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Scab_fruits",
  "Sample_size"
)
excluded_pseudo_trait_for_save_ <- c("Weight_sample", "Sample_size", "Scab_fruits")

# define analyses to be performed or not
perform_var_comp_analyses_ <- T
perform_signif_testing_gcym_ <- T
perform_signif_testing_genv_ <- T
parallelize_signif_testing_ <- T # should significance testing using LRT be parallelized ?

# define minimum number of observations for lmer
min_obs_lmer_ <- 5

# read phenotype data without md and knowledge rule based outliers
pheno_df_ <- as.data.frame(fread(paste0(
  pheno_dir_path_,
  "raw_phenotype_data_correc_manage_type_no_md_outliers.csv"
)))
pheno_df_ <- pheno_df_[pheno_df_$Management != "BUFFER", ]
colnames(pheno_df_)[match("Envir", colnames(pheno_df_))] <- "Environment"

# define traits
vect_traits <- setdiff(selected_traits_, excluded_pseudo_trait_for_save_)

# get environments, per trait, for which h2 are inliers
inlier_env_ <- as.data.frame(fread(paste0(
  h2_outlier_path,
  "envir_per_trait_retained_based_on_inlier_h2_for_adj_phenotypes.csv"
)))

# make graphic for mean pheno per management type in refpop
mean_pheno_per_manage_graphic_done_ <- create_mean_pheno_graphic_per_manage(
  pheno_df_,
  h2_outlier_path,
  selected_traits_,
  excluded_pseudo_trait_for_save_,
  output_gem_graphics_path
)
print(mean_pheno_per_manage_graphic_done_)

# make graphic for mean pheno per country type in refpop
mean_pheno_per_country_graphic_done_ <- create_mean_pheno_graphic_per_country(
  pheno_df_,
  h2_outlier_path,
  selected_traits_,
  excluded_pseudo_trait_for_save_,
  output_gem_graphics_path
)
print(mean_pheno_per_country_graphic_done_)

# make graphic for mean pheno per country type in refpop
mean_pheno_per_year_graphic_done_ <- create_mean_pheno_graphic_per_year(
  pheno_df_,
  h2_outlier_path,
  selected_traits_,
  excluded_pseudo_trait_for_save_,
  output_gem_graphics_path
)
print(mean_pheno_per_year_graphic_done_)

# variance components analyses
if (perform_var_comp_analyses_) {
  # estimate variance component contributions for the M_ge model per trait

  # define variance components for the model
  var_comp_names_ <- c("Environment", "Genotype", "Genotype:Environment")

  # get their contributions
  var_comp_contrib_per_trait <- estimate_var_comp_contrib_per_trait(
    pheno_df_,
    inlier_env_,
    var_comp_names_,
    vect_traits
  )$var_comp_contrib_per_trait

  # create required format for plot
  df_percent_var_contrib_per_trait <- do.call(cbind, var_comp_contrib_per_trait)
  df_plot <- melt(df_percent_var_contrib_per_trait,
    variable.name = "Variance component", value.name = "Percentage"
  )
  colnames(df_plot) <- c("Variance_component", "Trait", "Percentage_contribution")

  # make stacked barplots
  var_comp_genv_plot <- ggplot(df_plot, aes(
    x = Trait,
    y = Percentage_contribution, fill = Variance_component
  )) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(
      aes(label = ifelse(Percentage_contribution > 3,
        round(Percentage_contribution, 1), ""
      )),
      position = position_stack(vjust = 0.5),
      size = 3.5
    ) +
    labs(
      title = bquote("Percentage contribution of variance components in model" ~ M[g.env] ~
        "to total phenotypic variance for each trait in REFPOP (note: only the percentages superior to 3% are annotated)"),
      x = NULL,
      y = "Variance component percentage contribution (%)",
      fill = "Variance component"
    ) +
    scale_fill_manual(
      values = c(
        "Environment" = "coral3",
        "Genotype" = "green",
        "Genotype:Environment" = "gold4",
        "Residual" = "pink1"
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  # save plot
  ggsave(
    paste0(
      output_gem_graphics_path,
      "variance_decomposition_genv_refpop.png"
    ),
    plot = var_comp_genv_plot, width = 16, height = 8, dpi = 300
  )

  # estimate variance component contributions for the M_gcym model per trait

  # define variance components for the model
  var_comp_names_ <- c(
    "Country",
    "Management",
    "Year",
    "Genotype",
    "Country:Management",
    "Country:Year",
    "Management:Year",
    "Genotype:Country",
    "Genotype:Management",
    "Genotype:Year",
    "Genotype:Country:Management",
    "Genotype:Country:Year",
    "Country:Management:Year",
    "Genotype:Management:Year"
  )

  # get their contributions
  var_comp_contrib_per_trait <- estimate_var_comp_contrib_per_trait(
    pheno_df_,
    inlier_env_,
    var_comp_names_,
    vect_traits
  )$var_comp_contrib_per_trait

  # create required format for plot
  df_percent_var_contrib_per_trait <- do.call(cbind, var_comp_contrib_per_trait)
  df_plot <- melt(df_percent_var_contrib_per_trait,
    variable.name = "Variance component", value.name = "Percentage"
  )
  colnames(df_plot) <- c("Variance_component", "Trait", "Percentage_contribution")

  # make stacked barplots
  var_comp_gcym_plot <- ggplot(df_plot, aes(
    x = Trait, y =
      Percentage_contribution,
    fill = Variance_component
  )) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(
      aes(label = ifelse(Percentage_contribution > 3,
        round(Percentage_contribution, 1), ""
      )),
      position = position_stack(vjust = 0.5),
      size = 3.5
    ) +
    labs(
      title = bquote("Percentage contribution of variance components in model" ~ M[g.c.y.m] ~
        "to total phenotypic variance for each trait in REFPOP (note: only the percentages superior to 3% are annotated)"),
      x = NULL,
      y = "Variance component percentage contribution (%)",
      fill = "Variance component"
    ) +
    scale_fill_manual(
      values = c(
        "Country" = "red",
        "Management" = "skyblue",
        "Year" = "orange",
        "Genotype" = "green",
        "Country:Management" = "brown3",
        "Country:Year" = "#FF5722",
        "Management:Year" = "darkblue",
        "Genotype:Country" = "darkgreen",
        "Genotype:Management" = "grey",
        "Genotype:Year" = "#FFEB3B",
        "Genotype:Country:Management" = "#8E24AA",
        "Genotype:Country:Year" = "#046d7a",
        "Country:Management:Year" = "gold4",
        "Genotype:Management:Year" = "#8BC34A",
        "Residual" = "pink1"
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  # save plot
  ggsave(
    paste0(
      output_gem_graphics_path,
      "variance_decomposition_gcym_refpop.png"
    ),
    plot = var_comp_gcym_plot, width = 16, height = 8, dpi = 300
  )

  # estimate variance component contributions for the M_ymg model per trait

  # define countries
  vect_countries_ <- unique(pheno_df_$Country)

  for (country_ in vect_countries_) {
    print(country_)
    # define variance components for the model
    var_comp_names_ <- c(
      "Management",
      "Year",
      "Genotype",
      "Management:Year",
      "Genotype:Management",
      "Genotype:Year",
      "Genotype:Management:Year"
    )

    # define phenotype data frame for the specific country_
    pheno_df_country_ <- pheno_df_[pheno_df_$Country %in% country_, ]

    # get their contributions
    var_comp_contrib_per_trait <- estimate_var_comp_contrib_per_trait(
      pheno_df_ = pheno_df_country_,
      inlier_env_,
      var_comp_names_,
      vect_traits,
    )$var_comp_contrib_per_trait

    # create required format for plot
    df_percent_var_contrib_per_trait <- do.call(cbind, var_comp_contrib_per_trait)
    df_plot <- melt(df_percent_var_contrib_per_trait,
      variable.name = "Variance component", value.name = "Percentage"
    )
    colnames(df_plot) <- c("Variance_component", "Trait", "Percentage_contribution")

    # make stacked barplots
    var_comp_ymg_plot <- ggplot(df_plot, aes(
      x = Trait, y =
        Percentage_contribution,
      fill = Variance_component
    )) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(
        aes(label = ifelse(Percentage_contribution > 3,
          round(Percentage_contribution, 1), ""
        )),
        position = position_stack(vjust = 0.5),
        size = 3.5
      ) +
      labs(
        title = bquote("Percentage contribution of variance components in model" ~ M[g.y.m] ~
          "to total phenotypic variance for each trait in" ~ .(country_) ~
          "(note: only the percentages superior to 3% are annotated)"),
        x = NULL,
        y = "Variance component percentage contribution (%)",
        fill = "Variance component"
      ) +
      scale_fill_manual(
        values = c(
          "Management" = "skyblue",
          "Year" = "orange",
          "Genotype" = "green",
          "Management:Year" = "darkblue",
          "Genotype:Management" = "grey",
          "Genotype:Year" = "#FFEB3B",
          "Genotype:Management:Year" = "#8BC34A",
          "Residual" = "pink1"
        )
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    var_comp_ymg_plot

    # save plot
    ggsave(
      paste0(
        output_gem_graphics_path,
        paste0("variance_decomposition_ymg_", country_, ".png")
      ),
      plot = var_comp_ymg_plot, width = 16, height = 8, dpi = 300
    )
  }

  # estimate variance component contributions for the M_cmg model per trait

  # define countries
  vect_years_ <- unique(pheno_df_$Year)

  for (year_ in vect_years_) {
    print(year_)
    # define variance components for the model
    var_comp_names_ <- c(
      "Country",
      "Management",
      "Genotype",
      "Country:Management",
      "Genotype:Country",
      "Genotype:Management",
      "Genotype:Country:Management"
    )

    # define phenotype data frame for the specific year_
    pheno_df_country_ <- pheno_df_[pheno_df_$Year %in% year_, ]

    # get their contributions
    var_comp_contrib_per_trait <- estimate_var_comp_contrib_per_trait(
      pheno_df_ = pheno_df_country_,
      inlier_env_,
      var_comp_names_,
      vect_traits,
    )$var_comp_contrib_per_trait

    # create required format for plot
    df_percent_var_contrib_per_trait <- do.call(cbind, var_comp_contrib_per_trait)
    df_plot <- melt(df_percent_var_contrib_per_trait,
      variable.name = "Variance component", value.name = "Percentage"
    )
    colnames(df_plot) <- c("Variance_component", "Trait", "Percentage_contribution")

    # make stacked barplots
    var_comp_cmg_plot <- ggplot(df_plot, aes(
      x = Trait, y =
        Percentage_contribution,
      fill = Variance_component
    )) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(
        aes(label = ifelse(Percentage_contribution > 3,
          round(Percentage_contribution, 1), ""
        )),
        position = position_stack(vjust = 0.5),
        size = 3.5
      ) +
      labs(
        title = bquote("Percentage contribution of variance components in model" ~ M[g.c.m] ~
          "to total phenotypic variance for each trait in" ~ .(year_) ~
          "(note: only the percentages superior to 3% are annotated)"),
        x = NULL,
        y = "Variance component percentage contribution (%)",
        fill = "Variance component"
      ) +
      scale_fill_manual(
        values = c(
          "Country" = "red",
          "Management" = "skyblue",
          "Genotype" = "green",
          "Country:Management" = "brown3",
          "Genotype:Country" = "darkgreen",
          "Genotype:Management" = "grey",
          "Genotype:Country:Management" = "#8E24AA",
          "Residual" = "pink1"
        )
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    var_comp_ymg_plot

    # save plot
    ggsave(
      paste0(
        output_gem_graphics_path,
        paste0("variance_decomposition_cmg_", year_, ".png")
      ),
      plot = var_comp_cmg_plot, width = 16, height = 8, dpi = 300
    )
  }
}


# get environments for which h2 are inliers, for phenotypic means,
# and create vector of variables to keep
inlier_env_ <- as.data.frame(fread(paste0(
  h2_outlier_path,
  "envir_per_trait_retained_based_on_inlier_h2_for_adj_phenotypes.csv"
)))
var_to_keep <- c()
for (i in 1:nrow(inlier_env_)) {
  var_to_keep <- c(
    var_to_keep,
    paste0(
      inlier_env_$trait[i], "_",
      unlist(str_split(
        inlier_env_$env_retained_based_on_inlier_h2[i],
        ", "
      ))
    )
  )
}

# filter selected traits from pheno_df_ for environments whose h2 are inliers
pheno_filtered <- pheno_df_ %>%
  select(Country, Management, Year, Genotype, Environment, all_of(vect_traits)) %>%
  pivot_longer(
    cols = -c(Country, Management, Year, Genotype, Environment),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))
pheno_filtered <- pheno_filtered[
  paste0(pheno_filtered$Trait, "_", pheno_filtered$Environment) %in% var_to_keep,
]
df_ <- pheno_filtered %>%
  group_by(Genotype, Environment, Trait) %>%
  mutate(Rep = row_number()) %>%
  pivot_wider(names_from = Trait, values_from = Value) %>%
  select(-Rep)

# testing the statistical significance of main variance components using a chi-square distribution
# with q degrees of freedom, where q represents the number of coefficients tested for nullity under H0.
# This serves as an approximation for the Likelihood Ratio Test (LRT), though the true distribution
# is known to be a mixture of chi-square distributions with unknown degrees of freedom.
if (perform_signif_testing_genv_) {
  if (parallelize_signif_testing_) {
    df_$Genotype <- as.factor(paste0("_", df_$Genotype))
    df_$Environment <- as.factor(paste0("_", df_$Environment))

    mod_vars_ <- c("Genotype", "Environment", "Genotype:Environment")

    # function to check if a variable has multiple level:
    has_multiple_levels <- function(var) {
      return(length(unique(df_[[var]])) > 1)
    }

    # function to process a single variable-trait combination
    process_var_trait <- function(var_, trait_, df_, mod_vars_) {
      df_trait_ <- df_[!is.na(df_[[trait_]]), ]

      genotype_counts <- table(df_trait_$Genotype)
      genotypes_to_keep <- names(genotype_counts[genotype_counts > 1])

      df_trait_ <- df_trait_[df_trait_$Genotype %in% genotypes_to_keep, ]
      df_trait_ <- droplevels(df_trait_)

      if (sum(str_detect(mod_vars_, ":")) == 0) {
        lin_mod_vars_ <- mod_vars_[sapply(
          mod_vars_,
          function(v) {
            length(unique(df_trait_[[v]])) > 1
          }
        )]
      } else {
        lin_mod_vars_ <- mod_vars_
      }

      if (var_ %in% lin_mod_vars_) {
        if (str_detect(var_, pattern = ":")) {
          lmer_mod_trait_full <- lmer(
            as.formula(paste0(
              trait_,
              " ~ 1 + ",
              paste0(setdiff(lin_mod_vars_, var_), collapse = " + "),
              " + (1|", var_, ")"
            )),
            data = df_trait_
          )

          lmer_mod_trait_null <- lm(
            as.formula(paste0(
              trait_,
              " ~ 1 + ",
              paste0(setdiff(lin_mod_vars_, var_), collapse = " + ")
            )),
            data = df_trait_
          )

          p_val_ <- as.numeric(anova(
            lmer_mod_trait_full,
            lmer_mod_trait_null
          )$`Pr(>Chisq)`[2])
        } else {
          lmer_mod_trait_full <- lmer(
            as.formula(paste0(
              trait_,
              " ~ 1 + ",
              paste0(setdiff(lin_mod_vars_, c(var_, "Genotype:Environment")), collapse = " + "),
              " + (1|", var_, ")"
            )),
            data = df_trait_
          )

          lmer_mod_trait_null <- lm(
            as.formula(paste0(
              trait_,
              " ~ 1 + ",
              paste0(setdiff(lin_mod_vars_, c(var_, "Genotype:Environment")), collapse = " + ")
            )),
            data = df_trait_
          )

          p_val_ <- as.numeric(anova(
            lmer_mod_trait_full,
            lmer_mod_trait_null
          )$`Pr(>Chisq)`[2])
        }
        return(format(p_val_, digits = 3, nsmall = 3))
      }
      return(NA)
    }

    # run parallel processing for all variable-trait combinations
    results_list <- future_lapply(mod_vars_, function(var_) {
      sapply(vect_traits, function(trait_) {
        process_var_trait(var_, trait_, df_, mod_vars_)
      })
    })

    # convert list output to a data frame
    df_result_ <- as.data.frame(do.call(cbind, results_list))
    rownames(df_result_) <- vect_traits
    colnames(df_result_) <- mod_vars_

    # reset back to sequential execution
    plan(sequential)
    print(df_result_)

    # transform data frame to plot p-values
    df_melted <- melt(as.matrix(df_result_), na.rm = TRUE)
    colnames(df_melted) <- c("trait", "variance_component", "p_value")
    df_melted$p_value <- as.numeric(df_melted$p_value)
    df_melted$log_p <- -log10(df_melted$p_value)

    # add "(ns)" to log_p values below 3
    df_melted$log_p_label <- ifelse(df_melted$log_p < 3,
      paste0(round(df_melted$log_p, 2), " (non-significant)"),
      round(df_melted$log_p, 2)
    )

    # plot with "(ns)" added for values below 3
    lrt_signif_plot <- ggplot(df_melted, aes(x = variance_component, y = trait, fill = log_p)) +
      geom_tile() +
      geom_text(aes(label = log_p_label), color = "black", size = 3) +
      scale_fill_gradientn(
        colours = c("cyan3", "orange", "red2"),
        values = scales::rescale(c(0, 1, 2, 5, max(df_melted$log_p, na.rm = TRUE))),
        name = "-log10(p)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title = bquote("Statistical significance of variance components associated to" ~ M[g.env] ~
          "model, for each trait in REFPOP, based on the likelihood-ratio test (threshold: -log10(0.001) = 3)"),
        x = "Variance component",
        y = NULL
      )

    # save plot
    ggsave(
      paste0(
        output_gem_graphics_path,
        "variance_comp_genv_model_lrt_signif_refpop.png"
      ),
      plot = lrt_signif_plot, width = 16, height = 8, dpi = 300
    )
  }
}

if (perform_signif_testing_gcym_) {
  if (parallelize_signif_testing_) {
    # enable parallelization using available cores
    plan(multisession, workers = parallel::detectCores() - 1)

    df_$Year <- as.factor(paste0("_", df_$Year))
    df_$Country <- as.factor(paste0("_", df_$Country))
    df_$Genotype <- as.factor(paste0("_", df_$Genotype))
    df_$Management <- as.factor(paste0("_", df_$Management))

    mod_vars_ <- c("Country", "Management", "Year", "Genotype")

    # function to check if a variable has multiple levels
    has_multiple_levels <- function(var) {
      return(length(unique(df_[[var]])) > 1)
    }

    # function to process a single variable-trait combination
    process_var_trait <- function(var_, trait_, df_, mod_vars_) {
      df_trait_ <- df_[!is.na(df_[[trait_]]), ]

      genotype_counts <- table(df_trait_$Genotype)
      genotypes_to_keep <- names(genotype_counts[genotype_counts > 1])

      df_trait_ <- df_trait_[df_trait_$Genotype %in% genotypes_to_keep, ]
      df_trait_ <- droplevels(df_trait_)

      lin_mod_vars_ <- mod_vars_[sapply(
        mod_vars_,
        function(v) {
          length(unique(df_trait_[[v]])) > 1
        }
      )]

      if (var_ %in% lin_mod_vars_) {
        lmer_mod_trait_full <- lmer(
          as.formula(paste0(
            trait_,
            " ~ 1 + ",
            paste0(setdiff(lin_mod_vars_, var_), collapse = " + "),
            " + (1|", var_, ")"
          )),
          data = df_trait_
        )

        lmer_mod_trait_null <- lm(
          as.formula(paste0(
            trait_,
            " ~ 1 + ",
            paste0(setdiff(lin_mod_vars_, var_), collapse = " + ")
          )),
          data = df_trait_
        )

        p_val_ <- as.numeric(anova(
          lmer_mod_trait_full,
          lmer_mod_trait_null
        )$`Pr(>Chisq)`[2])

        return(format(p_val_, digits = 3, nsmall = 3))
      }
      return(NA)
    }

    # run parallel processing for all variable-trait combinations
    results_list <- future_lapply(mod_vars_, function(var_) {
      sapply(vect_traits, function(trait_) {
        process_var_trait(var_, trait_, df_, mod_vars_)
      })
    })

    # convert list output to a data frame
    df_result_ <- as.data.frame(do.call(cbind, results_list))
    rownames(df_result_) <- vect_traits
    colnames(df_result_) <- mod_vars_

    # reset back to sequential execution
    plan(sequential)
    print(df_result_)
  } else {
    df_$Year <- as.factor(paste0("_", df_$Year))
    df_$Country <- as.factor(paste0("_", df_$Country))
    df_$Genotype <- as.factor(paste0("_", df_$Genotype))
    df_$Management <- as.factor(paste0("_", df_$Management))

    mod_vars_ <- c("Country", "Management", "Year", "Genotype")

    df_result_ <- data.frame(matrix(NA,
      nrow = length(vect_traits),
      ncol = length(mod_vars_)
    ))
    rownames(df_result_) <- vect_traits
    colnames(df_result_) <- mod_vars_

    for (var_ in mod_vars_) {
      print(paste0("var_: ", var_))

      for (trait_ in vect_traits) {
        print(paste0("trait_: ", trait_))

        df_trait_ <- df_[!is.na(df_[, trait_]), ]

        genotype_counts <- table(df_trait_$Genotype)
        genotypes_to_keep <- names(genotype_counts[genotype_counts > 1])

        df_trait_ <- df_trait_[df_trait_$Genotype %in% genotypes_to_keep, ]
        df_trait_ <- droplevels(df_trait_)
        lin_mod_vars_ <- mod_vars_[sapply(mod_vars_, has_multiple_levels)]


        print(var_ %in% lin_mod_vars_)
        if (var_ %in% lin_mod_vars_) {
          lmer_mod_trait_full <- lmer(
            as.formula(paste0(
              trait_,
              "~ 1 + ",
              paste0(setdiff(lin_mod_vars_, var_), collapse = " + "),
              " + (1|", var_, ")"
            )),
            data = df_trait_
          )

          lmer_mod_trait_null <- lm(
            as.formula(paste0(
              trait_,
              "~ 1 + ",
              paste0(setdiff(lin_mod_vars_, var_), collapse = " + ")
            )),
            data = df_trait_
          )

          p_val_ <- as.numeric(anova(
            lmer_mod_trait_full,
            lmer_mod_trait_null
          )$`Pr(>Chisq)`[2])

          df_result_[trait_, var_] <- format(p_val_, digits = 3, nsmall = 3)
        }
      }
    }
  }

  # transform data frame to plot p-values
  df_melted <- melt(as.matrix(df_result_), na.rm = TRUE)
  colnames(df_melted) <- c("trait", "variance_component", "p_value")
  df_melted$p_value <- as.numeric(df_melted$p_value)
  df_melted$log_p <- -log10(df_melted$p_value)

  # add "(ns)" to log_p values below 3
  df_melted$log_p_label <- ifelse(df_melted$log_p < 3,
    paste0(round(df_melted$log_p, 2), " (non-significant)"),
    round(df_melted$log_p, 2)
  )

  # plot with "(ns)" added for values below 3
  lrt_signif_plot <- ggplot(df_melted, aes(x = variance_component, y = trait, fill = log_p)) +
    geom_tile() +
    geom_text(aes(label = log_p_label), color = "black", size = 3) +
    scale_fill_gradientn(
      colours = c("cyan3", "orange", "red2"),
      values = scales::rescale(c(0, 1, 2, 5, max(df_melted$log_p, na.rm = TRUE))),
      name = "-log10(p)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = bquote("Statistical significance of variance components associated to" ~ M[g.c.y.m] ~
        "model, for each trait in REFPOP, based on the likelihood-ratio test (threshold: -log10(0.001) = 3)"),
      x = "Variance component",
      y = NULL
    )

  # save plot
  ggsave(
    paste0(
      output_gem_graphics_path,
      "variance_comp_gcym_no_inter_model_lrt_signif_refpop.png"
    ),
    plot = lrt_signif_plot, width = 16, height = 8, dpi = 300
  )
}

# make correlation plots between traits for refpop and for each country

# make corrplots for refpop

# scale the data for correlation plots as it can affect the real correlation values
# due to variation in scales (e.g. higher standard deviation in one variable due to 
# different scales)
df_traits_refpop_ <- as.data.frame(scale(apply(
  df_[, vect_traits],
  2, as.numeric
), center = T, scale = T))

# remove columns that are entirely NA
df_traits_refpop_ <- df_traits_refpop_ %>% select(where(~ any(!is.na(.))))

# compute the correlation matrix, handling NA values properly
cor_matrix <- cor(df_traits_refpop_, use = "pairwise.complete.obs", method = "pearson")

# plot the heatmap with squares, showing both upper and lower triangles
corrplot_refpop_ <- ggcorrplot(cor_matrix,
  method = "square", # use squares instead of circles
  type = "full", # show both upper and lower triangular parts
  lab = TRUE, # display correlation coefficients
  lab_size = 3, # set label size
  colors = c("blue", "white", "red"), # color gradient: negative (blue) to positive (red)
  title = "Pearson correlations between raw phenotypic values of traits (without outliers) in REFPOP.",
  tl.cex = 10
) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(), # remove axis titles
    strip.background = element_blank() # remove strip background
  )

# save correlation plot
ggsave(
  paste0(
    output_gem_graphics_path,
    "correlation_plot_of_traits_refpop.png"
  ),
  plot = corrplot_refpop_, width = 16, height = 8, dpi = 300
)

# make corrplots for each country
for (country_ in vect_sites_) {
  df_traits_country_ <- df_[df_$Country %in% country_, vect_traits]

  # scale the data for correlation plots as it can affect the real correlation values
  # due to variation in scales (e.g. higher standard deviation in one variable due to 
  # different scales)
  df_traits_country_ <- as.data.frame(scale(apply(
    df_traits_country_,
    2, as.numeric
  ), center = T, scale = T))
  
  # remove columns that are entirely NA
  df_traits_country_ <- df_traits_country_ %>% select(where(~ any(!is.na(.))))

  # compute the correlation matrix, handling NA values properly
  cor_matrix <- cor(df_traits_country_, use = "pairwise.complete.obs", method = "pearson")

  # plot the heatmap with squares, showing both upper and lower triangles
  corrplot_country_ <- ggcorrplot(cor_matrix,
    method = "square", # use squares instead of circles
    type = "full", # show both upper and lower triangular parts
    lab = TRUE, # display correlation coefficients
    lab_size = 3, # set label size
    colors = c("blue", "white", "red"), # color gradient: negative (blue) to positive (red)
    title = paste0("Pearson correlations between raw phenotypic values of traits (without outliers) in ", country_, "."),
    tl.cex = 10
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title = element_blank(), # remove axis titles
      strip.background = element_blank() # remove strip background
    )

  # save correlation plot
  ggsave(
    paste0(
      output_gem_graphics_path,
      "correlation_plot_of_traits_", country_, ".png"
    ),
    plot = corrplot_country_, width = 16, height = 8, dpi = 300
  )
}
