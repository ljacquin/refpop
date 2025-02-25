# script meant to perform genotype x environment x management (GEM) analyses
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
library(car)
library(stringr)
library(grid)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(lme4)
library(rstatix)
library(dplyr)
library(anytime)
library(foreach)
library(parallel)
library(doParallel)

# define computation mode, i.e. "local" or "cluster"
computation_mode <- "local"

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
excluded_pseudo_trait_for_save_ <- c("Weight_sample", "Sample_size")

# define minimum number of observations for lmer
min_obs_lmer_ <- 5

# read phenotype data without md and knowledge rule based outliers
pheno_df_ <- as.data.frame(fread(paste0(
  pheno_dir_path_,
  "raw_phenotype_data_correc_manage_type_no_md_outliers.csv"
)))
pheno_df_ <- pheno_df_[pheno_df_$Management != "BUFFER", ]

# make graphic for mean h2 per management across all environments in refpop
mean_h2_graphic_done_ <- create_mean_h2_graphic_per_manage_and_all_env(
  h2_file_path,
  selected_traits_,
  excluded_pseudo_trait_for_save_,
  vect_alpha_mean_h2 = c(0.01),
  output_gem_graphics_path
)
print(mean_h2_graphic_done_)

for (site_ in vect_sites_) {
  mean_h2_site_graphic_done_ <- create_mean_h2_graphic_per_manage_and_site(
    site_,
    h2_file_path,
    selected_traits_,
    excluded_pseudo_trait_for_save_,
    vect_alpha_mean_h2 = c(0.01),
    output_gem_graphics_path
  )
  print(paste0(
    "Graphic done for h2 associated to ", site_, ": ",
    as.character(mean_h2_site_graphic_done_)
  ))
}

# make graphic for mean pheno per management across all environments in refpop
mean_pheno_graphic_done_ <- create_mean_pheno_graphic_per_manage_and_all_env(
  pheno_df_,
  h2_outlier_path,
  selected_traits_,
  excluded_pseudo_trait_for_save_,
  vect_alpha_mean_y = c(0.01),
  output_gem_graphics_path
)
print(mean_pheno_graphic_done_)

for (site_ in vect_sites_) {
  mean_pheno_site_graphic_done_ <- create_mean_pheno_graphic_per_manage_and_site(
    site_,
    pheno_df_,
    h2_outlier_path,
    selected_traits_,
    excluded_pseudo_trait_for_save_,
    vect_alpha_mean_y = c(0.01),
    output_gem_graphics_path
  )
  print(paste0(
    "Graphic done for phenotypes in ", site_, ": ",
    as.character(mean_pheno_site_graphic_done_)
  ))
}


# get environments, per trait, for which h2 are inliers
inlier_env_ <- as.data.frame(fread(paste0(
  h2_outlier_path,
  "envir_per_trait_retained_based_on_inlier_h2_for_adj_phenotypes.csv"
)))

# variance components analyzes
vect_traits <- setdiff(selected_traits_, excluded_pseudo_trait_for_save_)

# compute variance components and their contributions for all components 
list_var_contrib <- list()
for (trait_ in vect_traits) {
  
  # specifiy inliers environments for trait_
  inlier_env_trait_ <- unlist(str_split(
    inlier_env_$env_retained_based_on_inlier_h2[
      inlier_env_$trait == trait_
    ], ", "
  ))
  pheno_df_trait_ <- pheno_df_[pheno_df_$Envir %in% inlier_env_trait_, ]
  pheno_df_trait_ <- pheno_df_trait_[!is.na(pheno_df_trait_[, trait_]), ]
  
  if (nrow(pheno_df_trait_) > min_obs_lmer_) {
    tryCatch(
      {
        # fit mixed model for variance components estimation
        lmer_mod_ <-
          lmer(
            as.formula(paste0(
              trait_,
              " ~ 1 + (1 | Management) + (1 | Country) + (1 | Year) + (1 | Genotype) +
              (1 | Country:Year) + (1 | Country:Management) +
              (1 | Genotype:Management) + (1 | Genotype:Country)"
            )),
            data = pheno_df_trait_
          )
        
        # is model singular 
        if ( isSingular(lmer_mod_) ){
          print(names(which(test<1e-8)))
        }
        
        # get variance components
        var_cov <- VarCorr(lmer_mod_)
        
        sigma2_management <- as.numeric(attr(var_cov$Management, "stddev")^2)
        sigma2_country <- as.numeric(attr(var_cov$Country, "stddev")^2)
        sigma2_year <- as.numeric(attr(var_cov$Year, "stddev")^2)
        sigma2_genotype <- as.numeric(attr(var_cov$Genotype, "stddev")^2)
        
        sigma2_country_year_inter <- as.numeric(attr(var_cov$`Country:Year`, "stddev")^2) /
          length(unique(pheno_df_trait_$Year))
        sigma2_country_manage_inter <- as.numeric(attr(var_cov$`Country:Management`, "stddev")^2) /
          length(unique(pheno_df_trait_$Management))
        
        sigma2_geno_country_inter <- as.numeric(attr(var_cov$`Genotype:Country`, "stddev")^2) /
          length(unique(pheno_df_trait_$Country))
        sigma2_geno_management_inter <- as.numeric(attr(var_cov$`Genotype:Management`, "stddev")^2) /
          length(unique(pheno_df_trait_$Management))

        
        # get residual variance
        df_group_by <- pheno_df_trait_ %>%
          group_by(Genotype, Envir) %>%
          summarise(frequency = n()) %>%
          data.frame()
        nr_bar_ <- mean(df_group_by$frequency)
        
        n_ <- nr_bar_ * length(unique(pheno_df_trait_$Country)) *
          length(unique(pheno_df_trait_$Management)) *
          length(unique(pheno_df_trait_$Year))
        
        sigma2_resid <- as.numeric(sigma(lmer_mod_)^2) / n_
        
        # get total variance
        var_y <- sigma2_country + sigma2_year + sigma2_management + sigma2_genotype +
          sigma2_country_year_inter + sigma2_country_manage_inter +
          sigma2_geno_country_inter + sigma2_geno_management_inter + 
          sigma2_resid
        
        # get contributions in percentage
        list_var_contrib[[trait_]] <- c(
          "Country" =
            100 * (sigma2_country / var_y),
          "Year" =
            100 * (sigma2_year / var_y),
          "Management" =
            100 * (sigma2_management / var_y),
          "Country x Year" =
            100 * (sigma2_country_year_inter / var_y),
          "Country x Management" =
            100 * (sigma2_country_manage_inter / var_y),  
          "Genotype" =
            100 * (sigma2_genotype / var_y),
          "Genotype x Country" =
            100 * (sigma2_geno_country_inter / var_y),
          "Genotype x Management" =
            100 * (sigma2_geno_management_inter / var_y),
          "Residual" =
            100 * (sigma2_resid / var_y)
        )
      },
      error = function(e) {
        cat(
          "Error with : ", conditionMessage(e), "\n"
        )
      }
    )
  } else {
    print(paste0("Not enough data for ", trait_, ": ", nrow(pheno_df_trait_)))
  }
}

# create required format for plot
df_percent_var_contrib_per_trait <- do.call(cbind, list_var_contrib)
df_plot <- melt(df_percent_var_contrib_per_trait,
                variable.name = "Variance component", value.name = "Percentage"
)
colnames(df_plot) <- c("Variance_component", "Trait", "Percentage_contribution")

# make stacked barplots
pf <- ggplot(df_plot, aes(x = Trait, y = Percentage_contribution, fill = Variance_component)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(Percentage_contribution > 5, round(Percentage_contribution, 1), "")),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  labs(
    title = "Percentage contribution of variance components to phenotypic variance for each trait in REFPOP
    (only the percentages superior to 5% are annotated)",
    x = "Trait",
    y = "Variance component percentage contribution (%)",
    fill = "Variance component"
  ) +
  scale_fill_manual(
    values = c(
      "Country" = "#1f78b4",
      "Year" = "#33a02c",
      "Management" = "#e31a1c",
      "Country x Year" = "#6a3d9a",
      "Country x Management" = "darkred",
      "Genotype" = "#ff7f00",
      "Genotype x Country" = "#b15928",
      "Genotype x Management" = "#a6cee3",
      "Residual" = "#fb9a99"
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
    "variance_decomposition_cymg_refpop_penalized.png"
  ),
  plot = pf, width = 16, height = 8, dpi = 300
)

# compute variance components and their contributions for genotype and environment
# components 
list_var_contrib <- list()
for (trait_ in vect_traits) {
  
  # specifiy inliers environments for trait_
  inlier_env_trait_ <- unlist(str_split(
    inlier_env_$env_retained_based_on_inlier_h2[
      inlier_env_$trait == trait_
    ], ", "
  ))
  pheno_df_trait_ <- pheno_df_[pheno_df_$Envir %in% inlier_env_trait_, ]
  pheno_df_trait_ <- pheno_df_trait_[!is.na(pheno_df_trait_[, trait_]), ]
  
  if (nrow(pheno_df_trait_) > min_obs_lmer_) {
    tryCatch(
      {
        # fit mixed model for variance components estimation
        lmer_mod_ <-
          lmer(
            as.formula(paste0(
              trait_,
              " ~ 1 + (1 | Genotype) + (1 | Envir) + (1 | Genotype:Envir)"
            )),
            data = pheno_df_trait_
          )
        
        # get variance components
        var_cov <- VarCorr(lmer_mod_)
        
        sigma2_genotype <- as.numeric(attr(var_cov$Genotype, "stddev")^2)
        sigma2_envir <- as.numeric(attr(var_cov$Envir, "stddev")^2)
        sigma2_geno_envir_inter <- as.numeric(attr(var_cov$`Genotype:Envir`, "stddev")^2) /
          length(unique(pheno_df_trait_$Envir))
        
        # get residual variance
        df_group_by <- pheno_df_trait_ %>%
          group_by(Genotype, Envir) %>%
          summarise(frequency = n()) %>%
          data.frame()
        nr_bar_ <- mean(df_group_by$frequency)
        
        n_ <- nr_bar_ * length(unique(pheno_df_trait_$Envir)) 
        
        sigma2_resid <- as.numeric(sigma(lmer_mod_)^2) / n_
        
        # get total variance
        var_y <- sigma2_genotype + sigma2_envir + sigma2_geno_envir_inter + 
          sigma2_resid
        
        # get contributions in percentage
        list_var_contrib[[trait_]] <- c(
          "Environment" =
            100 * (sigma2_envir / var_y),
          "Genotype" =
            100 * (sigma2_genotype / var_y),
          "Genotype x Environment" =
            100 * (sigma2_geno_envir_inter / var_y),
          "Residual" =
            100 * (sigma2_resid / var_y)
        )
      },
      error = function(e) {
        cat(
          "Error with : ", conditionMessage(e), "\n"
        )
      }
    )
  } else {
    print(paste0("Not enough data for ", trait_, ": ", nrow(pheno_df_trait_)))
  }
}

# create required format for plot
df_percent_var_contrib_per_trait <- do.call(cbind, list_var_contrib)
df_plot <- melt(df_percent_var_contrib_per_trait,
                variable.name = "Variance component", value.name = "Percentage"
)
colnames(df_plot) <- c("Variance_component", "Trait", "Percentage_contribution")

# make stacked barplots
pge <- ggplot(df_plot, aes(x = Trait, y = Percentage_contribution, fill = Variance_component)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(Percentage_contribution > 5, round(Percentage_contribution, 1), "")),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  labs(
    title = "Percentage contribution of variance components of Mge model to phenotypic variance for each trait in REFPOP
    (only the percentages superior to 5% are annotated)",
    x = "Trait",
    y = "Variance component percentage contribution (%)",
    fill = "Variance component"
  ) +
  scale_fill_manual(
    values = c(
      "Environment" = "darkturquoise",
      "Genotype" = "#ff7f00",
      "Genotype x Environment" = "darkgoldenrod",
      "Residual" = "#fb9a99"
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
pge

# save plot
ggsave(
  paste0(
    output_gem_graphics_path,
    "variance_decomposition_ge_refpop_penalized.png"
  ),
  plot = pge, width = 16, height = 8, dpi = 300
)
