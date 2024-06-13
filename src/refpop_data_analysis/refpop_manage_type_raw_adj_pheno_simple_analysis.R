# script meant to analyse refpop pedigree data
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(data.table)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(rstudioapi)
library(stringr)
library(lme4)
# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../functions.R")

# set paths
pheno_dir_path <- "../../data/phenotype_data/"
# output result path for phenotype graphics
output_pheno_graphics_path <- "../../results/graphics/phenotype_graphics/"

# define analyzed trait_
trait_ <- "Harvest_date"

# get raw and adjusted phenotypes for analyzed trait_
df_raw <- droplevels(as.data.frame(fread(paste0(
  pheno_dir_path,
  "phenotype_raw_data_no_outliers.csv"
)))[, c("Genotype", "Management", "Envir", trait_)])
table(df_raw$Envir)
length(unique(df_raw$Envir))

df_adj <- droplevels(as.data.frame(fread(paste0(
  pheno_dir_path,
  "spats_per_env_adjusted_phenotypes/",
  paste0(trait_, "_spats_adjusted_phenotypes_long_format.csv")
)))[, c(
  "Genotype", "Management", "Envir",
  paste0(trait_, "_spats_adj_pheno")
)])

# show number of environments for after drop_na
df_raw <- df_raw %>% drop_na(all_of(trait_))

# slice raw and adj phenotypes for management type
df_raw_type_1 <- df_raw[df_raw$Management == 1, ]
df_raw_type_2 <- df_raw[df_raw$Management == 2, ]
df_adj_type_1 <- df_adj[df_adj$Management == 1, ]
df_adj_type_2 <- df_adj[df_adj$Management == 2, ]

# raw data
print(paste0(
  "Number of environments for ", trait_,
  " raw data : ",
  length(unique(df_raw$Envir))
))
print(paste0(
  "Number of environments for ", trait_,
  " raw data associated to management type 1 : ",
  length(unique(df_raw_type_1$Envir))
))
print(paste0(
  "Number of environments for ", trait_,
  " raw data associated to management type 2 : ",
  length(unique(df_raw_type_2$Envir))
))

# adj data
print(paste0(
  "Number of environments for ", trait_,
  " adj data : ",
  length(unique(df_adj$Envir))
))
print(paste0(
  "Number of environments for ", trait_,
  " adj data associated to management type 1 : ",
  length(unique(df_adj_type_1$Envir))
))
print(paste0(
  "Number of environments for ", trait_,
  " adj data associated to management type 2 : ",
  length(unique(df_adj_type_2$Envir))
))

# get common environments and associated slices
common_env_ <- intersect(
  unique(df_raw$Envir),
  unique(df_adj$Envir)
)
df_raw <- df_raw[df_raw$Envir %in% common_env_, ]
df_adj <- df_adj[df_adj$Envir %in% common_env_, ]

# initialize data frame to save results
row_names <- c("h2", "sigma2G", "sigma2E", "nr_bar_")

df_raw_h2_var_comp <- as.data.frame(matrix(NA,
  nrow = 4, ncol = length(common_env_)
))
rownames(df_raw_h2_var_comp) <- row_names
colnames(df_raw_h2_var_comp) <- common_env_

df_adj_h2_var_comp <- as.data.frame(matrix(NA,
  nrow = 4, ncol = length(common_env_)
))
rownames(df_adj_h2_var_comp) <- row_names
colnames(df_adj_h2_var_comp) <- common_env_

for (env_ in common_env_) {
  # get sliced data frames per env_
  df_raw_env_ <- df_raw[df_raw$Envir == env_, ]
  df_adj_env_ <- df_adj[df_adj$Envir == env_, ]

  # compute mixed random effect model to compute heritabilities and variance components

  # raw data computations
  trait_lmer_mod_ <- lmer(as.formula(paste0(trait_, " ~ 1 + (1 | Genotype)")),
    data = df_raw_env_
  )

  nr_bar_ <- mean(table(df_raw_env_$Genotype))

  list_h2_var_comp_ <- compute_indiv_location_clonal_mean_h2_var_comp(
    trait_lmer_mod_, nr_bar_
  )
  df_raw_h2_var_comp["h2", env_] <- list_h2_var_comp_$h2
  df_raw_h2_var_comp["sigma2G", env_] <- list_h2_var_comp_$sigma2G
  df_raw_h2_var_comp["sigma2E", env_] <- list_h2_var_comp_$sigma2E
  df_raw_h2_var_comp["nr_bar_", env_] <- list_h2_var_comp_$nr_bar_

  # adjusted data computations
  trait_spats_lmer_mod_ <- lmer(
    as.formula(paste0(
      trait_,
      "_spats_adj_pheno ~ 1 + (1 | Genotype)"
    )),
    data = df_adj_env_
  )

  nr_bar_spats_ <- mean(table(df_adj_env_$Genotype))
  list_spats_h2_var_comp_ <- compute_indiv_location_clonal_mean_h2_var_comp(
    trait_spats_lmer_mod_, nr_bar_spats_
  )

  df_adj_h2_var_comp["h2", env_] <- list_spats_h2_var_comp_$h2
  df_adj_h2_var_comp["sigma2G", env_] <- list_spats_h2_var_comp_$sigma2G
  df_adj_h2_var_comp["sigma2E", env_] <- list_spats_h2_var_comp_$sigma2E
  df_adj_h2_var_comp["nr_bar_", env_] <- list_spats_h2_var_comp_$nr_bar_
}

df_raw_h2_var_comp[, 1:5]
df_adj_h2_var_comp[, 1:5]
