# script meant to analyse refpop phenotypic data
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
library(devtools)
if ("refpop_env" %in% conda_list()$name) {
  use_condaenv("refpop_env")
}
install_other_requirements <- F
if (install_other_requirements) {
  install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install("M3C")
  install.packages("remotes")
  remotes::install_github("hemstrow/snpR")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
}
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(dplyr)
library(htmlwidgets)
library(stringr)
library(lme4)
library(tidyr)
library(lsmeans)
library(lme4)

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
emm_options(rg.limit = 50000)

# set maximum number of principal components to be tested using akaike
# information criterion
max_n_comp_ <- 10

# define number of snp for sampling
# according to a uniform distribution
snp_sample_size_ <- 50e3

# threshold for removing columns with too much na
col_na_thresh_ <- 0.3

# set paths
pheno_dir_path <- "../../data/phenotype_data/"

spats_adj_pheno_path <- paste0(
  pheno_dir_path,
  "spats_per_env_adjusted_phenotypes/"
)
genom_dir_path <- "../../data/genotype_data/"

# output result path for phenotype graphics
output_pheno_graphics_path <- "../../results/graphics/phenotype_graphics/"

# define selected variables
vars_to_keep_ <- c("Envir", "Management", "Genotype")
excluded_pseudo_trait_for_save_ <- "Sample_size"

# adjusted phenotype data analysis

# get file names for spats adjusted phenotypes and replace pattern
# "_spats_adjusted_.*" with "" for trait names
files_names_spats_adj_pheno <- list.files(spats_adj_pheno_path)
trait_names_ <- str_replace_all(files_names_spats_adj_pheno,
  "_spats_adjusted_.*",
  replacement = ""
)
# initialize lists for lsmeans and blups associated to all traits
list_ls_means_adj_pheno_per_geno <- vector("list", length(trait_names_))
names(list_ls_means_adj_pheno_per_geno) <- trait_names_
blup_list_ <- vector("list", length(trait_names_))
names(blup_list_) <- trait_names_

# initialize vector for aic values
aic_ <- rep(0, max_n_comp_)

# compute principal components to be used as fixed covariates for population
# structure correction
geno_df <- as.data.frame(fread(paste0(genom_dir_path, "genotype_data.csv")))

# sample markers, according to a uniform distribution using identical as for
# genomic prediction (in order for all things to be equal)
set.seed(123)
idx_snp_sample_size_ <- sample(
  2:ncol(geno_df),
  size = snp_sample_size_, replace = F
)
geno_df <- geno_df[, c(
  match("Genotype", colnames(geno_df)),
  idx_snp_sample_size_
)]
geno_names_ <- geno_df[, 1]
geno_df <- geno_df[, -1]
geno_pca <- mixOmics::pca(apply(geno_df, 2, as.numeric), ncomp = max_n_comp_)
pc_coord_df_ <- as.data.frame(geno_pca$variates)[, 1:max_n_comp_]
pc_var_names_ <- colnames(pc_coord_df_)
pc_coord_df_$Genotype <- geno_names_

# intialize lists for multi env h2 and ls-means of genotypes for each trait
multi_env_clonal_h2_traits_ <- rep(0, length(trait_names_))
names(multi_env_clonal_h2_traits_) <- trait_names_

# for each trait, compute blups and lsmeans across all environments
# for each genotype
for (file_ in files_names_spats_adj_pheno) {
  print(paste0("computation for file : ", file_))

  df_ <- as.data.frame(fread(paste0(spats_adj_pheno_path, file_)))
  df_ <- df_[, c(vars_to_keep_, colnames(df_)[str_detect(
    colnames(df_),
    "spats_adj_pheno"
  )])]
  Y <- colnames(df_)[str_detect(colnames(df_), "spats_adj_pheno")]
  # for isolated genotypes with single phenotypic values (i.e. non repeated),
  # repeat their unique phenotypic values one time to make blup computation
  # possible
  idx_non_repeat_geno <- which(table(df_$Genotype) < 2)
  if (length(idx_non_repeat_geno) > 0) {
    df_ <- rbind(df_, df_[idx_non_repeat_geno, ])
  }
  # merge principal components and adjusted phenotypes based on genotype key
  df_ <- merge(df_, pc_coord_df_,
    by = "Genotype", all = TRUE
  )
  # remove any NA for trait adjusted phenotype before blup and lsmeans computation
  df_ <- df_[!is.na(df_[, Y]), ]

  if (length(unique(df_$Envir)) > 1) {
    # compute blups for genotypes using a linear mixed model (LMM) which fits
    # pc as fixed effects inorder to account for population structure

    # compute aic values in order to select number of pcs
    for (n_comp_ in 1:max_n_comp_) {
      lmer_model_ <- lmer(
        as.formula(paste0(
          Y,
          " ~ 1 + Envir + ", paste(pc_var_names_[1:n_comp_],
            collapse = " + "
          ),
          " + (1 | Genotype)"
        )),
        data = df_
      )
      aic_[n_comp_] <- AIC(lmer_model_)
    }
    n_opt_comp_aic_ <- which.min(aic_)
    print(paste0("number of pc selected: ", n_opt_comp_aic_))

    # estimate model based on selected number of pcs which minimize aic
    lmer_model_ <- lmer(
      as.formula(paste0(
        Y,
        " ~ 1 + Envir + ", paste(pc_var_names_[1:n_opt_comp_aic_],
          collapse = " + "
        ),
        " + (1 | Genotype)"
      )),
      data = df_
    )
    blup_list_[[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )]] <- data.frame(
      "Genotype" = rownames(ranef(lmer_model_)$Genotype),
      "blup" = as.numeric(unlist(ranef(lmer_model_)$Genotype))
    )

    # compute adjusted ls-means for genotypes across environments
    lm_model <- lm(formula(paste0(Y, "~ 1 + Genotype + Envir")), data = df_)
    ls_means <- as.data.frame(
      lsmeans(lm_model, ~Genotype)
    )[, c("Genotype", "lsmean")]

    # if the model coefficients are not estimable due to collinearity, use lm_()
    # to compute the coefficients with a pseudo-inverse approach. Then, use
    # lsmeans_() for calculating least squares means, which is based on group_by()
    # to replace the original lsmeans() function, as the latter requires a standard lm object.
    if (sum(is.na(lm_model$coefficients)) > 0) {
      lm_model <- lm_(formula(paste0(Y, "~ 1 + Genotype + Envir")), data = df_)
      ls_means <- as.data.frame(
        lsmeans_(lm_model, df_)
      )
    }
    colnames(ls_means)[match("lsmean", colnames(ls_means))] <-
      paste0(
        str_replace_all(file_, "_spats_adjusted_.*", replacement = ""),
        "_lsmean"
      )

    list_ls_means_adj_pheno_per_geno[[
      str_replace_all(file_, "_spats_adjusted_.*",
        replacement = ""
      )
    ]] <- ls_means

    # compute multi-location clonal mean heritability
    lmer_model_multi_loc <- lmer(
      as.formula(paste0(
        Y,
        " ~ 1 + Envir + (1 | Genotype) + (1 | Genotype:Envir)"
      )),
      data = df_
    )
    nl <- length(unique(df_$Envir))

    df_group_by <- df_ %>%
      group_by(Genotype, Envir) %>%
      summarise(frequency = n()) %>%
      data.frame()

    nr_bar_ <- mean(df_group_by$frequency)

    multi_env_clonal_h2_traits_[str_replace_all(
      file_, "_spats_adjusted_.*",
      replacement = ""
    )] <- compute_multi_location_clonal_mean_h2(
      lmer_model_multi_loc,
      nr_bar_, nl
    )
  } else {
    # compute blups for genotypes using a linear mixed model (LMM) which fits
    # pc as fixed effects inorder to account for population structure

    # compute aic values in order to select number of pcs
    for (n_comp_ in 1:max_n_comp_) {
      lmer_model_ <- lmer(
        as.formula(paste0(
          Y, " ~ 1 + ", paste(pc_var_names_[1:n_comp_],
            collapse = " + "
          ),
          " + (1 | Genotype)"
        )),
        data = df_
      )
      aic_[n_comp_] <- AIC(lmer_model_)
    }
    n_opt_comp_aic_ <- which.min(aic_)
    print(paste0("number of pc selected: ", n_opt_comp_aic_))

    # estimate model based on selected number of pcs which minimize aic
    lmer_model_ <- lmer(
      as.formula(paste0(
        Y, " ~ 1 + ", paste(pc_var_names_[1:n_opt_comp_aic_],
          collapse = " + "
        ),
        " + (1 | Genotype)"
      )),
      data = df_
    )
    blup_list_[[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )]] <- data.frame(
      "Genotype" = rownames(ranef(lmer_model_)$Genotype),
      "blup" = as.numeric(unlist(ranef(lmer_model_)$Genotype))
    )

    # compute adjusted ls-means for genotypes for unique environment
    lm_model <- lm(formula(paste0(Y, "~ 1 + Genotype")), data = df_)
    ls_means <- as.data.frame(
      lsmeans(lm_model, ~Genotype)
    )[, c("Genotype", "lsmean")]

    # if the model coefficients are not estimable due to collinearity, use lm_()
    # to compute the coefficients with a pseudo-inverse approach. Then, utilize
    # lsmeans_() for calculating least squares means, which is based on group_by()
    # to replace the original lsmeans() function, as the latter requires a standard lm object.
    if (sum(is.na(lm_model$coefficients)) > 0) {
      lm_model <- lm_(formula(paste0(Y, "~ 1 + Genotype")), data = df_)
      ls_means <- as.data.frame(
        lsmeans_(lm_model, df_)
      )
    }
    colnames(ls_means)[match("lsmean", colnames(ls_means))] <-
      paste0(
        str_replace_all(file_, "_spats_adjusted_.*", replacement = ""),
        "_lsmean"
      )

    list_ls_means_adj_pheno_per_geno[[
      str_replace_all(file_, "_spats_adjusted_.*",
        replacement = ""
      )
    ]] <- ls_means

    # compute single-location clonal mean heritability for unique environment
    lmer_model_single_loc_ <- lmer(as.formula(paste0(
      Y,
      " ~ 1 + (1 | Genotype)"
    )), data = df_)
    nr_bar_ <- mean(table(df_$Genotype))

    multi_env_clonal_h2_traits_[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )] <- compute_indiv_location_clonal_mean_h2(
      lmer_model_single_loc_, nr_bar_
    )
  }
}

# reduce blup list
blup_df <- Reduce(
  function(x, y) {
    merge(x, y,
      by = "Genotype",
      all = T
    )
  },
  blup_list_
)
colnames(blup_df) <- c("Genotype", trait_names_)

# write blups
fwrite(blup_df, file = paste0(
  pheno_dir_path,
  "blup_phenotypes.csv"
))

# merge list of ls_means into a single data frame for genotypes
lsmean_df <- Reduce(
  function(x, y) merge(x, y, by = "Genotype", all = T),
  list_ls_means_adj_pheno_per_geno
)
colnames(lsmean_df) <- c("Genotype", trait_names_)

# remove columns with excess NA if any
na_count <- colSums(is.na(lsmean_df))
idx_col_to_drop <- which(na_count / nrow(lsmean_df) > col_na_thresh_)
if (length(idx_col_to_drop) > 0) {
  lsmean_df <- lsmean_df[, -idx_col_to_drop]
}

# write lsmeans
fwrite(lsmean_df,
  file = paste0(
    pheno_dir_path,
    "adjusted_ls_mean_phenotypes.csv"
  )
)

# compute correlation matrix
cor_matrix <- cor(na.omit(lsmean_df[, -na.omit(match(
  c(
    "Genotype",
    excluded_pseudo_trait_for_save_
  ), colnames(lsmean_df)
))]))

# create an interactive heatmap
heatmap_colors_ <- colorRampPalette(c("red", "black"))(100)

heatmap_pearson <- plot_ly(
  z = cor_matrix,
  x = colnames(cor_matrix),
  y = colnames(cor_matrix),
  type = "heatmap",
  colorscale = list(list(seq(0, 1, length.out = length(heatmap_colors_)), heatmap_colors_))
)

# add annotations
heatmap_pearson <- heatmap_pearson %>%
  layout(
    title = "Pearson correlation between adjusted phenotypic LS-means of traits, per genotype, across all environments",
    xaxis = list(title = ""),
    yaxis = list(title = "")
  )

# display the heatmap_pearson
saveWidget(heatmap_pearson, file = paste0(
  output_pheno_graphics_path,
  "pearson_cor_adj_pheno_ls_means_per_geno_all_env.html"
))

# create a bar chart for the multi-location clonal mean h2
idx_exclu_pseudo_trait_ <- match(
  excluded_pseudo_trait_for_save_,
  names(multi_env_clonal_h2_traits_)
)
if (length(idx_exclu_pseudo_trait_) > 0) {
  multi_env_clonal_h2_traits <- multi_env_clonal_h2_traits_[
    -idx_exclu_pseudo_trait_
  ]
}
traits <- names(multi_env_clonal_h2_traits)
values <- as.numeric(multi_env_clonal_h2_traits)

# create the scatter plot with diamond markers
scatter_plot <- plot_ly(
  x = traits,
  y = values,
  type = "scatter",
  mode = "markers",
  marker = list(
    symbol = "diamond",
    size = 10,
    color = "rgb(158,202,225)",
    line = list(
      color = "rgb(8,48,107)",
      width = 1.5
    )
  )
) %>%
  layout(
    title = "Multi-location clonal mean heritability (h2) computed from adjusted phenotypes",
    xaxis = list(title = "Trait"),
    yaxis = list(title = "Heritability (h2)", range = c(0, 1.0))
  )

saveWidget(scatter_plot, file = paste0(
  output_pheno_graphics_path,
  "multi_location_clonal_h2.html"
))
