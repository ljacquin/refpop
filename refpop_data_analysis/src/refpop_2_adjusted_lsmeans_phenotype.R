# script meant to analyse refpop phenotypic data
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
install_other_requirements <- FALSE
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
library(rstudioapi)
library(stringr)
library(lme4)
library(lsmeans)

# set options to increase memory and suppress warnings
options(expressions = 5e5)
options(warn = -1)
emm_options(rg.limit = 20000)

# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../../functions.R")

# set paths and parameters
path_spats_adj_pheno <- "../../data/phenotype_data/spats_per_env_adjusted_phenotypes/"
output_pheno_graphics_path <- "../../data/graphics/pheno_graphics/"

# define selected variables
vars_to_keep_ <- c("Envir", "Management", "Genotype")

# threshold for removing columns with too much na
threshold_na <- 0.3

# adjusted phenotype data analysis
files_names_spats_adj_pheno <- list.files(path_spats_adj_pheno)
trait_names_ <- str_replace_all(files_names_spats_adj_pheno,
  "_spats_adjusted_.*",
  replacement = ""
)

multi_env_clonal_h2_traits_ <- rep(0, length(trait_names_))
names(multi_env_clonal_h2_traits_) <- trait_names_

list_ls_means_adj_pheno_per_geno <- vector("list", length(trait_names_))
names(list_ls_means_adj_pheno_per_geno) <- trait_names_


for (file_ in files_names_spats_adj_pheno) {
  print(paste0("computation for file : ", file_))

  df_ <- as.data.frame(fread(paste0(path_spats_adj_pheno, file_)))
  df_ <- df_[, c(vars_to_keep_, colnames(df_)[str_detect(
    colnames(df_),
    "spats_adj_pheno"
  )])]
  Y <- colnames(df_)[str_detect(colnames(df_), "spats_adj_pheno")]

  if (length(unique(df_$Envir)) > 1) {
    # compute multi-location clonal mean heritability
    lmer_model_ <- lmer(
      as.formula(paste0(
        Y,
        " ~ 1 + (1 | Genotype) + Envir + (1 | Genotype:Envir)"
      )),
      data = df_
    )
    nl <- length(unique(df_$Envir))

    df_group_by <- df_ %>%
      group_by(Genotype, Envir) %>%
      summarise(frequency = n()) %>%
      data.frame()

    nr_bar_ <- mean(df_group_by$frequency)

    multi_env_clonal_h2_traits_[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )] <-
      compute_multi_location_clonal_mean_h2(lmer_model_, nr_bar_, nl)

    # compute adjusted ls-means for genotypes across environments
    lm_model <- lm(formula(paste0(Y, "~ Genotype + Envir")), data = df_)
    ls_means <- as.data.frame(lsmeans(lm_model, ~Genotype))[, c("Genotype", "lsmean")]
    colnames(ls_means)[match("lsmean", colnames(ls_means))] <- paste0(
      str_replace_all(file_, "_spats_adjusted_.*", replacement = ""), "_lsmean"
    )

    list_ls_means_adj_pheno_per_geno[[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )]] <- ls_means
  } else {
    # compute multi-location clonal mean heritability for unique environment
    lmer_model_ <- lmer(as.formula(paste0(
      Y,
      " ~ 1 + (1 | Genotype)"
    )), data = df_)
    nr_bar_ <- mean(table(df_$Genotype))

    multi_env_clonal_h2_traits_[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )] <- compute_indiv_location_clonal_mean_h2(
      lmer_model_, nr_bar_
    )

    # compute adjusted ls-means for genotypes for unique environment
    lm_model <- lm(formula(paste0(Y, "~ Genotype")), data = df_)
    ls_means <- as.data.frame(lsmeans(lm_model, ~Genotype))[, c("Genotype", "lsmean")]

    colnames(ls_means)[match("lsmean", colnames(ls_means))] <- paste0(
      str_replace_all(file_, "_spats_adjusted_.*", replacement = ""), "_lsmean"
    )

    list_ls_means_adj_pheno_per_geno[[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )]] <- ls_means
  }
}

#  merge list of ls_means into a single data frame for genotypes
merged_df <- Reduce(
  function(x, y) merge(x, y, by = "Genotype", all = TRUE),
  list_ls_means_adj_pheno_per_geno
)

# convert merge object to data.frame
merged_df <- as.data.frame(merged_df)
na_count <- colSums(is.na(merged_df))
filtered_df <- merged_df[, na_count / nrow(merged_df) <= threshold_na]
colnames(filtered_df)[str_detect(colnames(filtered_df), "_lsmean")] <-
  str_replace_all(colnames(filtered_df)[str_detect(colnames(filtered_df), "_lsmean")],
    pattern = "_lsmean", replacement = ""
  )

# compute correlation matrix
cor_matrix <- cor(filtered_df[, -1])

# create an interactif heatmap
my_colors <- colorRampPalette(c("red", "black"))(100)

heatmap_pearson <- plot_ly(
  z = cor_matrix,
  x = colnames(cor_matrix),
  y = colnames(cor_matrix),
  type = "heatmap",
  colorscale = list(list(seq(0, 1, length.out = length(my_colors)), my_colors))
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
traits <- names(multi_env_clonal_h2_traits_)
values <- as.numeric(multi_env_clonal_h2_traits_)

# create the bar chart
bar_plot <- plot_ly(
  x = traits,
  y = values,
  type = "bar",
  marker = list(
    color = "rgb(158,202,225)",
    line = list(
      color = "rgb(8,48,107)",
      width = 1.5
    )
  )
) %>%
  layout(
    title = "Multi-location clonal mean heritability (h2) computed from adjusted phenotypes",
    xaxis = list(title = "Traits"),
    yaxis = list(title = "Heritability (h2)")
  )

saveWidget(bar_plot, file = paste0(
  output_pheno_graphics_path,
  "multi_location_clonal_h2_bar_plot_multi_location_clonal_mean_h2.html"
))
