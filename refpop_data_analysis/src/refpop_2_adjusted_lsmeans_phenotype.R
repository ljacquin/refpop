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
library(tidyr)
library(lsmeans)

# set options to increase memory and suppress warnings
options(expressions = 5e5)
options(warn = -1)
emm_options(rg.limit = 20000)

# umap parameters, most sensitive ones
random_state_umap_ <- 15 # 15, 30 and 50
n_neighbors_umap_ <- 15 # 15, 30, 50
min_dist_ <- 0.1

# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../../functions.R")

# set paths
pheno_dir_path <- "../../data/phenotype_data/"
output_pheno_graphics_path <- "../../data/graphics/phenotype_graphics/"
spats_adj_pheno_path <- paste0(pheno_dir_path, "spats_per_env_adjusted_phenotypes/")

# define selected variables
vars_to_keep_ <- c("Envir", "Management", "Genotype")
excluded_pseudo_trait_for_save_ <- "Sample_size"

# threshold for removing columns with too much na
col_na_thresh_ <- 0.3

# set color palette for genotype families
color_palette_family <- c(
  "black",
  "lightblue",
  colorRampPalette(c("blue", "deepskyblue"))(10),
  colorRampPalette(c("orange", "orange2"))(2),
  colorRampPalette(c("aquamarine"))(1),
  colorRampPalette(c("magenta", "magenta2"))(2),
  colorRampPalette(c("firebrick3", "firebrick4"))(2),
  colorRampPalette(c("green", "darkgreen"))(5),
  colorRampPalette(c("darkorchid4"))(1),
  colorRampPalette(c("gold1", "gold2"))(2),
  colorRampPalette(c("lightcoral"))(1)
)

# set new color palette for genotype origins
color_palette_origin <- c(
  "magenta",
  "lightblue",
  "blue",
  "orange",
  "aquamarine",
  "red",
  "green",
  "darkorchid4",
  "gold2",
  "lightcoral",
  "yellow"
)

# adjusted phenotype data analysis

# get file names for spats adjusted phenotypes and replace pattern 
# "_spats_adjusted_.*" with "" for trait names
files_names_spats_adj_pheno <- list.files(spats_adj_pheno_path)
trait_names_ <- str_replace_all(files_names_spats_adj_pheno,
  "_spats_adjusted_.*",
  replacement = ""
)

# intialize lists for multi env h2 and ls-means of genotypes for each trait
multi_env_clonal_h2_traits_ <- rep(0, length(trait_names_))
names(multi_env_clonal_h2_traits_) <- trait_names_

list_ls_means_adj_pheno_per_geno <- vector("list", length(trait_names_))
names(list_ls_means_adj_pheno_per_geno) <- trait_names_


for (file_ in files_names_spats_adj_pheno) {
  print(paste0("computation for file : ", file_))

  df_ <- as.data.frame(fread(paste0(spats_adj_pheno_path, file_)))
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
    )] <- compute_multi_location_clonal_mean_h2(lmer_model_, nr_bar_, nl)

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
pheno_df <- Reduce(
  function(x, y) merge(x, y, by = "Genotype", all = TRUE),
  list_ls_means_adj_pheno_per_geno
)

# convert merge object to data.frame
pheno_df <- as.data.frame(pheno_df)
na_count <- colSums(is.na(pheno_df))
idx_col_to_drop <- which(na_count / nrow(pheno_df) > col_na_thresh_)
if (length(idx_col_to_drop) > 0) {
  pheno_df <- pheno_df[, -idx_col_to_drop]
}
colnames(pheno_df)[str_detect(colnames(pheno_df), "_lsmean")] <-
  str_replace_all(colnames(pheno_df)[str_detect(colnames(pheno_df), "_lsmean")],
    pattern = "_lsmean", replacement = ""
  )
fwrite(pheno_df, file=paste0(pheno_dir_path,'adjusted_ls_means_phenotypes.csv'))


# compute correlation matrix
cor_matrix <- cor(na.omit(pheno_df[, -na.omit(match(
  c(
    "Genotype",
    excluded_pseudo_trait_for_save_
  ), colnames(pheno_df)
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
    xaxis = list(title = "Trait"),
    yaxis = list(title = "Heritability (h2)")
  )

saveWidget(bar_plot, file = paste0(
  output_pheno_graphics_path,
  "multi_location_clonal_h2_bar_plot_multi_location_clonal_mean_h2.html"
))

# get family and origin informations
geno_fam_orig_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "genotype_family_origin_information.csv"
)))
pheno_df <- merge(pheno_df, geno_fam_orig_df, by = "Genotype", all = TRUE)
pheno_df <- pheno_df %>% select(c(Genotype, Family, Origin), everything())
pheno_df <- drop_na(pheno_df, any_of("Family"))

# apply umap
pheno_umap <- pheno_df[, -match(
  c("Genotype", "Family", "Origin"),
  colnames(pheno_df)
)]
# make sure their is no missing data, replace few na with column means
pheno_umap <- apply(pheno_umap, 2, impute_mean)

# standardize the data also before umap
pheno_umap <- scale(pheno_umap, center = T, scale = T)

pheno_umap_2d <- data.frame(umap(pheno_umap,
  n_components = 2,
  random_state = random_state_umap_,
  n_neighbors = n_neighbors_umap_,
  min_dist = min_dist_
)[["layout"]])

# plot umap with family as label
pheno_umap_2d$label <- pheno_df$Family
labels_ <- unique(pheno_df$Family)
n_family <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_family[1:n_family]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP phenotype data"
label_title_ <- "Family (except accession)"

fig_x_y <- plot_ly(
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = fig_title_,
    xaxis = list(title = "first component"),
    yaxis = list(title = "second component")
  )
# regroup by label
for (label_ in unique(pheno_umap_2d$label)) {
  data_subset <- pheno_umap_2d[pheno_umap_2d$label == label_, ]
  fig_x_y <- fig_x_y %>%
    add_trace(
      data = data_subset,
      x = ~X1, y = ~X2,
      type = "scatter", mode = "markers",
      marker = list(color = color_labels_[label_]),
      name = label_
    )
}
fig_x_y <- fig_x_y %>% layout(
  legend = list(title = list(text = label_title_))
)
saveWidget(fig_x_y,
  file = paste0(
    output_pheno_graphics_path,
    "complete_phenotype_refpop_umap_2d_family_as_label.html"
  )
)

# plot umap with origin as label
pheno_umap_2d$label <- pheno_df$Origin
labels_ <- unique(pheno_df$Origin)
n_origin <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_origin[1:n_origin]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP phenotype data"
label_title_ <- "Origin"

fig_x_y <- plot_ly(
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = fig_title_,
    xaxis = list(title = "first component"),
    yaxis = list(title = "second component")
  )
# regroup by label
for (label_ in unique(pheno_umap_2d$label)) {
  data_subset <- pheno_umap_2d[pheno_umap_2d$label == label_, ]
  fig_x_y <- fig_x_y %>%
    add_trace(
      data = data_subset,
      x = ~X1, y = ~X2,
      type = "scatter", mode = "markers",
      marker = list(color = color_labels_[label_]),
      name = label_
    )
}
fig_x_y <- fig_x_y %>% layout(
  legend = list(title = list(text = label_title_))
)
saveWidget(fig_x_y,
  file = paste0(
    output_pheno_graphics_path,
    "complete_phenotype_refpop_umap_2d_origin_as_label.html"
  )
)
