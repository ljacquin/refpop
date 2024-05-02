# script meant to analyse refpop pedigree data
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(rstudioapi)
library(stringr)

# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../../functions.R")

# set paths
pedig_dir_path <- "../../data/pedigree_data/"
pheno_dir_path <- "../../data/phenotype_data/"
geno_dir_path <- "../../data/genotype_data/"
output_pedig_graphics_path <- "../../data/graphics/pedig_graphics/"

# umap parameters, most sensitive ones
random_state_umap_ <- 15 
min_dist_ <- 0.1
n_neighbors_umap_ <- 50 # Note : as for K-NN, the number of neighbors 
                        # must be increased for sparse data

# define selected_traits_ 
selected_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Weight_sample", "Sample_size"
)

# set color palette for families
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

# set new color palette for origin
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

# get geographical origins of genotypes
geno_fam_orig_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_family_origin_information.csv"
)))

# get pedigree data
pedig_df <- as.data.frame(fread(paste0(
  pedig_dir_path,
  "updated_pedigree.txt"
)))

# rename to a common key, i.e. genotype
colnames(pedig_df)[
  colnames(pedig_df) == "Index"
] <- "Genotype"

# merge pedig_df with geno_fam_orig_df
pedig_df <- merge(pedig_df, geno_fam_orig_df, by = "Genotype", all = TRUE)
pedig_df <- drop_na(pedig_df, any_of("Family"))
pedig_df <- pedig_df[match(unique(pedig_df$Genotype), pedig_df$Genotype), ]

# replace missing parent value by genotype value
pedig_df <- replace_missing_parent_by_genotype(pedig_df)
pedig_incid_mat <- create_pedig_incid_mat(pedig_df)

# umap pedigree plots
pedig_umap_2d <- data.frame(umap(
  pedig_incid_mat[
    , -match("Genotype", colnames(pedig_incid_mat))
  ],
  n_components = 2,
  random_state = random_state_umap_,
  n_neighbors = n_neighbors_umap_,
  min_dist = min_dist_
)[["layout"]])

# plot umap pedigree with family as label
pedig_umap_2d$label <- pedig_df$Family
labels_ <- unique(pedig_df$Family)
n_family <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_family[1:n_family]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP pedigree data"
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
for (label_ in unique(pedig_umap_2d$label)) {
  data_subset <- pedig_umap_2d[pedig_umap_2d$label == label_, ]
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
    output_pedig_graphics_path,
    "complete_pedigree_refpop_umap_2d_family_as_label.html"
  )
)

# plot umap pedigree with origin as label
pedig_umap_2d$label <- pedig_df$Origin
labels_ <- unique(pedig_df$Origin)
n_origin <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_origin[1:n_origin]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP pedigree data"
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
for (label_ in unique(pedig_umap_2d$label)) {
  data_subset <- pedig_umap_2d[pedig_umap_2d$label == label_, ]
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
    output_pedig_graphics_path,
    "complete_pedigree_refpop_umap_2d_origin_as_label.html"
  )
)

# apply umap to pedigree and phenotype data
pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "adjusted_ls_means_phenotypes.csv"
)))
# center and scale data for umap and replace few na (13 in dataframe)
pheno_df[,selected_traits_] <- apply(scale(pheno_df[,selected_traits_], center = T,
                                     scale = T),2, impute_mean)
pedig_pheno_df <- merge(pheno_df, pedig_incid_mat, by = 'Genotype', all = T)

# umap pedigree phenotype plots
pedig_pheno_umap_2d <- data.frame(umap(
  pedig_pheno_df[
    , -match(c("Genotype", "Family", "Origin"), colnames(pedig_pheno_df))
  ],
  n_components = 2,
  random_state = random_state_umap_,
  n_neighbors = n_neighbors_umap_,
  min_dist = min_dist_
)[["layout"]])

# plot umap pedigree phenotype with family as label
pedig_pheno_umap_2d$label <- pedig_df$Family
labels_ <- unique(pedig_df$Family)
n_family <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_family[1:n_family]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP pedigree and phenotype data"
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
for (label_ in unique(pedig_pheno_umap_2d$label)) {
  data_subset <- pedig_pheno_umap_2d[pedig_pheno_umap_2d$label == label_, ]
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
             output_pedig_graphics_path,
             "complete_pedigree_phenotype_refpop_umap_2d_family_as_label.html"
           )
)

# plot umap pedigree phenotype with origin as label
pedig_pheno_umap_2d$label <- pedig_df$Origin
labels_ <- unique(pedig_df$Origin)
n_origin <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_origin[1:n_origin]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP pedigree and phenotype data"
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
for (label_ in unique(pedig_pheno_umap_2d$label)) {
  data_subset <- pedig_pheno_umap_2d[pedig_pheno_umap_2d$label == label_, ]
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
             output_pedig_graphics_path,
             "complete_pedigree_phenotype_refpop_umap_2d_origin_as_label.html"
           )
)

