# script meant to analyse refpop pedigree and phenotye data structure
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
library(devtools)
if ("refpop_env" %in% conda_list()$name) {
  use_condaenv("refpop_env")
}
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(stringr)
library(mixOmics)

# define computation mode, i.e. "local" or "cluster"
computation_mode <- 'cluster'

# if comutations are local in rstudio, detect and set script path
# automatically using rstudioapi
if ( identical(computation_mode, 'local') ){
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}

# source functions
source("../functions.R")

# set paths
pedig_dir_path <- "../../data/pedigree_data/"
pheno_dir_path <- "../../data/phenotype_data/"
geno_dir_path <- "../../data/genotype_data/"
# output result path for pedigree graphics
output_pedig_graphics_path <- "../../results/graphics/pedigree_graphics/"
# output result path for phenotype graphics
output_pheno_graphics_path <- "../../results/graphics/phenotype_graphics/"

# define selected_traits_
selected_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Weight_sample", "Sample_size"
)

# umap parameters, most sensitive ones
random_state_umap_ <- 15
min_dist_ <- 0.1
n_neighbors_umap_ <- 50 # Note : as for K-NN, the number of neighbors
# must be increased for sparse data

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
pedig_df <- merge(pedig_df, geno_fam_orig_df, by = "Genotype", all = T)
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
    "complete_refpop_pedigree_umap_2d_family_as_label.html"
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
    "complete_refpop_pedigree_umap_2d_origin_as_label.html"
  )
)

# apply umap to pedigree and breeding value data
bv_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "adjusted_ls_mean_breeding_values.csv"
)))

# center and scale data for umap and replace few na (13 in dataframe)
bv_df[, selected_traits_] <- apply(scale(bv_df[, selected_traits_],
  center = T,
  scale = T
), 2, impute_mean)


# merge bv_df and geno_fam_orig_df according to genotype
bv_df <- merge(bv_df, geno_fam_orig_df, by = "Genotype", all = TRUE)
bv_df <- drop_na(bv_df, any_of("Family"))

# remove non numeric variables and apply umap
bv_umap <- bv_df[, -match(
  c("Genotype", "Family", "Origin"),
  colnames(bv_df)
)]

bv_umap_2d <- data.frame(umap(bv_umap,
  n_components = 2,
  random_state = random_state_umap_,
  n_neighbors = n_neighbors_umap_,
  min_dist = min_dist_
)[["layout"]])

# plot umap with family as label
bv_umap_2d$label <- bv_df$Family
labels_ <- unique(bv_df$Family)
n_family <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_family[1:n_family]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP breeding value data"
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
for (label_ in unique(bv_umap_2d$label)) {
  data_subset <- bv_umap_2d[bv_umap_2d$label == label_, ]
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
    "complete_refpop_breeding_value_umap_2d_family_as_label.html"
  )
)

# plot umap with origin as label
bv_umap_2d$label <- bv_df$Origin
labels_ <- unique(bv_df$Origin)
n_origin <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_origin[1:n_origin]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP breeding value data"
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
for (label_ in unique(bv_umap_2d$label)) {
  data_subset <- bv_umap_2d[bv_umap_2d$label == label_, ]
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
    "complete_refpop_breeding_value_umap_2d_origin_as_label.html"
  )
)

# merge breeding value and pedig data
pedig_bv_df <- merge(bv_df, pedig_incid_mat, by = "Genotype", all = F)
pedig_bv_df <- pedig_bv_df[, -match(
    c("Genotype", "Family", "Origin"),
    colnames(pedig_bv_df)
  )
]
# umap pedigree breeding value plots
pedig_bv_umap_2d <- data.frame(umap(pedig_bv_df,
  n_components = 2,
  random_state = random_state_umap_,
  n_neighbors = n_neighbors_umap_,
  min_dist = min_dist_
)[["layout"]])

# plot umap pedigree breeding value with family as label
pedig_bv_umap_2d$label <- pedig_df$Family
labels_ <- unique(pedig_df$Family)
n_family <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_family[1:n_family]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP pedigree and breeding value data"
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
for (label_ in unique(pedig_bv_umap_2d$label)) {
  data_subset <- pedig_bv_umap_2d[pedig_bv_umap_2d$label == label_, ]
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
    "complete_refpop_pedigree_breeding_value_umap_2d_family_as_label.html"
  )
)

# plot umap pedigree breeding value with origin as label
pedig_bv_umap_2d$label <- pedig_df$Origin
labels_ <- unique(pedig_df$Origin)
n_origin <- length(labels_)

# define colors for labels
color_labels_ <- color_palette_origin[1:n_origin]
names(color_labels_) <- labels_

fig_title_ <- "UMAP 2D plot for REFPOP pedigree and breeding value data"
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
for (label_ in unique(pedig_bv_umap_2d$label)) {
  data_subset <- pedig_bv_umap_2d[pedig_bv_umap_2d$label == label_, ]
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
    "complete_refpop_pedigree_breeding_value_umap_2d_origin_as_label.html"
  )
)

# perform pca for pedig_incid_mat, bv_df and pedig_bv_df

# pedig_incid_mat
pedig_pca_obj_ <- pca(pedig_incid_mat,
  ncomp = 2, center = T, scale = T
)
pedig_pca_mat_ <- as.data.frame(pedig_pca_obj_$variates$X)
pedig_pca_exp_var_ <- pedig_pca_obj_$prop_expl_var$X

# plot coordinates of individuals on two first pcs :

# plot for family

# define colors for labels
labels_ <- unique(pedig_df$Family)
n_family <- length(labels_)
color_labels_ <- color_palette_family[1:n_family]
names(color_labels_) <- labels_
pedig_pca_mat_$label <- pedig_df$Family

# create plot
fig_x_y <- plot_ly(
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for REFPOP pedigree data",
    xaxis = list(title = paste0(
      names(pedig_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(pedig_pca_exp_var_)[1], 2), "%"
    )),
    yaxis = list(title = paste0(
      names(pedig_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(pedig_pca_exp_var_)[2], 2), "%"
    ))
  )
# regroup by label
for (label_ in unique(pedig_pca_mat_$label)) {
  data_subset <- pedig_pca_mat_[pedig_pca_mat_$label == label_, ]
  fig_x_y <- fig_x_y %>%
    add_trace(
      data = data_subset,
      x = ~PC1, y = ~PC2,
      type = "scatter", mode = "markers",
      marker = list(color = color_labels_[label_]),
      name = label_
    )
}
fig_x_y <- fig_x_y %>% layout(
  legend = list(title = list(text = "<b> Family (except accession) </b>"))
)
# save graphics
saveWidget(fig_x_y, file = paste0(
  output_pedig_graphics_path,
  "complete_refpop_pedigree_pca_family_as_label.html"
))

# plot for origin

# define colors for labels
labels_ <- unique(pedig_df$Origin)
n_origin <- length(labels_)
color_labels_ <- color_palette_origin[1:n_origin]
names(color_labels_) <- labels_
pedig_pca_mat_$label <- pedig_df$Origin

# create plot
fig_x_y <- plot_ly(
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for REFPOP pedigree data",
    xaxis = list(title = paste0(
      names(pedig_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(pedig_pca_exp_var_)[1], 2), "%"
    )),
    yaxis = list(title = paste0(
      names(pedig_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(pedig_pca_exp_var_)[2], 2), "%"
    ))
  )
# regroup by label
for (label_ in unique(pedig_pca_mat_$label)) {
  data_subset <- pedig_pca_mat_[pedig_pca_mat_$label == label_, ]
  fig_x_y <- fig_x_y %>%
    add_trace(
      data = data_subset,
      x = ~PC1, y = ~PC2,
      type = "scatter", mode = "markers",
      marker = list(color = color_labels_[label_]),
      name = label_
    )
}
fig_x_y <- fig_x_y %>% layout(
  legend = list(title = list(text = "<b> Origin </b>"))
)
# save graphics
saveWidget(fig_x_y, file = paste0(
  output_pedig_graphics_path,
  "complete_refpop_pedigree_pca_origin_as_label.html"
))

# bv_df
bv_df <- pedig_bv_df[, selected_traits_]
bv_pca_obj_ <- pca(bv_df,
  ncomp = 2, center = T, scale = T
)
bv_pca_mat_ <- as.data.frame(bv_pca_obj_$variates$X)
bv_pca_exp_var_ <- bv_pca_obj_$prop_expl_var$X

# plot coordinates of individuals on two first pcs :

# plot for family

# define colors for labels
labels_ <- unique(pedig_df$Family)
n_family <- length(labels_)
color_labels_ <- color_palette_family[1:n_family]
names(color_labels_) <- labels_
bv_pca_mat_$label <- pedig_df$Family

# create plot
fig_x_y <- plot_ly(
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for REFPOP breeding value data",
    xaxis = list(title = paste0(
      names(bv_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(bv_pca_exp_var_)[1], 4), "%"
    )),
    yaxis = list(title = paste0(
      names(bv_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(bv_pca_exp_var_)[2], 4), "%"
    ))
  )
# regroup by label
for (label_ in unique(bv_pca_mat_$label)) {
  data_subset <- bv_pca_mat_[bv_pca_mat_$label == label_, ]
  fig_x_y <- fig_x_y %>%
    add_trace(
      data = data_subset,
      x = ~PC1, y = ~PC2,
      type = "scatter", mode = "markers",
      marker = list(color = color_labels_[label_]),
      name = label_
    )
}
fig_x_y <- fig_x_y %>% layout(
  legend = list(title = list(text = "<b> Family (except accession) </b>"))
)
# save graphics
saveWidget(fig_x_y, file = paste0(
  output_pheno_graphics_path,
  "complete_refpop_breeding_value_pca_family_as_label.html"
))

# plot for origin

# define colors for labels
labels_ <- unique(pedig_df$Origin)
n_origin <- length(labels_)
color_labels_ <- color_palette_origin[1:n_origin]
names(color_labels_) <- labels_
bv_pca_mat_$label <- pedig_df$Origin

# create plot
fig_x_y <- plot_ly(
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for REFPOP breeding value data",
    xaxis = list(title = paste0(
      names(bv_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(bv_pca_exp_var_)[1], 4), "%"
    )),
    yaxis = list(title = paste0(
      names(bv_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(bv_pca_exp_var_)[2], 4), "%"
    ))
  )
# regroup by label
for (label_ in unique(bv_pca_mat_$label)) {
  data_subset <- bv_pca_mat_[bv_pca_mat_$label == label_, ]
  fig_x_y <- fig_x_y %>%
    add_trace(
      data = data_subset,
      x = ~PC1, y = ~PC2,
      type = "scatter", mode = "markers",
      marker = list(color = color_labels_[label_]),
      name = label_
    )
}
fig_x_y <- fig_x_y %>% layout(
  legend = list(title = list(text = "<b> Origin </b>"))
)
# save graphics
saveWidget(fig_x_y, file = paste0(
  output_pheno_graphics_path,
  "complete_refpop_breeding_value_pca_origin_as_label.html"
))

# pedig_bv_df
pedig_bv_pca_obj_ <- pca(
  pedig_bv_df,
  ncomp = 2, center = T, scale = T
)
pedig_bv_pca_mat_ <- as.data.frame(pedig_bv_pca_obj_$variates$X)
pedig_bv_pca_exp_var_ <- pedig_bv_pca_obj_$prop_expl_var$X

# plot coordinates of individuals on two first pcs :

# plot for family

# define colors for labels
labels_ <- unique(pedig_df$Family)
n_family <- length(labels_)
color_labels_ <- color_palette_family[1:n_family]
names(color_labels_) <- labels_
pedig_bv_pca_mat_$label <- pedig_df$Family

# create plot
fig_x_y <- plot_ly(
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for REFPOP pedigree and breeding value data",
    xaxis = list(title = paste0(
      names(pedig_bv_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(pedig_bv_pca_exp_var_)[1], 4), "%"
    )),
    yaxis = list(title = paste0(
      names(pedig_bv_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(pedig_bv_pca_exp_var_)[2], 4), "%"
    ))
  )
# regroup by label
for (label_ in unique(pedig_bv_pca_mat_$label)) {
  data_subset <- pedig_bv_pca_mat_[pedig_bv_pca_mat_$label == label_, ]
  fig_x_y <- fig_x_y %>%
    add_trace(
      data = data_subset,
      x = ~PC1, y = ~PC2,
      type = "scatter", mode = "markers",
      marker = list(color = color_labels_[label_]),
      name = label_
    )
}
fig_x_y <- fig_x_y %>% layout(
  legend = list(title = list(text = "<b> Family (except accession) </b>"))
)
# save graphics
saveWidget(fig_x_y, file = paste0(
  output_pedig_graphics_path,
  "complete_refpop_pedigree_breeding_value_pca_family_as_label.html"
))

# plot for origin

# define colors for labels
labels_ <- unique(pedig_df$Origin)
n_origin <- length(labels_)
color_labels_ <- color_palette_origin[1:n_origin]
names(color_labels_) <- labels_
pedig_bv_pca_mat_$label <- pedig_df$Origin

# create plot
fig_x_y <- plot_ly(
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for REFPOP pedigree and breeding value data",
    xaxis = list(title = paste0(
      names(pedig_bv_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(pedig_bv_pca_exp_var_)[1], 2), "%"
    )),
    yaxis = list(title = paste0(
      names(pedig_bv_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(pedig_bv_pca_exp_var_)[2], 2), "%"
    ))
  )
# regroup by label
for (label_ in unique(pedig_bv_pca_mat_$label)) {
  data_subset <- pedig_bv_pca_mat_[pedig_bv_pca_mat_$label == label_, ]
  fig_x_y <- fig_x_y %>%
    add_trace(
      data = data_subset,
      x = ~PC1, y = ~PC2,
      type = "scatter", mode = "markers",
      marker = list(color = color_labels_[label_]),
      name = label_
    )
}
fig_x_y <- fig_x_y %>% layout(
  legend = list(title = list(text = "<b> Origin </b>"))
)
# save graphics
saveWidget(fig_x_y, file = paste0(
  output_pedig_graphics_path,
  "complete_refpop_pedigree_breeding_value_pca_origin_as_label.html"
))
