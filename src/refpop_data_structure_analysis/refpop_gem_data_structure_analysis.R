# script meant to gem structure analyses
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
library(devtools)
if ("refpop_env" %in% conda_list()$name) {
  use_condaenv("refpop_env")
}
library(bigsnpr)
install_other_requirements <- F
if (install_other_requirements) {
  install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install("mixOmicsTeam/mixOmics")
  install.packages("remotes")
  remotes::install_github("hemstrow/snpR")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
}
library(snpR)
library(mixOmics)
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(stringr)
library(foreach)
library(doParallel)

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
# detect the number of cores to use
num_cores <- detectCores()

# set paths and booleans
pheno_dir_path <- "../../data/phenotype_data/"
spats_adj_pheno_path <- paste0(pheno_dir_path, "spats_per_env_adjusted_phenotypes/")

# output result path for genotype graphics
output_gem_graphics_path <- "../../results/graphics/gem_graphics/"

# set vector of labels to use
vect_country_or_year_as_label <- c("country", "year")

# define umap logical, training and plot parameters

# should umap and/or pca be performed ?
perform_umap_ <- T
perform_pca_ <- T

# umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 15
min_dist_ <- 1e-6

# color function for countries and years

# color function for countries
color_func_countries <- colorRampPalette(
  c(
    "magenta",
    "blue",
    "orange",
    "red",
    "green",
    "darkorchid4"
  )
)

# color function for years
color_func_year <- colorRampPalette(
  c(
    "#E0AD8B",
    "#2D9047",
    "#D3208A",
    "#D36720",
    "#A7CD6C",
    "#6CC1CD"
  )
)

# get file names for spats adjusted phenotypes and replace pattern
# "_spats_adjusted_.*" with "" for trait names
files_names_spats_adj_pheno <- list.files(spats_adj_pheno_path)
trait_names_ <- str_replace_all(files_names_spats_adj_pheno,
  "_spats_adjusted_.*",
  replacement = ""
)

for (file_ in files_names_spats_adj_pheno) {
  print(file_)
  try(
    {
      df_ <- as.data.frame(fread(paste0(spats_adj_pheno_path, file_)))
      df_$env_manag <- paste0(df_$Envir, "_", df_$Management)
      trait_name_ <- paste0(str_replace_all(file_, "_spats_adjusted_.*",
        replacement = ""
      ), "_spats_adj_pheno")
      df_ <- df_ %>% dplyr::select(env_manag, everything())
      df_ <- df_[, c("env_manag", "Genotype", trait_name_)]

      # transform data to wide format for trait_name_
      df_ <- df_ %>%
        pivot_wider(
          names_from = Genotype,
          values_from = trait_name_
        )
      df_ <- as.data.frame(df_)

      # remove env with buffer management type
      df_ <- df_[!str_detect(df_$env_manag, "BUFFER"), ]

      # compute cell means
      num_cols_ <- colnames(df_)[-1]
      for (col_ in num_cols_) {
        df_[[col_]] <- sapply(df_[[col_]], compute_cell_mean)
      }

      # split env_manag according to country, year and management type
      split_cols <- str_split(df_$env_manag, pattern = "_")

      # create new cols in data frame
      df_$country <- sapply(split_cols, `[`, 1)
      df_$year <- sapply(split_cols, `[`, 2)
      df_$management_type <- sapply(split_cols, `[`, 3)

      df_ <- df_ %>% dplyr::select(
        c(country, year, management_type),
        everything()
      )
      # delete management 2 and drop manangement column
      df_ <- df_[df_$management_type != "2", ]
      df_ <- df_[, !colnames(df_) %in% "management_type"]

      # remove cols with a rate of NA above threshold and impute those having a
      # maximum rate of threshold
      df_ <- remove_col_with_na_thresh(df_, threshold = 0.3)$filtered_df
      df_ <- as.data.frame(sapply(df_, impute_mean))

      # remove env_manag column
      df_ <- df_[, -match("env_manag", colnames(df_))]

      if (perform_umap_) {
        # compute umap in 2d and set vector of labels
        if (n_neighbors_umap_ < nrow(df_)) {
          umap_2d <- data.frame(umap(
            apply(df_[, -match(
              c("country", "year"),
              colnames(df_)
            )], 2, as.numeric),
            n_components = 2, random_state = random_state_umap_,
            n_neighbors = n_neighbors_umap_, min_dist = min_dist_
          )[["layout"]])
        } else {
          umap_2d <- data.frame(umap(
            apply(df_[, -match(
              c("country", "year"),
              colnames(df_)
            )], 2, as.numeric),
            n_components = 2, random_state = random_state_umap_,
            n_neighbors = floor(2 / 3 * nrow(df_)), min_dist = min_dist_
          )[["layout"]])
        }

        for (use_country_or_year_as_label in
          vect_country_or_year_as_label) {
          # umap plots

          # set palette of colors according to label used
          set.seed(123)

          # modify trait_name_ into trait_
          trait_ <- str_replace_all(trait_name_,
            pattern = "_spats_adj_pheno",
            replacement = ""
          )

          # set arbitrary colors for unique countries
          if (identical(use_country_or_year_as_label, "country")) {
            color_palette <- color_func_countries(length(unique(df_$country)))
            vect_labels <- df_$country
            label_title_ <- "Country"
            title_ <- paste0(
              "REFPOP UMAP 2D plot for mean adjusted phenotype of genotype x country associated to ",
              trait_
            )
            # set arbitrary colors for unique years
          } else {
            color_palette <- color_func_year(length(unique(df_$year)))
            vect_labels <- df_$year
            label_title_ <- "Year"
            title_ <- paste0(
              "REFPOP UMAP 2D plot for mean adjusted phenotype of genotype x year associated to ",
              trait_
            )
          }

          # get unique labels, their number, color palette and names,
          # and set vector of labels to umap
          labels_ <- unique(vect_labels)
          n_labels <- length(labels_)
          color_labels_ <- color_palette[1:n_labels]
          names(color_labels_) <- labels_
          umap_2d$label <- vect_labels

          # make plot with plot_ly
          fig_x_y <- plot_ly(
            type = "scatter", mode = "markers"
          ) %>%
            layout(
              plot_bgcolor = "#e5ecf6",
              title = title_,
              xaxis = list(title = "first component"),
              yaxis = list(title = "second component")
            )
          # regroup by label
          for (label_ in unique(umap_2d$label)) {
            data_subset <- umap_2d[umap_2d$label == label_, ]
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

          # test if folder exist for trait_
          if (!dir.exists(paste0(output_gem_graphics_path, "umap/", trait_))) {
            dir.create(paste0(output_gem_graphics_path, "umap/", trait_))
          }

          # save graphics
          saveWidget(fig_x_y,
            file = paste0(
              output_gem_graphics_path, "umap/",
              trait_, "/gem_umap_", tolower(trait_), "_",
              str_replace_all(tolower(label_title_),
                pattern = " ", replacement = "_"
              ),
              ".html"
            )
          )
        }
      }

      if (perform_pca_) {
        # compute pca in 2d and set vector of labels
        pca_obj_ <- pca(
          apply(df_[, -match(
            c("country", "year"),
            colnames(df_)
          )], 2, as.numeric),
          ncomp = 2, center = T, scale = T
        )
        pca_mat_ <- as.data.frame(pca_obj_$variates$X)
        pca_exp_var_ <- pca_obj_$prop_expl_var$X

        for (use_country_or_year_as_label in
          vect_country_or_year_as_label) {
          # umap plots

          # set palette of colors according to label used
          set.seed(123)

          # modify trait_name_ into trait_
          trait_ <- str_replace_all(trait_name_,
            pattern = "_spats_adj_pheno",
            replacement = ""
          )

          # set arbitrary colors for unique countries
          if (identical(use_country_or_year_as_label, "country")) {
            color_palette <- color_func_countries(length(unique(df_$country)))
            vect_labels <- df_$country
            label_title_ <- "Country"
            title_ <- paste0(
              "REFPOP PCA 2D plot for mean adjusted phenotype of genotype x country associated to ",
              trait_
            )
            # set arbitrary colors for unique years
          } else {
            color_palette <- color_func_year(length(unique(df_$year)))
            vect_labels <- df_$year
            label_title_ <- "Year"
            title_ <- paste0(
              "REFPOP PCA 2D plot for mean adjusted phenotype of genotype x year associated to ",
              trait_
            )
          }

          # get unique labels, their number, color palette and names,
          # and set vector of labels to umap
          labels_ <- unique(vect_labels)
          n_labels <- length(labels_)
          color_labels_ <- color_palette[1:n_labels]
          names(color_labels_) <- labels_
          pca_mat_$label <- vect_labels

          # make plot with plot_ly
          fig_x_y <- plot_ly(
            type = "scatter", mode = "markers"
          ) %>%
            layout(
              plot_bgcolor = "#e5ecf6",
              title = title_,
              xaxis = list(title = paste0(
                names(pca_exp_var_)[1], ": ",
                signif(100 * as.numeric(pca_exp_var_)[1], 2), "%"
              )),
              yaxis = list(title = paste0(
                names(pca_exp_var_)[2], ": ",
                signif(100 * as.numeric(pca_exp_var_)[2], 2), "%"
              ))
            )
          # regroup by label
          for (label_ in unique(pca_mat_$label)) {
            data_subset <- pca_mat_[pca_mat_$label == label_, ]
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
            legend = list(title = list(text = label_title_))
          )

          # test if folder exist for trait_
          if (!dir.exists(paste0(output_gem_graphics_path, "pca/", trait_))) {
            dir.create(paste0(output_gem_graphics_path, "pca/", trait_))
          }

          # save graphics
          saveWidget(fig_x_y,
            file = paste0(
              output_gem_graphics_path, "pca/",
              trait_, "/gem_pca_", tolower(trait_), "_",
              str_replace_all(tolower(label_title_),
                pattern = " ", replacement = "_"
              ),
              ".html"
            )
          )
        }
      }
    },
    silent = T
  )
}
