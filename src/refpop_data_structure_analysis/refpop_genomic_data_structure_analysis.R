# script meant to analyse refpop genotypic data
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
library(htmlwidgets)
library(stringr)
library(foreach)
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
# detect the number of cores to use
num_cores <- detectCores()

# set paths and booleans
geno_dir_path <- "../../data/genotype_data/"
pheno_dir_path <- "../../data/phenotype_data/"
progeny_data_path <- "../../data/progeny_data/"
# output result path for genotype graphics
output_geno_graphics_path <- "../../results/graphics/genotype_graphics/"

bed_file <- paste0(geno_dir_path, "refpop_genotype.bed")
bim_file <- paste0(geno_dir_path, "refpop_genotype.bim")
fam_file <- paste0(geno_dir_path, "refpop_genotype.fam")

read_with_bigsnpr <- TRUE
perform_umap_ <- FALSE

# define umap training and plot paraemeters

# umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 15
min_dist_ <- 0.1

# define refpop train data for umap : "complete", "accessions", "progeny"
umap_refpop_train_data <- "complete"
pca_refpop_train_data <- "complete"

# define label data for umap :  origin, family or genotype (genotype not recommended)
vect_origin_family_or_genotype_as_label_ <- c("origin", "family")

# if umap_refpop_train_data = "accessions", should progenies be projected using umap ?
predict_umap_progeny_ <- FALSE

# color function and palettes for genotypes, families and origins

# color function for genotypes
col_func_genotype <- colorRampPalette(c(
  "black", "blue", "red", "orange",
  "yellow", "green"
))

# color palette for families (28 counts)
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

# color palette for origins (11 counts)
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

# read plink genotype data
if (read_with_bigsnpr) {
  # reading the bedfile for the first time or its associated rds file if exists
  # note .bed and .fam files should be in the same directory
  if (!file.exists(paste0(geno_dir_path, "refpop_genotype.rds"))) {
    df_ <- snp_readBed(bed_file)
    df_ <- readRDS(paste0(geno_dir_path, "refpop_genotype.rds"))
  } else {
    df_ <- readRDS(paste0(geno_dir_path, "refpop_genotype.rds"))
  }
  fam_df <- df_[["fam"]]
  map_df <- df_[["map"]]
  geno_df <- as.data.frame(df_[["genotypes"]][
    1:nrow(fam_df),
    1:nrow(map_df)
  ])
  rm(df_)
  colnames(geno_df) <- map_df$marker.ID
  geno_df$Genotype <- fam_df$sample.ID
  geno_df <- geno_df %>% dplyr::select(Genotype, everything())
  if (!file.exists(paste0(geno_dir_path, "genotype_data.csv"))) {
    fwrite(geno_df, paste0(geno_dir_path, "genotype_data.csv"))
  }
}

# get geographical origins of genotypes
geno_origin_ <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_geographical_origins.csv"
)))

# rename to a common key, i.e. genotype
colnames(geno_origin_)[
  colnames(geno_origin_) == "Genotype Code"
] <- "Genotype"

# get families of genotypes
geno_fam_ <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_family_patterns.csv"
)))

# detect and set families for genotypes
geno_origin_$Family <- NA
for (i in 1:length(geno_fam_$Pattern)) {
  idx_geno_fam <- which(str_detect(
    geno_origin_$Genotype,
    pattern = geno_fam_$Pattern[i]
  ))
  geno_origin_$Family[idx_geno_fam] <- geno_fam_$Family[i]
}
geno_origin_$Family[is.na(geno_origin_$Family)] <- "Accession"

# define geno_fam_orig_vect_
geno_fam_orig_vect_ <- c("Genotype", "Family", "Origin")

# merge geno_df with geno_origin
geno_df <- merge(geno_df, geno_origin_[, geno_fam_orig_vect_],
  by = "Genotype", all = TRUE
)
geno_df <- geno_df %>% dplyr::select(c(Genotype, Family, Origin), everything())
idx_origin_geno_names <- which(colnames(geno_df) %in% c(
  "Genotype",
  "Family",
  "Origin"
))

# useful : save genotype, family and origin information to genotype and phenotype
# data folders
fwrite(geno_df[, geno_fam_orig_vect_],
  file = paste0(geno_dir_path, "genotype_family_origin_information.csv"),
  sep = ","
)

fwrite(geno_df[, geno_fam_orig_vect_],
  file = paste0(pheno_dir_path, "genotype_family_origin_information.csv"),
  sep = ","
)

for (use_origin_family_or_genotype_as_label_ in
  vect_origin_family_or_genotype_as_label_) {
  # umap plots

  if (perform_umap_) {
    # set palette of colors according to label used
    set.seed(123)

    # set arbitrary colors for unique genotypes
    if (identical(use_origin_family_or_genotype_as_label_, "genotype")) {
      color_palette <- col_func_genotype(length(unique(geno_df$Genotype)))

      # set color labels for families (28 counts)
    } else if (identical(use_origin_family_or_genotype_as_label_, "family")) {
      color_palette <- color_palette_family

      # if population is not complete select the corresponding color for subset
      if (identical(umap_refpop_train_data, "accessions")) {
        color_palette <- color_palette[1]
      } else if (identical(umap_refpop_train_data, "progeny")) {
        color_palette <- color_palette[-1]
      }

      # set color labels for origins (11 counts)
    } else {
      color_palette <- color_palette_origin
    }

    if (identical(umap_refpop_train_data, "complete")) {
      # compute umap in 2D
      geno_umap_2d <- data.frame(umap(geno_df[, -idx_origin_geno_names],
        n_components = 2, random_state = random_state_umap_,
        n_neighbors = n_neighbors_umap_, min_dist = min_dist_
      )[["layout"]])

      # compute umap in 3D
      geno_umap_3d <- data.frame(umap(geno_df[, -idx_origin_geno_names],
        n_components = 3, random_state = random_state_umap_,
        n_neighbors = n_neighbors_umap_, min_dist = min_dist_
      )[["layout"]])

      # save umap models for this case if predict_umap_progeny_ is true
    } else if (identical(umap_refpop_train_data, "accessions")) {
      sub_geno_df <- geno_df[-which(geno_df$Origin %in% "P"), ]

      # compute umap in 2D
      geno_umap_2d_model <- umap(sub_geno_df[, -idx_origin_geno_names],
        n_components = 2, random_state = random_state_umap_,
        n_neighbors = n_neighbors_umap_, min_dist = min_dist_
      )
      geno_umap_2d <- data.frame(geno_umap_2d_model[["layout"]])

      # compute umap in 3D
      geno_umap_3d_model <- umap(sub_geno_df[, -idx_origin_geno_names],
        n_components = 3, random_state = random_state_umap_,
        n_neighbors = n_neighbors_umap_, min_dist = min_dist_
      )
      geno_umap_3d <- data.frame(geno_umap_3d_model[["layout"]])
    } else if (identical(umap_refpop_train_data, "progeny")) {
      sub_geno_df <- geno_df[which(geno_df$Origin %in% "P"), ]

      # compute umap in 2D
      geno_umap_2d <- data.frame(umap(sub_geno_df[, -idx_origin_geno_names],
        n_components = 2, random_state = random_state_umap_,
        n_neighbors = n_neighbors_umap_, min_dist = min_dist_
      )[["layout"]])
      # compute umap in 3D
      geno_umap_3d <- data.frame(umap(sub_geno_df[, -idx_origin_geno_names],
        n_components = 3, random_state = random_state_umap_,
        n_neighbors = n_neighbors_umap_, min_dist = min_dist_
      )[["layout"]])
    }

    if (identical(use_origin_family_or_genotype_as_label_, "origin")) {
      # define colors for labels
      labels_ <- unique(geno_df$Origin)
      n_origins <- length(labels_)
      color_labels_ <- color_palette[1:n_origins]
      names(color_labels_) <- labels_

      if (identical(umap_refpop_train_data, "complete")) {
        # define label according to origin
        geno_umap_2d$label <- geno_df$Origin
        geno_umap_3d$label <- geno_df$Origin
      } else {
        geno_umap_2d$label <- sub_geno_df$Origin
        geno_umap_3d$label <- sub_geno_df$Origin
      }
    } else if (identical(use_origin_family_or_genotype_as_label_, "family")) {
      # define colors for labels
      labels_ <- unique(geno_df$Family)
      n_family <- length(labels_)
      color_labels_ <- color_palette[1:n_family]
      names(color_labels_) <- labels_

      if (identical(umap_refpop_train_data, "complete")) {
        # define label according to origin
        geno_umap_2d$label <- geno_df$Family
        geno_umap_3d$label <- geno_df$Family
      } else {
        geno_umap_2d$label <- sub_geno_df$Family
        geno_umap_3d$label <- sub_geno_df$Family
      }
    } else if (identical(use_origin_family_or_genotype_as_label_, "genotype")) {
      # define colors for labels
      labels_ <- unique(geno_df$Genotype)
      n_geno <- length(labels_)
      color_labels_ <- color_palette[1:n_geno]
      names(color_labels_) <- labels_

      if (identical(umap_refpop_train_data, "complete")) {
        # define label according to origin
        geno_umap_2d$label <- geno_df$Genotype
        geno_umap_3d$label <- geno_df$Genotype
      } else {
        geno_umap_2d$label <- sub_geno_df$Genotype
        geno_umap_3d$label <- sub_geno_df$Genotype
      }
    }

    # make umap prediction for progenies based on unsupervised learning
    # from accessions as a special case
    if (identical(umap_refpop_train_data, "accessions") && predict_umap_progeny_) {
      sub_geno_df <- geno_df[which(geno_df$Origin %in% "P"), ]

      # 2d umap prediction for progenies
      geno_umap_2d_progeny <- as.data.frame(predict(
        geno_umap_2d_model, sub_geno_df[, -idx_origin_geno_names]
      ))
      if (identical(use_origin_family_or_genotype_as_label_, "family")) {
        geno_umap_2d_progeny$label <- sub_geno_df$Family
      } else if (identical(use_origin_family_or_genotype_as_label_, "origin")) {
        geno_umap_2d_progeny$label <- sub_geno_df$Origin
      }
      colnames(geno_umap_2d_progeny) <- colnames(geno_umap_2d)
      geno_umap_2d <- rbind(geno_umap_2d, geno_umap_2d_progeny)

      # 3d umap prediction for progenies
      geno_umap_3d_progeny <- as.data.frame(predict(
        geno_umap_3d_model, sub_geno_df[, -idx_origin_geno_names]
      ))
      if (identical(use_origin_family_or_genotype_as_label_, "family")) {
        geno_umap_3d_progeny$label <- sub_geno_df$Family
      } else if (identical(use_origin_family_or_genotype_as_label_, "origin")) {
        geno_umap_3d_progeny$label <- sub_geno_df$Origin
      }
      colnames(geno_umap_3d_progeny) <- colnames(geno_umap_3d)
      geno_umap_3d <- rbind(geno_umap_3d, geno_umap_3d_progeny)
    }

    # 2D plot
    # create base graphic
    if (identical(umap_refpop_train_data, "accessions") && predict_umap_progeny_) {
      umap_2d_title_ <- "UMAP 2D plot for REFPOP genotype data with umap trained
    on accessions, and progenies projected using trained model"
      output_path_2d_umap <- paste0(
        output_geno_graphics_path,
        umap_refpop_train_data, "/",
        umap_refpop_train_data, "_genotype_refpop_progeny_projected_umap_2d_",
        use_origin_family_or_genotype_as_label_,
        "_as_label.html"
      )
    } else {
      umap_2d_title_ <- "UMAP 2D plot for REFPOP genotype data"
      output_path_2d_umap <- paste0(
        output_geno_graphics_path,
        umap_refpop_train_data, "/",
        umap_refpop_train_data, "_genotype_refpop_umap_2d_",
        use_origin_family_or_genotype_as_label_,
        "_as_label.html"
      )
    }

    # modify label_title_ for family special case
    if (identical(use_origin_family_or_genotype_as_label_, "family") &&
      !identical(umap_refpop_train_data, "progeny")) {
      label_title_ <- paste0(
        "<b> ",
        str_to_title(use_origin_family_or_genotype_as_label_),
        " (except accession)",
        "</b>"
      )
    } else {
      label_title_ <- paste0(
        "<b> ",
        str_to_title(use_origin_family_or_genotype_as_label_),
        "</b>"
      )
    }

    fig_x_y <- plot_ly(
      type = "scatter", mode = "markers"
    ) %>%
      layout(
        plot_bgcolor = "#e5ecf6",
        title = umap_2d_title_,
        xaxis = list(title = "first component"),
        yaxis = list(title = "second component")
      )
    # regroup by label
    for (label_ in unique(geno_umap_2d$label)) {
      data_subset <- geno_umap_2d[geno_umap_2d$label == label_, ]
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
    # save graphics
    saveWidget(fig_x_y, file = output_path_2d_umap)

    # 3D plot
    # create base graphic
    if (identical(umap_refpop_train_data, "accessions") && predict_umap_progeny_) {
      umap_3d_title_ <- "UMAP 3D plot for REFPOP genotype data with umap trained
    on accessions, and progenies projected using trained model"
      output_path_3d_umap <- paste0(
        output_geno_graphics_path,
        umap_refpop_train_data, "/",
        umap_refpop_train_data, "_genotype_refpop_progeny_projected_umap_3d_",
        use_origin_family_or_genotype_as_label_,
        "_as_label.html"
      )
    } else {
      umap_3d_title_ <- "UMAP 3D plot for REFPOP genotype data"
      output_path_3d_umap <- paste0(
        output_geno_graphics_path,
        umap_refpop_train_data, "/",
        umap_refpop_train_data, "_genotype_refpop_umap_3d_",
        use_origin_family_or_genotype_as_label_,
        "_as_label.html"
      )
    }
    fig_x_y_z <- plot_ly(
      type = "scatter3d",
      mode = "markers"
    ) %>%
      layout(
        plot_bgcolor = "#e5ecf6",
        title = umap_3d_title_,
        xaxis = list(title = "first component"),
        yaxis = list(title = "second component"),
        zaxis = list(title = "third component")
      )
    # regroup by label
    for (label_ in unique(geno_umap_3d$label)) {
      data_subset <- geno_umap_3d[geno_umap_3d$label == label_, ]
      fig_x_y_z <- fig_x_y_z %>%
        add_trace(
          data = data_subset,
          x = ~X1, y = ~X2, z = ~X3,
          type = "scatter3d",
          mode = "markers",
          marker = list(color = color_labels_[label_]),
          name = label_
        )
    }
    fig_x_y_z <- fig_x_y_z %>% layout(
      legend = list(title = list(text = label_title_))
    )
    # save graphics
    saveWidget(fig_x_y_z, file = output_path_3d_umap)
  }

  # perform pca for geno_df
  geno_fam_orig_df_ <- geno_df[, c("Family", "Origin")]
  geno_df_ <- geno_df[, -match(c("Family", "Origin"), colnames(geno_df))]
  geno_df_ <- remove_monomorphic_markers(geno_df_)$filtered_df
  geno_pca_obj_ <- pca(geno_df_[, -match("Genotype", colnames(geno_df_))],
    ncomp = 500, center = TRUE, scale = TRUE
  )
  geno_pca_mat_ <- as.data.frame(geno_pca_obj_$variates$X)
  geno_pca_exp_var_ <- geno_pca_obj_$prop_expl_var$X
  geno_pca_cum_exp_var_ <- geno_pca_obj_$cum.var
  plot(geno_pca_cum_exp_var_)

  # plot coordinates of individuals on two first pcs :
  if (identical(use_origin_family_or_genotype_as_label_, "family")) {
    # define colors for labels
    labels_ <- unique(geno_fam_orig_df_$Family)
    n_family <- length(labels_)
    color_labels_ <- color_palette_family[1:n_family]
    names(color_labels_) <- labels_
    geno_pca_mat_$label <- geno_fam_orig_df_$Family

    # create plot
    fig_x_y <- plot_ly(
      type = "scatter", mode = "markers"
    ) %>%
      layout(
        plot_bgcolor = "#e5ecf6",
        title = "PCA 2D plot for REFPOP genotype data",
        xaxis = list(title = paste0(
          names(geno_pca_exp_var_)[1], ": ",
          signif(100 * as.numeric(geno_pca_exp_var_)[1], 2), "%"
        )),
        yaxis = list(title = paste0(
          names(geno_pca_exp_var_)[2], ": ",
          signif(100 * as.numeric(geno_pca_exp_var_)[2], 2), "%"
        ))
      )
    # regroup by label
    for (label_ in unique(geno_pca_mat_$label)) {
      data_subset <- geno_pca_mat_[geno_pca_mat_$label == label_, ]
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
      output_geno_graphics_path, pca_refpop_train_data, "/",
      pca_refpop_train_data, "_genotype_pca_family_as_label.html"
    ))
  } else if (identical(use_origin_family_or_genotype_as_label_, "origin")) {
    # define colors for labels
    labels_ <- unique(geno_fam_orig_df_$Origin)
    n_origin <- length(labels_)
    color_labels_ <- color_palette_origin[1:n_origin]
    names(color_labels_) <- labels_
    geno_pca_mat_$label <- geno_fam_orig_df_$Origin

    # create plot
    fig_x_y <- plot_ly(
      type = "scatter", mode = "markers"
    ) %>%
      layout(
        plot_bgcolor = "#e5ecf6",
        title = "PCA 2D plot for REFPOP genotype data",
        xaxis = list(title = paste0(
          names(geno_pca_exp_var_)[1], ": ",
          signif(100 * as.numeric(geno_pca_exp_var_)[1], 2), "%"
        )),
        yaxis = list(title = paste0(
          names(geno_pca_exp_var_)[2], ": ",
          signif(100 * as.numeric(geno_pca_exp_var_)[2], 2), "%"
        ))
      )
    # regroup by label
    for (label_ in unique(geno_pca_mat_$label)) {
      data_subset <- geno_pca_mat_[geno_pca_mat_$label == label_, ]
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
      output_geno_graphics_path, pca_refpop_train_data, "/",
      pca_refpop_train_data, "_genotype_pca_origin_as_label.html"
    ))
  }
}

# # get phased genotype data and split their columns according to each phase
# phased_geno_df <- readRDS(paste0(geno_dir_path, "phased_data/phased_genotypes.RDS"))
# genotype_names <- colnames(phased_geno_df)[-match(
#   "chromosome",
#   colnames(phased_geno_df)
# )]
# chromosome_num_col_ <- phased_geno_df$chromosome
#
# fam_df$sample.ID[fam_df$sample.ID %in% genotype_names]
#
# # initialize the cluster
# cl <- makeCluster(num_cores)
#
# # register the cluster
# registerDoParallel(cl)
#
# # apply the function to each column of phased_geno_df in parallel
# phased_geno_split_list <- foreach(
#   col = phased_geno_df[, genotype_names],
#   .combine = cbind
# ) %dopar% {
#   split_column(col)
# }
#
# # stop the cluster
# stopCluster(cl)

# convert the list to a dataframe
# phased_geno_split_df <- do.call(cbind, phased_geno_split_list)
# phased_geno_split_df <- as.data.frame(phased_geno_split_df)

# rename the columns of the result
# colnames(phased_geno_split_df) <- rep(genotype_names, each = 2)

# add chromosome number column
# phased_geno_split_df <- cbind(chromosome_num_col_, phased_geno_split_df)
# colnames(phased_geno_split_df)[1] <- "chromosome"
