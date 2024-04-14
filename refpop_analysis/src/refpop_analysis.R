# script meant to analyse refpop genotypic and phenotypic data
# note: text is formatted from Addins using Style active file from styler package

# clear memory, set options to increase memory and suppress warnings
rm(list = ls())
options(expressions = 5e5)
options(warn = -1)

# source libraries
library(bigsnpr)
library(reticulate)
install_other_requirements <- FALSE
if (install_other_requirements) {
  install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install("snpStats")
  BiocManager::install("M3C")
  install.packages("remotes")
  remotes::install_github("hemstrow/snpR")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
}
library(snpR)
library(M3C)
library(Rtsne)
library(snpStats)
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(dplyr)
library(htmlwidgets)
library(rstudioapi)
library(stringr)
setwd(dirname(getActiveDocumentContext()$path))
source("../../functions.R")

# set paths and booleans
geno_data_path <- "../../data/genotype_data/"
pheno_data_path <- "../../data/phenotype_data/"
progeny_data_path <- "../../data/progeny_data/"
output_geno_graphics_path <- "../../data/graphics/geno_graphics/"
read_with_bigsnpr <- TRUE
read_with_snpStats <- FALSE # possible issue with the package or wrong usage
use_plot_tsne_ <- FALSE # computationally too demanding
use_ggplot_umap_ <- FALSE
use_plotly_umap_ <- TRUE # I prefer plotly, its a matter of personal taste

# set path for .bim, .bed and .fam files
bed_file <- paste0(geno_data_path, "refpop_genotype.bed")
bim_file <- paste0(geno_data_path, "refpop_genotype.bim")
fam_file <- paste0(geno_data_path, "refpop_genotype.fam")

# ---------- genotype data analysis

# read plink genotype data
if (read_with_bigsnpr) {
  # reading the bedfile for the first time or its associated rds file if exists
  # note .bed and .fam files should be in the same directory
  if (!file.exists(paste0(geno_data_path, "refpop_genotype.rds"))) {
    df_ <- snp_readBed(bed_file)
    df_ <- readRDS(paste0(geno_data_path, "refpop_genotype.rds"))
  } else {
    df_ <- readRDS(paste0(geno_data_path, "refpop_genotype.rds"))
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
  geno_df <- geno_df %>% select(Genotype, everything())
}

# possible other way to read plink genotype data but
# the usage might be wrong or may be read.plink() has an issue
if (read_with_snpStats) {
  # read plik files
  df_ <- read.plink(bed_file, bim_file, fam_file)
  geno_df <- df_[["genotypes"]]
  fam_df <- df_[["fam"]]
  map_df <- df_[["map"]]
}

# get geographical origins of genotypes
geno_origin_ <- as.data.frame(fread(paste0(
  geno_data_path,
  "geographical_origins.csv"
)))

# rename to a common key, i.e. genotype
colnames(geno_origin_)[
  colnames(geno_origin_) == "Genotype Code"
] <- "Genotype"

# get families of genotypes
geno_fam_ <- as.data.frame(fread(paste0(
  geno_data_path,
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

# merge geno_df with geno_origin
geno_df <- merge(geno_df, geno_origin_[, c("Genotype", "Family", "Origin")],
  by = "Genotype", all = TRUE
)
geno_df <- geno_df %>% select(c(Genotype, Family, Origin), everything())
idx_origin_geno_names <- which(colnames(geno_df) %in% c(
  "Genotype",
  "Family",
  "Origin"
))

# t-SNE plots, sample columns if necessary to prevent stack overflow
if (use_plot_tsne_) {
  set.seed(42)
  sample_size_ <- 15000

  sel_cols_ <- sort(sample(
    x = 1:ncol(geno_df),
    size = sample_size_,
    replace = FALSE
  ))

  sel_cols_ <- c(
    -idx_origin_geno_names,
    sel_cols_
  )

  # Rtsne with snpR
  Rtsne_out <- Rtsne(
    geno_df[, sel_cols_],
    perplexity = 30
  )
  plot(Rtsne_out$Y)

  # tsne with M3C
  tsne_out <- tsne(
    unique(geno_df[, sel_cols_]),
    perplex = 30
  )
  plot(tsne_out$Y)
}

# umap plots, much less computationally demanding than tsne and preferable

# define refpop train data for umap : complete, accessions, progeny
umap_refpop_train_data <- "complete"
# define label data for umap :  origin, family or genotype (genotype not recommended)
use_origin_family_or_genotype_as_label_ <- "origin"

# set palette of colors according to label used
set.seed(123)

# set arbitrary colors for unique genotypes
if (identical(use_origin_family_or_genotype_as_label_, "genotype")) {
  colfunc <- colorRampPalette(c(
    "black", "blue", "red", "orange",
    "yellow", "green"
  ))
  color_palette <- colfunc(length(unique(geno_df$Genotype)))

  # set color labels for families (28 counts)
} else if (identical(use_origin_family_or_genotype_as_label_, "family")) {
  color_palette <- c(
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
  # if population is not complete select the corresponding color for subset
  if (identical(umap_refpop_train_data, "accessions")) {
    color_palette <- color_palette[1]
  } else if (identical(umap_refpop_train_data, "progeny")) {
    color_palette <- color_palette[-1]
  }

  # set color labels for origins (11 counts)
} else {
  color_palette <- c(
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
}

if (use_plotly_umap_) {
  if (identical(umap_refpop_train_data, "complete")) {
    # compute umap in 2D
    geno_umap_2d <- data.frame(umap(geno_df[, -idx_origin_geno_names],
      n_components = 2, random_state = 15
    )[["layout"]])
    # compute umap in 3D
    geno_umap_3d <- data.frame(umap(geno_df[, -idx_origin_geno_names],
      n_components = 3, random_state = 15
    )[["layout"]])
  } else if (identical(umap_refpop_train_data, "accessions")) {
    sub_geno_df <- geno_df[-which(geno_df$Origin %in% "P"), ]

    # compute umap in 2D
    geno_umap_2d <- data.frame(umap(sub_geno_df[, -idx_origin_geno_names],
      n_components = 2, random_state = 15
    )[["layout"]])
    # compute umap in 3D
    geno_umap_3d <- data.frame(umap(sub_geno_df[, -idx_origin_geno_names],
      n_components = 3, random_state = 15
    )[["layout"]])
  } else if (identical(umap_refpop_train_data, "progeny")) {
    sub_geno_df <- geno_df[which(geno_df$Origin %in% "P"), ]

    # compute umap in 2D
    geno_umap_2d <- data.frame(umap(sub_geno_df[, -idx_origin_geno_names],
      n_components = 2, random_state = 15
    )[["layout"]])
    # compute umap in 3D
    geno_umap_3d <- data.frame(umap(sub_geno_df[, -idx_origin_geno_names],
      n_components = 3, random_state = 15
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

  # 2D plot
  # create base graphic
  fig_x_y <- plot_ly(
    type = "scatter", mode = "markers"
  ) %>%
    layout(
      plot_bgcolor = "#e5ecf6",
      title = "UMAP 2D plot for REFPOP genotype data",
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
  # save graphics
  saveWidget(fig_x_y, file = paste0(
    output_geno_graphics_path,
    umap_refpop_train_data, "/",
    umap_refpop_train_data, "_refpop_umap_2d_",
    use_origin_family_or_genotype_as_label_,
    "_as_label.html"
  ))

  # 3D plot
  # create base graphic
  fig_x_y_z <- plot_ly(
    type = "scatter3d",
    mode = "markers"
  ) %>%
    layout(
      plot_bgcolor = "#e5ecf6",
      title = "UMAP 3D plot for REFPOP genotype data",
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
  # save graphics
  saveWidget(fig_x_y_z, file = paste0(
    output_geno_graphics_path,
    umap_refpop_train_data, "/",
    umap_refpop_train_data, "_refpop_umap_3d_",
    use_origin_family_or_genotype_as_label_,
    "_as_label.html"
  ))
}

if (use_ggplot_umap_) {
  geno_umap_2d$origin <- geno_df$Origin
  # umap 2D plot
  ggplot(geno_umap_2d, aes(x = X1, y = X2, color = origin)) +
    geom_point() +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    labs(x = "", y = "", title = "UMAP 2D plot for REFPOP genotype data") +
    scale_fill_manual(
      values = c(
        "WCE" = "blue", "SE" = "red", "USA" = "green",
        "ANZ" = "orange", "JPN" = "purple", "P" = "pink",
        "SEE" = "yellow", "NEE" = "brown", "CAN" = "cyan",
        "U" = "gray", "ZAF" = "magenta"
      ),
      labels = c(
        "WCE" = "WCE", "SE" = "SE", "USA" = "USA",
        "ANZ" = "ANZ", "JPN" = "JPN", "P" = "P",
        "SEE" = "SEE", "NEE" = "NEE", "CAN" = "CAN",
        "U" = "U", "ZAF" = "ZAF"
      )
    ) +
    theme(
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      plot.title = element_text(hjust = 0.5)
    )
}


# ---------- phenotype data analysis
