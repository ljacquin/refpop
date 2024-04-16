# script meant to analyse refpop genotypic and phenotypic data
# note: text is formatted from Addins using Style active file from styler package

# clear memory, set options to increase memory and suppress warnings
rm(list = ls())
options(expressions = 5e5)
options(warn = -1)
emm_options(rg.limit = 20000)

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
library(lsmeans)
setwd(dirname(getActiveDocumentContext()$path))
source("../../functions.R")

# set paths and booleans
geno_data_path <- "../../data/genotype_data/"
pheno_data_path <- "../../data/phenotype_data/"
progeny_data_path <- "../../data/progeny_data/"
output_geno_graphics_path <- "../../data/graphics/geno_graphics/"
output_pheno_graphics_path <- "../../data/graphics/pheno_graphics/"

bed_file <- paste0(geno_data_path, "refpop_genotype.bed")
bim_file <- paste0(geno_data_path, "refpop_genotype.bim")
fam_file <- paste0(geno_data_path, "refpop_genotype.fam")

read_with_bigsnpr <- TRUE
read_with_snpStats <- FALSE # possible issue with the package or wrong usage
use_plot_tsne_ <- FALSE # computationally too demanding
use_ggplot_umap_ <- FALSE
use_plotly_umap_ <- TRUE # I prefer plotly, its a matter of personal taste

# threshold for removing columns with too much na
threshold_na <- 0.3

# define umap training and plot paraemeters

# define refpop train data for umap : complete, accessions, progeny
umap_refpop_train_data <- "progeny"

# define label data for umap :  origin, family or genotype (genotype not recommended)
use_origin_family_or_genotype_as_label_ <- "family"

# predict umap for progeny
predict_umap_progeny_ <- FALSE

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

# umap plots

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

    # save umap models for this case if predict_umap_progeny_ is true
  } else if (identical(umap_refpop_train_data, "accessions")) {
    
    sub_geno_df <- geno_df[-which(geno_df$Origin %in% "P"), ]

    # compute umap in 2D
    geno_umap_2d_model <- umap(sub_geno_df[, -idx_origin_geno_names],
      n_components = 2, random_state = 15
    )
    geno_umap_2d <- data.frame(geno_umap_2d_model[["layout"]])

    # compute umap in 3D
    geno_umap_3d_model <- umap(sub_geno_df[, -idx_origin_geno_names],
      n_components = 3, random_state = 15
    )
    geno_umap_3d <- data.frame(geno_umap_3d_model[["layout"]])
    
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

  # treat umap prediction for progenies based on unsupervised learning
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
      umap_refpop_train_data, "_refpop_progeny_projected_umap_2d_",
      use_origin_family_or_genotype_as_label_,
      "_as_label.html"
    )
  }else{
    umap_2d_title_ <- "UMAP 2D plot for REFPOP genotype data"
    output_path_2d_umap <- paste0(
      output_geno_graphics_path,
      umap_refpop_train_data, "/",
      umap_refpop_train_data, "_refpop_umap_2d_",
      use_origin_family_or_genotype_as_label_,
      "_as_label.html"
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
  fig_x_y <- fig_x_y %>%layout(  
    legend=list(title=list(text=paste0('<b> ',
   str_to_title(use_origin_family_or_genotype_as_label_), '</b>'))))
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
      umap_refpop_train_data, "_refpop_progeny_projected_umap_3d_",
      use_origin_family_or_genotype_as_label_,
      "_as_label.html"
    )
  }else{
    umap_3d_title_ <- "UMAP 3D plot for REFPOP genotype data"
    output_path_3d_umap <- paste0(
      output_geno_graphics_path,
      umap_refpop_train_data, "/",
      umap_refpop_train_data, "_refpop_umap_3d_",
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
  fig_x_y_z <- fig_x_y_z %>%layout(  
    legend=list(title=list(text=paste0('<b> ',
   str_to_title(use_origin_family_or_genotype_as_label_), '</b>'))))
  # save graphics
  saveWidget(fig_x_y_z, file = output_path_3d_umap)
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

# phased_geno <- readRDS(paste0(geno_data_path,'phased_geno.RDS'))

# ---------- phenotype data analysis

path_spats_adj_pheno <- paste0(pheno_data_path,'spats_adjusted_phenotypes/')
files_names_spats_adj_pheno <- list.files(path_spats_adj_pheno)
trait_names_ <- str_replace_all(files_names_spats_adj_pheno,
                                '_spats_adjusted_.*', 
                                replacement = '')

multi_env_clonal_h2_traits_ <- rep(0, length(trait_names_))
names(multi_env_clonal_h2_traits_) <- trait_names_

list_ls_means_adj_pheno_per_geno <- vector('list', length(trait_names_))
names(list_ls_means_adj_pheno_per_geno) <- trait_names_


for ( file_ in files_names_spats_adj_pheno){
  print(paste0("computation for file : ", file_))
  
  df_ <- as.data.frame(fread(paste0(path_spats_adj_pheno, file_)))
  df_ <- df_[, c('Genotype','Envir',colnames(df_)[str_detect(colnames(df_), 
                                           'spats_adj_pheno_')])]
  Y <- colnames(df_)[str_detect(colnames(df_),'spats_adj_pheno_')]

  if (length(unique(df_$Envir))>1){
    # compute multi-location clonal mean heritability
    lmer_model_ <- lmer(as.formula(paste0(Y, 
                  " ~ 1 + (1 | Genotype) + Envir + (1 | Genotype:Envir)")),
                  data = df_)
    nl <- length(unique(df_$Envir))
    
    df_group_by <- df_ %>%
      group_by(Genotype, Envir) %>%
      summarise(frequency = n()) %>% 
      data.frame()
    
    nr_bar_ <- mean(df_group_by$frequency)
    
    multi_env_clonal_h2_traits_[str_replace_all(file_, '_spats_adjusted_.*', 
                                                replacement = '')] <- 
      compute_multi_location_clonal_mean_h2(lmer_model_, nr_bar_, nl)
    
    # compute adjusted ls-means for genotypes across environments 
    lm_model <- lm(formula(paste0(Y, '~ Genotype + Envir')), data = df_)
    ls_means <- as.data.frame(lsmeans(lm_model, ~ Genotype, 
                ref.grid = list(rg.limit = 20000)))[ ,c('Genotype', 'lsmean')]
    colnames(ls_means)[match('lsmean', colnames(ls_means))] <- paste0(
    str_replace_all(file_, '_spats_adjusted_.*', replacement = ''), '_lsmean')
    
    list_ls_means_adj_pheno_per_geno[[str_replace_all(file_, '_spats_adjusted_.*', 
                                  replacement = '')]]<-ls_means
  
  }else{
    
    # compute multi-location clonal mean heritability for unique environment
    lmer_model_ <- lmer(as.formula(paste0(Y, 
                    " ~ 1 + (1 | Genotype)")), data = df_)
    nr_bar_ <- mean(table(df_$Genotype))
    
    multi_env_clonal_h2_traits_[str_replace_all(file_, '_spats_adjusted_.*', 
              replacement = '')]<- compute_indiv_location_clonal_mean_h2(
                lmer_model_, nr_bar_ )
    
    # compute adjusted ls-means for genotypes for unique environment
    lm_model <- lm(formula(paste0(Y, '~ Genotype')), data = df_)
    ls_means <- as.data.frame(lsmeans(lm_model, ~ Genotype, 
                ref.grid = list(rg.limit = 20000)))[ ,c('Genotype', 'lsmean')]
    
    colnames(ls_means)[match('lsmean', colnames(ls_means))] <- paste0(
      str_replace_all(file_, '_spats_adjusted_.*', replacement = ''), '_lsmean')
    
    list_ls_means_adj_pheno_per_geno[[str_replace_all(file_, '_spats_adjusted_.*', 
                                              replacement = '')]]<-ls_means
  }
}

#  merge list of ls_means into a single data frame for genotypes
merged_df <- Reduce(function(x, y) merge(x, y, by = "Genotype", all = TRUE),
                    list_ls_means_adj_pheno_per_geno)

# convert merge object to data.frame
merged_df <- as.data.frame(merged_df)
na_count <- colSums(is.na(merged_df))
filtered_df <- na.omit(merged_df[, na_count/nrow(merged_df) <= threshold])
colnames(filtered_df)[str_detect(colnames(filtered_df), "_lsmean")] <-
  str_replace_all(colnames(filtered_df)[str_detect(colnames(filtered_df), "_lsmean")],
                  pattern = "_lsmean", replacement = '')

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
saveWidget(heatmap_pearson, file = paste0(output_pheno_graphics_path,
                     'pearson_cor_adj_pheno_ls_means_per_geno_all_env.html'))

# create a bar chart for the multi-location clonal mean h2
traits <- names(multi_env_clonal_h2_traits_)
values <- as.numeric(multi_env_clonal_h2_traits_)

# create the bar chart
bar_plot <- plot_ly(x = traits, 
                y = values, 
                type = 'bar', 
                marker = list(color = 'rgb(158,202,225)',
                              line = list(color = 'rgb(8,48,107)', 
                                          width = 1.5))) %>%
  layout(title = "Multi-location clonal mean heritability (h2) computed from adjusted phenotypes",
         xaxis = list(title = "Traits"),
         yaxis = list(title = "Heritability (h2)"))

saveWidget(bar_plot, file = paste0(output_pheno_graphics_path,
              'multi_location_clonal_h2_bar_plot_multi_location_clonal_mean_h2.html'))

