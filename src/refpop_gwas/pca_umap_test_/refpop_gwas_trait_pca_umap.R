# script meant to perform wiser, adjusted ls-means, and population structure
# corrected adjusted ls-means mlmm gwas for the refpop
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
if ("refpop_env" %in% conda_list()$name) {
  print("using refpop_env")
  use_condaenv("refpop_env")
}
# install other requirements from github if necessary
install_other_requirements <- F
if (install_other_requirements) {
  # reticulate::install_miniconda()
  conda_create("refpop_env")
  use_condaenv("refpop_env")
  library(devtools)
  devtools::install_github("ljacquin/KRMM")
  devtools::install_github("rstudio/tensorflow")
  library(tensorflow)
  install_tensorflow(envname = "refpop_env")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
  install.packages("umap")
}
use_tensorflow_ <- F
if (use_tensorflow_) {
  # leave tensorflow and keras for later use
  library(tensorflow)
  library(keras3)
  tensorflow::tf$random$set_seed(0)
  py_module_available("keras") # must return TRUE
  py_module_available("tensorflow") # must return TRUE
  py_discover_config("keras") # more info on the python env, tf and keras
}
library(umap)
library(MASS)
library(data.table)
library(stringr)
library(lme4)
library(tidyr)
library(FactoMineR)
library(doParallel)
library(doRNG)
library(robustbase)
library(foreach)
library(parallel)
library(missForest)
library(Matrix)
library(matrixcalc)
library(rgl)
library(Rfast)
library(cvTools)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(KRMM)
library(kernlab)
library(whitening)
library(glmnet)
library(ranger)
library(mixOmics)
library(future)
library(future.apply)
# library(GAPIT)
library(mlmm.gwas)
library(fastqq)
# get all mlmm.gwas from namespace and make them globally visible
list_mlmm_gwas_fcts <- ls(getNamespace("mlmm.gwas"), all.names = TRUE)
for (f_ in list_mlmm_gwas_fcts) {
  assign(f_, get(f_, envir = getNamespace("mlmm.gwas")),
    envir = .GlobalEnv
  )
}

# define computation mode, i.e. local or cluster
computation_mode <- "local"

# if comutations are local in rstudio, detect and set script path
# automatically using rstudioapi
if (identical(computation_mode, "local")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}

# source functions
source("../functions.R")

# set options
options(future.globals.maxSize = 60 * 1024^3)
options(expressions = 5e5)
options(warn = -1)

# define mlmm maximum number of steps for forward stepwise regression
max_step_mlmm_fwd_reg_ <- 3

# define alpha for corrected Bonferroni threshold
alpha_bonf_ <- 0.05

# define proportion of explained variance required to defined number of pcs
prop_var_ <- 0.95

# define number of snp to sample for pca
snp_sample_size_ <- 50e3

# define number of retained components for mixOmics (NB. all these components
# will not be used)
n_comp_pca_ <- 500

# define number of retained components for umap
n_comp_umap_ <- 100

# define umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 15
min_dist_ <- 0.1

# should a design matrix be computed for the groups accession and progeny, and
# families ?
compute_fix_eff_design_fam_ <- F

# should pca and/or umap be computed for population structure consideration in
# mlmm ?
compute_pca_for_mlmm_ <- F
compute_umap_for_mlmm_ <- F

# define prefix of markers
marker_prefix_ <- "AX-"

# set input paths
geno_dir_path <- "../../data/genotype_data/"
pheno_dir_path <- "../../data/phenotype_data/"

# set path for wiser phenotypes estimated using whitening
wiser_pheno_dir_path <- "../../data/phenotype_data/wiser_phenotype_estimates/"

# output result path for genotype graphics
# output_pred_results_path <- "../../results/gwas/"
output_pred_graphics_path <- "../../results/graphics/gwas_graphics/"

# define list of all traits
list_all_traits_ <- c(
  "Harvest_date", "Fruit_weight", "Fruit_number",
  "Fruit_weight_single", "Color_over", "Russet_freq_all",
  "Trunk_diameter", "Trunk_increment", "Flowering_intensity",
  "Flowering_begin", "Flowering_full", "Flowering_end",
  "Scab", "Powdery_mildew", "Weight_sample"
)

# should wiser and genomic prediction be performed for a subset of traits ?
# this can be useful for managing computation resources
use_subset_list <- F

# define traits_ for wiser
if (use_subset_list) {
  traits_ <- c(
    "Color_over", "Weight_sample"
  )
} else {
  traits_ <- list_all_traits_
}

# # get trait arguments
# args <- commandArgs(trailingOnly = T)
# trait_num <- as.integer(args[1])

trait_num <- 1

# define trait_
trait_ <- traits_[trait_num]
print(paste0("trait: ", trait_))

# create result repository for trait_ if it does not exist
if (!dir.exists(paste0(output_pred_graphics_path, trait_))) {
  dir.create(paste0(output_pred_graphics_path, trait_))
}

# get wiser and adjusted ls-means phenotypes, and genotype data (genomic data)

# get trait wiser phenotypes
wiser_trait_ <- readRDS(paste0(
  wiser_pheno_dir_path,
  paste0("wiser_obj_linear_kernel_", trait_)
))$wiser_phenotypes
colnames(wiser_trait_)[2] <- trait_

# get trait ls-means phenotypes
ls_means_trait_ <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "adjusted_ls_mean_phenotypes.csv"
)))[, c("Genotype", trait_)]

# get marker data
marker_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_data.csv"
)))

# get map data
map_df <- as.data.frame(readRDS(paste0(
  geno_dir_path,
  "phased_data/chr_map.RDS"
)))

# remove monomorphic markers
marker_df <- remove_monomorphic_markers(marker_df)
monomorphic_markers_list_ <- marker_df$monomorphic_markers
marker_df <- marker_df$filtered_df

# merge ls_means_trait_ and marker_df for integrity of analyses
merged_df <- merge(wiser_trait_, marker_df, by = "Genotype")
wiser_trait_ <- merged_df[, c("Genotype", trait_)]

# merge ls_means_trait_ and marker_df for integrity of analyses
merged_df <- merge(ls_means_trait_, marker_df, by = "Genotype")
ls_means_trait_ <- merged_df[, c("Genotype", trait_)]

# rename rows and remove genotype and trait_ columns
geno_names <- merged_df$Genotype
marker_df <- merged_df[, -match(c("Genotype", trait_), colnames(merged_df))]
rownames(marker_df) <- geno_names

if (compute_fix_eff_design_fam_) {
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
  colnames(geno_fam_)[
    match("Family", colnames(geno_fam_))
  ] <- "Family_"

  # detect and set families for genotypes
  geno_origin_$Family_ <- NA
  for (i in 1:length(geno_fam_$Pattern)) {
    idx_geno_fam <- which(str_detect(
      geno_origin_$Genotype,
      pattern = geno_fam_$Pattern[i]
    ))
    geno_origin_$Family_[idx_geno_fam] <- geno_fam_$Family_[i]
  }
  geno_origin_$Group_ <- geno_origin_$Family_
  geno_origin_$Group_[!is.na(geno_origin_$Group_)] <- "Progeny"
  geno_origin_$Group_[is.na(geno_origin_$Group_)] <- "Accession"
  geno_origin_$Family_[is.na(geno_origin_$Family_)] <- "Dummy_family"

  # merge geno_origin_ with ls_means_trait_ data frame for integrity of analyses,
  # w.r.t to population structure correction, and create design matrix for the
  # groups and families fixed effects (i.e. population structure correction)
  merged_df <- merge(ls_means_trait_, geno_origin_, by = "Genotype")
  fix_eff_design_mat_ <- model.matrix(
    as.formula(paste0(
      "~ ",
      paste0(c("Group_", "Family_"),
        collapse = "+"
      ), "- 1"
    )),
    data = merged_df
  )
  fix_eff_design_mat_ <- fix_eff_design_mat_[
    , !str_detect(colnames(fix_eff_design_mat_), "Dummy_family")
  ]
  rownames(fix_eff_design_mat_) <- merged_df$Genotype
  colnames(fix_eff_design_mat_) <- str_replace_all(
    colnames(fix_eff_design_mat_),
    pattern = " ", replacement = "_"
  )
}
rm(merged_df)

# get common markers between map and marker data, and order them
marker_df_snp_names_ <- data.frame(
  "index" = 1:ncol(marker_df),
  "snp.name" = colnames(marker_df)
)
common_snp_names_ <- merge(
  map_df, marker_df_snp_names_,
  by = "snp.name"
)
common_snp_names_ <- common_snp_names_[order(
  common_snp_names_$chromosome,
  common_snp_names_$position
), ]

# reorder the markers (if necessary), and verify equality for integrity
marker_df <- marker_df[, common_snp_names_$snp.name]
map_df <- map_df[match(common_snp_names_$snp.name, map_df$snp.name), ]
map_df <- map_df[, c("snp.name", "chromosome", "position")]
identical(map_df$snp.name, colnames(marker_df))

# define corrected Bonferroni threshold
threshold_ <- alpha_bonf_ / ncol(marker_df)

# compute gram matrix
marker_matrix <- scale(apply(marker_df, 2, as.numeric), scale = F)
k_mat <- marker_matrix %*% t(marker_matrix)
rownames(marker_matrix) <- rownames(marker_df)
colnames(marker_matrix) <- colnames(marker_df)
colnames(k_mat) <- rownames(k_mat) <- rownames(marker_df)

if (compute_pca_for_mlmm_ | compute_umap_for_mlmm_) {
  # sample snp to accelerate pca and umap computation (NB. >= 50e3 is enough to
  # capture population structure for the refpop, see Jung et al. 2020 or
  # Jacquin et al. 2025)
  set.seed(123)
  idx_snp_ <- sort(
    sample(1:ncol(marker_matrix), size = snp_sample_size_, replace = F),
    decreasing = F
  )
}

if (compute_pca_for_mlmm_) {
  # compute pca for marker data
  geno_pca_obj_ <- mixOmics::pca(
    marker_matrix[, idx_snp_],
    ncomp = n_comp_pca_
  )
  geno_pca_mat_ <- as.data.frame(geno_pca_obj_$variates$X)
  plot(geno_pca_obj_$cum.var)
  idx_explained_var_ <- as.numeric(which(geno_pca_obj_$cum.var >= prop_var_)[1])
  if (is.na(idx_explained_var_)) {
    idx_explained_var_ <- n_comp_pca_
  }
}

if (compute_umap_for_mlmm_) {
  # compute umap for marker data
  geno_umap_mat_ <- data.frame(umap(marker_matrix[, idx_snp_],
    n_components = n_comp_umap_, random_state = random_state_umap_,
    n_neighbors = n_neighbors_umap_, min_dist = min_dist_
  )[["layout"]])
}

# mlmm gwas for wiser

# reformat response for wiser
wiser_y <- wiser_trait_[, trait_]
names(wiser_y) <- wiser_trait_$Genotype

# fit mlmm forward stepwise regression algorithm for wiser phenotypes
mlmm_components_wiser_ <- mlmm_allmodels_(
  wiser_y, list(marker_df), list(k_mat),
  maxsteps = max_step_mlmm_fwd_reg_
)

# get mlmm algorithm fitted components (i.e. p-values and degrees of freedom)
mlmm_p_val_wiser_ <- mlmm_components_wiser_$p_val
mlmm_deg_freedom_wiser_ <- mlmm_components_wiser_$df_list

# get step associated to "optimal" model based on ebic minimization
res_eBIC <- as.data.frame(eBIC_allmodels(
  wiser_y,
  frommlmm_toebic(list(marker_matrix), mlmm_p_val_wiser_),
  list(k_mat), ncol(marker_matrix)
))
optimal_model_associated_step_ <- which.min(res_eBIC$eBIC_0.75[-1])

# make a manhattan plot
manhattan_plot_done_ <- manhattan_plot_(mlmm_p_val_wiser_,
  map = map_df,
  steps = optimal_model_associated_step_,
  output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
  file_name_ = paste0(trait_, "_wiser_manhattan_plot.png"),
  threshold_ = -log10(threshold_),
  main_ = paste0(
    "MLMM GWAS for ", trait_,
    " WISER phenotypes for REFPOP,
    using a corrected Bonferroni threshold for 0.05"
  )
)
manhattan_plot_done_

# get p-values associated with the step of the identified "optimal" model (i.e.,
# the model minimizing the extended-BIC). Note that these p-values do not correspond
# to the "optimal" model itself but to other models of the same form under H1 that
# were not significant at the step of the identified "optimal" model.
# The p-values specifically associated with the "optimal" model under H1 are those
# corresponding to the SNPs included in this model.
p_val_step_opt_model_ <- mlmm_p_val_wiser_[[
  optimal_model_associated_step_ + 1
]]
head(sort(p_val_step_opt_model_, decreasing = F))

# get degrees of freedom of "optimal" model (i.e. minimizing extended-BIC). In
# other words, get df1 = q = p1 (nb. param. under H1) - p0 (nb. param. under H0)
# which equals 1 in case of the nullity test of q = 1 variable (SNP) under H0,
# and df2 = n - p1 which decreases by 1 at each step (i.e. p1 increases by 1 at
# each step).
deg_freed_step_opt_model_ <- sort(unique(mlmm_deg_freedom_wiser_[[
  optimal_model_associated_step_ + 1
]]))

# make mlmm gwas qqplot for trait WISER estimated phenotypes
qq_plot_done_ <- qqplot_f_distrib_mlmm_gwas_(
  p_val_obs_ = p_val_step_opt_model_,
  df1 = as.numeric(deg_freed_step_opt_model_[1]),
  df2 = as.numeric(deg_freed_step_opt_model_[2]),
  seed_ = 42,
  main_ = paste0(
    "MLMM GWAS Q-Q plot for ",
    trait_, " WISER estimated phenotypes"
  ),
  output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
  file_name_ = paste0(trait_, "_wiser_qq_plot.png")
)
qq_plot_done_

# get significant snp according to defined corrected Bonferroni threshold
signif_snp_threshold_ <- as.data.frame(threshold_allmodels(
  threshold = threshold_, mlmm_p_val_wiser_
))
if (
  sum(grepl(marker_prefix_, signif_snp_threshold_$`p-value`)) > 0) {
  signif_snp_threshold_ <- correct_signif_snp_threshold(
    signif_snp_threshold_
  )
}

# make boxplots for trait phenotypes according to SNP allele combinations (
# i.e. homozygote and heterozygote ) of genotypes, and write significant markers
# infos
if (length(signif_snp_threshold_$SNP) > 0) {
  for (snp_ in signif_snp_threshold_$SNP) {
    # get marker p-value for boxplot
    p_val_ <- signif(
      as.numeric(
        signif_snp_threshold_$p.value[
          match(snp_, signif_snp_threshold_$SNP)
        ]
      ), 3
    )
    df_ <- data.frame(
      "wiser_y" = wiser_y, "marker" = marker_matrix[, snp_]
    )
    r2_ <- summary(
      lm(as.formula("wiser_y ~ 1 + marker"), data = df_)
    )$r.squared
    r2_ <- signif(as.numeric(r2_), 3)
    geno_box_plot_ <- genotypes_boxplot_(
      X = marker_matrix,
      Y = wiser_y,
      ylab_ = paste0(trait_, " wiser phenotypes"),
      xlab_ = paste0(
        "Genotype alleles combination on chromosome ",
        map_df$chromosome[match(snp_, map_df$snp.name)]
      ),
      markers = snp_,
      r2_ = r2_,
      p_val_ = p_val_,
      output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
      file_name_ = paste0(
        trait_, "_wiser_gwas_significant_marker_", snp_, ".png"
      )
    )
  }
  # get and write significant markers infos
  signif_marker_info <- data.frame(
    "SNP" = signif_snp_threshold_$SNP,
    "p_value" = as.numeric(signif(
      as.numeric(signif_snp_threshold_[, 2]), 3
    )),
    "chromosome_number" = map_df$chromosome[
      map_df$snp.name %in% signif_snp_threshold_$SNP
    ],
    "position_in_bp" = map_df$position[
      map_df$snp.name %in% signif_snp_threshold_$SNP
    ]
  )
  fwrite(signif_marker_info, file = paste0(
    paste0(output_pred_graphics_path, trait_, "/"),
    trait_, "_wiser_gwas_significant_marker_info.csv"
  ))
}

# mlmm gwas for adjusted ls-means

# reformat response for ls_means
ls_means_y <- ls_means_trait_[, trait_]
names(ls_means_y) <- ls_means_trait_$Genotype

# fit mlmm forward stepwise regression algorithm for adjusted ls-means phenotypes
mlmm_components_ls_means_ <- mlmm_allmodels_(
  ls_means_y, list(marker_df), list(k_mat),
  maxsteps = max_step_mlmm_fwd_reg_
)

# get mlmm algorithm fitted components (i.e. p-values and degrees of freedom)
mlmm_p_val_ls_means_ <- mlmm_components_ls_means_$p_val
mlmm_deg_freedom_ls_means_ <- mlmm_components_ls_means_$df_list

# get step associated to "optimal" model based on ebic minimization
res_eBIC <- as.data.frame(eBIC_allmodels(
  ls_means_y,
  frommlmm_toebic(list(marker_matrix), mlmm_p_val_ls_means_),
  list(k_mat), ncol(marker_matrix)
))
optimal_model_associated_step_ <- which.min(res_eBIC$eBIC_0.75[-1])

# make a manhattan plot
manhattan_plot_done_ <- manhattan_plot_(
  mlmm_p_val_ls_means_,
  map = map_df,
  steps = optimal_model_associated_step_,
  output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
  file_name_ = paste0(trait_, "_ls_means_manhattan_plot.png"),
  threshold_ = -log10(threshold_),
  main_ = paste0(
    "MLMM GWAS for ", trait_,
    " LS-means phenotypes for REFPOP,
    using a corrected Bonferroni threshold for 0.05"
  )
)
manhattan_plot_done_

# get p-values associated with the step of the identified "optimal" model (i.e.,
# the model minimizing the extended-BIC). Note that these p-values do not correspond
# to the "optimal" model itself but to other models of the same form under H1 that
# were not significant at the step of the identified "optimal" model.
# The p-values specifically associated with the "optimal" model under H1 are those
# corresponding to the SNPs included in this model.
p_val_step_opt_model_ <- mlmm_p_val_ls_means_[[
  optimal_model_associated_step_ + 1
]]
head(sort(p_val_step_opt_model_, decreasing = F))

# get degrees of freedom of "optimal" model (i.e. minimizing extended-BIC). In
# other words, get df1 = q = p1 (nb. param. under H1) - p0 (nb. param. under H0)
# which equals 1 in case of the nullity test of q = 1 variable (SNP) under H0,
# and df2 = n - p1 which decreases by 1 at each step (i.e. p1 increases by 1 at
# each step).
deg_freed_step_opt_model_ <- sort(unique(mlmm_deg_freedom_ls_means_[[
  optimal_model_associated_step_ + 1
]]))

# make mlmm gwas qqplot for trait adjsuted ls-means estimated phenotypes
qq_plot_done_ <- qqplot_f_distrib_mlmm_gwas_(
  p_val_obs_ = p_val_step_opt_model_,
  df1 = as.numeric(deg_freed_step_opt_model_[1]),
  df2 = as.numeric(deg_freed_step_opt_model_[2]),
  seed_ = 42,
  main_ = paste0(
    "MLMM GWAS Q-Q plot for ",
    trait_, " adjusted LS-means estimated phenotypes"
  ),
  output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
  file_name_ = paste0(trait_, "_ls_means_qq_plot.png")
)
qq_plot_done_

# get significant snp according to defined corrected Bonferroni threshold
signif_snp_threshold_ <- as.data.frame(threshold_allmodels(
  threshold = threshold_, mlmm_p_val_ls_means_
))
if (
  sum(grepl(marker_prefix_, signif_snp_threshold_$`p-value`)) > 0) {
  signif_snp_threshold_ <- correct_signif_snp_threshold(
    signif_snp_threshold_
  )
}

# make boxplots for trait phenotypes according to SNP allele combinations (
# i.e. homozygote and heterozygote ) of genotypes, and write significant markers
# infos
if (length(signif_snp_threshold_$SNP) > 0) {
  for (snp_ in signif_snp_threshold_$SNP) {
    # get marker p-value for boxplot
    p_val_ <- signif(
      as.numeric(
        signif_snp_threshold_$p.value[
          match(snp_, signif_snp_threshold_$SNP)
        ]
      ), 3
    )
    df_ <- data.frame(
      "ls_means_y" = ls_means_y, "marker" = marker_matrix[, snp_]
    )
    r2_ <- summary(
      lm(as.formula("ls_means_y ~ 1 + marker"), data = df_)
    )$r.squared
    r2_ <- signif(as.numeric(r2_), 3)
    geno_box_plot_ <- genotypes_boxplot_(
      X = marker_matrix,
      Y = ls_means_y,
      ylab_ = paste0(trait_, " adjusted LS-means phenotypes"),
      xlab_ = paste0(
        "Genotype alleles combination on chromosome ",
        map_df$chromosome[match(snp_, map_df$snp.name)]
      ),
      markers = snp_,
      r2_ = r2_,
      p_val_ = p_val_,
      output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
      file_name_ = paste0(
        trait_, "_ls_means_gwas_significant_marker_",
        snp_, ".png"
      )
    )
  }
  # get and write significant markers infos
  signif_marker_info <- data.frame(
    "SNP" = signif_snp_threshold_$SNP,
    "p_value" = as.numeric(signif(
      as.numeric(signif_snp_threshold_[, 2]), 3
    )),
    "chromosome_number" = map_df$chromosome[
      map_df$snp.name %in% signif_snp_threshold_$SNP
    ],
    "position_in_bp" = map_df$position[
      map_df$snp.name %in% signif_snp_threshold_$SNP
    ]
  )
  fwrite(signif_marker_info, file = paste0(
    paste0(output_pred_graphics_path, trait_, "/"),
    trait_, "_ls_means_gwas_significant_marker_info.csv"
  ))
}

if (compute_pca_for_mlmm_) {
  # pca population structure corrected mlmm gwas for adjusted ls-means

  # reformat response for ls_means
  ls_means_y <- ls_means_trait_[, trait_]
  names(ls_means_y) <- ls_means_trait_$Genotype

  # fit mlmm forward stepwise regression algorithm for adjusted ls-means phenotypes
  # while accounting for population structure
  mlmm_components_ls_means_ <- mlmm_allmodels_(
    ls_means_y, list(marker_df), list(k_mat),
    maxsteps = max_step_mlmm_fwd_reg_,
    cofs = geno_pca_mat_[, 1:idx_explained_var_]
  )

  # get mlmm algorithm fitted components (i.e. p-values and degrees of freedom)
  mlmm_p_val_ls_means_ <- mlmm_components_ls_means_$p_val
  mlmm_deg_freedom_ls_means_ <- mlmm_components_ls_means_$df_list

  # get step associated to "optimal" model based on ebic minimization
  res_eBIC <- as.data.frame(eBIC_allmodels(
    ls_means_y,
    frommlmm_toebic(list(marker_matrix), mlmm_p_val_ls_means_),
    list(k_mat), ncol(marker_matrix)
  ))
  optimal_model_associated_step_ <- which.min(res_eBIC$eBIC_0.75[-1])

  # make a manhattan plot
  manhattan_plot_done_ <- manhattan_plot_(
    mlmm_p_val_ls_means_,
    map = map_df,
    steps = optimal_model_associated_step_,
    output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
    file_name_ = paste0(trait_, "_pca_pop_struct_ls_means_manhattan_plot.png"),
    threshold_ = -log10(threshold_),
    main_ = paste0(
      "MLMM PCA corrected GWAS for ", trait_,
      " LS-means phenotypes for REFPOP,
    using a corrected Bonferroni threshold for 0.05"
    )
  )
  manhattan_plot_done_

  # get p-values associated with the step of the identified "optimal" model (i.e.,
  # the model minimizing the extended-BIC). Note that these p-values do not correspond
  # to the "optimal" model itself but to other models of the same form under H1 that
  # were not significant at the step of the identified "optimal" model.
  # The p-values specifically associated with the "optimal" model under H1 are those
  # corresponding to the SNPs included in this model.
  p_val_step_opt_model_ <- mlmm_p_val_ls_means_[[
    optimal_model_associated_step_ + 1
  ]]
  head(sort(p_val_step_opt_model_, decreasing = F))

  # get degrees of freedom of "optimal" model (i.e. minimizing extended-BIC). In
  # other words, get df1 = q = p1 (nb. param. under H1) - p0 (nb. param. under H0)
  # which equals 1 in case of the nullity test of q = 1 variable (SNP) under H0,
  # and df2 = n - p1 which decreases by 1 at each step (i.e. p1 increases by 1 at
  # each step).
  deg_freed_step_opt_model_ <- sort(unique(mlmm_deg_freedom_ls_means_[[
    optimal_model_associated_step_ + 1
  ]]))

  # make mlmm gwas qqplot for trait adjsuted ls-means estimated phenotypes
  qq_plot_done_ <- qqplot_f_distrib_mlmm_gwas_(
    p_val_obs_ = p_val_step_opt_model_,
    df1 = as.numeric(deg_freed_step_opt_model_[1]),
    df2 = as.numeric(deg_freed_step_opt_model_[2]),
    seed_ = 42,
    main_ = paste0(
      "MLMM PCA corrected GWAS Q-Q plot for ",
      trait_, " adjusted LS-means estimated phenotypes"
    ),
    output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
    file_name_ = paste0(trait_, "_pca_pop_struct_ls_means_qq_plot.png")
  )
  qq_plot_done_

  # get significant snp according to defined corrected Bonferroni threshold
  signif_snp_threshold_ <- as.data.frame(threshold_allmodels(
    threshold = threshold_, mlmm_p_val_ls_means_
  ))
  if (
    sum(grepl(marker_prefix_, signif_snp_threshold_$`p-value`)) > 0) {
    signif_snp_threshold_ <- correct_signif_snp_threshold(
      signif_snp_threshold_
    )
  }

  # make boxplots for trait phenotypes according to SNP allele combinations (
  # i.e. homozygote and heterozygote ) of genotypes, and write significant markers
  # infos
  if (length(signif_snp_threshold_$SNP) > 0) {
    for (snp_ in signif_snp_threshold_$SNP) {
      # get marker p-value for boxplot
      p_val_ <- signif(
        as.numeric(
          signif_snp_threshold_$p.value[
            match(snp_, signif_snp_threshold_$SNP)
          ]
        ), 3
      )
      df_ <- data.frame(
        "ls_means_y" = ls_means_y, "marker" = marker_matrix[, snp_]
      )
      r2_ <- summary(
        lm(as.formula("ls_means_y ~ 1 + marker"), data = df_)
      )$r.squared
      r2_ <- signif(as.numeric(r2_), 3)
      geno_box_plot_ <- genotypes_boxplot_(
        X = marker_matrix,
        Y = ls_means_y,
        ylab_ = paste0(trait_, " adjusted LS-means phenotypes"),
        xlab_ = paste0(
          "Genotype alleles combination on chromosome ",
          map_df$chromosome[match(snp_, map_df$snp.name)]
        ),
        markers = snp_,
        r2_ = r2_,
        p_val_ = p_val_,
        output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
        file_name_ = paste0(
          trait_, "_pca_pop_struct_ls_means_gwas_significant_marker_",
          snp_, ".png"
        )
      )
    }
    # get and write significant markers infos
    signif_marker_info <- data.frame(
      "SNP" = signif_snp_threshold_$SNP,
      "p_value" = as.numeric(signif(
        as.numeric(signif_snp_threshold_[, 2]), 3
      )),
      "chromosome_number" = map_df$chromosome[
        map_df$snp.name %in% signif_snp_threshold_$SNP
      ],
      "position_in_bp" = map_df$position[
        map_df$snp.name %in% signif_snp_threshold_$SNP
      ]
    )
    fwrite(signif_marker_info, file = paste0(
      paste0(output_pred_graphics_path, trait_, "/"),
      trait_, "_pca_pop_struct_ls_means_gwas_significant_marker_info.csv"
    ))
  }
}

# umap population structure corrected mlmm gwas for adjusted ls-means

if (compute_umap_for_mlmm_) {
  # reformat response for ls_means
  ls_means_y <- ls_means_trait_[, trait_]
  names(ls_means_y) <- ls_means_trait_$Genotype

  # fit mlmm forward stepwise regression algorithm for adjusted ls-means phenotypes
  # while accounting for population structure
  mlmm_components_ls_means_ <- mlmm_allmodels_(
    ls_means_y, list(marker_df), list(k_mat),
    maxsteps = max_step_mlmm_fwd_reg_,
    cofs = geno_umap_mat_[, 1:n_comp_umap_]
  )

  # get mlmm algorithm fitted components (i.e. p-values and degrees of freedom)
  mlmm_p_val_ls_means_ <- mlmm_components_ls_means_$p_val
  mlmm_deg_freedom_ls_means_ <- mlmm_components_ls_means_$df_list

  # get step associated to "optimal" model based on ebic minimization
  res_eBIC <- as.data.frame(eBIC_allmodels(
    ls_means_y,
    frommlmm_toebic(list(marker_matrix), mlmm_p_val_ls_means_),
    list(k_mat), ncol(marker_matrix)
  ))
  optimal_model_associated_step_ <- which.min(res_eBIC$eBIC_0.75[-1])

  # make a manhattan plot
  manhattan_plot_done_ <- manhattan_plot_(
    mlmm_p_val_ls_means_,
    map = map_df,
    steps = optimal_model_associated_step_,
    output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
    file_name_ = paste0(trait_, "_umap_pop_struct_ls_means_manhattan_plot.png"),
    threshold_ = -log10(threshold_),
    main_ = paste0(
      "MLMM UMAP corrected GWAS for ", trait_,
      " LS-means phenotypes for REFPOP,
    using a corrected Bonferroni threshold for 0.05"
    )
  )
  manhattan_plot_done_

  # get p-values associated with the step of the identified "optimal" model (i.e.,
  # the model minimizing the extended-BIC). Note that these p-values do not correspond
  # to the "optimal" model itself but to other models of the same form under H1 that
  # were not significant at the step of the identified "optimal" model.
  # The p-values specifically associated with the "optimal" model under H1 are those
  # corresponding to the SNPs included in this model.
  p_val_step_opt_model_ <- mlmm_p_val_ls_means_[[
    optimal_model_associated_step_ + 1
  ]]
  head(sort(p_val_step_opt_model_, decreasing = F))

  # get degrees of freedom of "optimal" model (i.e. minimizing extended-BIC). In
  # other words, get df1 = q = p1 (nb. param. under H1) - p0 (nb. param. under H0)
  # which equals 1 in case of the nullity test of q = 1 variable (SNP) under H0,
  # and df2 = n - p1 which decreases by 1 at each step (i.e. p1 increases by 1 at
  # each step).
  deg_freed_step_opt_model_ <- sort(unique(mlmm_deg_freedom_ls_means_[[
    optimal_model_associated_step_ + 1
  ]]))

  # make mlmm gwas qqplot for trait adjsuted ls-means estimated phenotypes
  qq_plot_done_ <- qqplot_f_distrib_mlmm_gwas_(
    p_val_obs_ = p_val_step_opt_model_,
    df1 = as.numeric(deg_freed_step_opt_model_[1]),
    df2 = as.numeric(deg_freed_step_opt_model_[2]),
    seed_ = 42,
    main_ = paste0(
      "MLMM UMAP corrected GWAS Q-Q plot for ",
      trait_, " adjusted LS-means estimated phenotypes"
    ),
    output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
    file_name_ = paste0(trait_, "_umap_pop_struct_ls_means_qq_plot.png")
  )
  qq_plot_done_

  # get significant snp according to defined corrected Bonferroni threshold
  signif_snp_threshold_ <- as.data.frame(threshold_allmodels(
    threshold = threshold_, mlmm_p_val_ls_means_
  ))
  if (
    sum(grepl(marker_prefix_, signif_snp_threshold_$`p-value`)) > 0) {
    signif_snp_threshold_ <- correct_signif_snp_threshold(
      signif_snp_threshold_
    )
  }

  # make boxplots for trait phenotypes according to SNP allele combinations (
  # i.e. homozygote and heterozygote ) of genotypes, and write significant markers
  # infos
  if (length(signif_snp_threshold_$SNP) > 0) {
    for (snp_ in signif_snp_threshold_$SNP) {
      # get marker p-value for boxplot
      p_val_ <- signif(
        as.numeric(
          signif_snp_threshold_$p.value[
            match(snp_, signif_snp_threshold_$SNP)
          ]
        ), 3
      )
      df_ <- data.frame(
        "ls_means_y" = ls_means_y, "marker" = marker_matrix[, snp_]
      )
      r2_ <- summary(
        lm(as.formula("ls_means_y ~ 1 + marker"), data = df_)
      )$r.squared
      r2_ <- signif(as.numeric(r2_), 3)
      geno_box_plot_ <- genotypes_boxplot_(
        X = marker_matrix,
        Y = ls_means_y,
        ylab_ = paste0(trait_, " adjusted LS-means phenotypes"),
        xlab_ = paste0(
          "Genotype alleles combination on chromosome ",
          map_df$chromosome[match(snp_, map_df$snp.name)]
        ),
        markers = snp_,
        r2_ = r2_,
        p_val_ = p_val_,
        output_path_ = paste0(output_pred_graphics_path, trait_, "/"),
        file_name_ = paste0(
          trait_, "_umap_pop_struct_ls_means_gwas_significant_marker_",
          snp_, ".png"
        )
      )
    }
    # get and write significant markers infos
    signif_marker_info <- data.frame(
      "SNP" = signif_snp_threshold_$SNP,
      "p_value" = as.numeric(signif(
        as.numeric(signif_snp_threshold_[, 2]), 3
      )),
      "chromosome_number" = map_df$chromosome[
        map_df$snp.name %in% signif_snp_threshold_$SNP
      ],
      "position_in_bp" = map_df$position[
        map_df$snp.name %in% signif_snp_threshold_$SNP
      ]
    )
    fwrite(signif_marker_info, file = paste0(
      paste0(output_pred_graphics_path, trait_, "/"),
      trait_, "_umap_pop_struct_ls_means_gwas_significant_marker_info.csv"
    ))
  }
}
