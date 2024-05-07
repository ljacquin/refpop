# script meant to perform genomic prediction and analyses for refpop
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(MASS)
library(data.table)
library(stringr)
library(FactoMineR)
library(rstudioapi)
library(doParallel)
library(doRNG)
library(robustbase)
library(foreach)
library(parallel)
library(missForest)
library(Matrix)
library(rgl)
library(cvTools)
library(ggplot2)
library(plotly)
library(ranger)
library(htmlwidgets)
library(dplyr)
library(reticulate)
library(devtools)
if ("refpop_geno_pred_env" %in% conda_list()$name) {
  use_condaenv("refpop_geno_pred_env")
}
# install other requirements from github if necessary
install_other_requirements <- F
if (install_other_requirements) {
  # reticulate::install_miniconda()
  conda_create("refpop_geno_pred_env")
  use_condaenv("refpop_geno_pred_env")
  devtools::install_github("ljacquin/KRMM")
  devtools::install_github("rstudio/tensorflow")
  library(tensorflow)
  install_tensorflow(envname = "refpop_geno_pred_env")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
  install.packages("umap")
}
library(KRMM)
library(kernlab)
library(tensorflow)
library(keras3)
library(umap)
tensorflow::tf$random$set_seed(0)
py_module_available("keras") # must return TRUE
py_module_available("tensorflow") # must return TRUE
py_discover_config("keras") # more info on the python env, tf and keras
setwd(dirname(getActiveDocumentContext()$path)) # set automated wd detection for script
source("../../functions.R")

# define function(s) and package(s) to export for parallelization
pkgs_to_export_ <- c(
  "ranger",
  "kernlab",
  "KRMM",
  "foreach"
)

# set paths
geno_dir_path <- "../../data/genotype_data/"
pheno_dir_path <- "../../data/phenotype_data/"
output_pred_results_path <- "../../data/genomic_pred_result_data/"
output_pred_graphics_path <- "../../data/graphics/genomic_pred_result_graphics/"

# define trait_
trait_ <- "Harvest_date" # "Flowering_full" # "Harvest_date"

# define shift seed value by
shift_seed_by_ <- 10

# define number of shuffles
n_shuff_ <- 5

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

# get phenotype and genotype data, and family and origin info
pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "adjusted_ls_means_phenotypes.csv"
)))

geno_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_data.csv"
)))

geno_fam_orig_df <- as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_family_origin_information.csv"
)))

# remove monomorphic markers
geno_df <- remove_monomorphic_markers(geno_df)
monomorphic_markers_list_ <- geno_df$monomorphic_markers
geno_df <- geno_df$filtered_df

# set genotypes in the same order between pheno and geno data
set.seed(123)
snp_sample_size_ <- 50e3
pheno_df <- pheno_df[match(geno_df$Genotype, pheno_df$Genotype), ]
idx_snp_sample_size_ <- sample(2:ncol(geno_df),
  size = snp_sample_size_, replace = F
)
geno_df_ <- geno_df[, idx_snp_sample_size_]

# remove na for analyzed trait_ and corresponding rows for marker data
idx_na_trait_ <- which(is.na(pheno_df[, trait_]))
if (length(idx_na_trait_) > 0) {
  pheno_df <- pheno_df[-idx_na_trait_, ]
  geno_df_ <- geno_df_[-idx_na_trait_, ]
}

# create indices for train set
n <- nrow(pheno_df)

# tune hyperparameters for gaussian svr and krmm, based on a small random sample,
# for original variables
set.seed(1234)
idx_tune_hyperpara <- sample(1:n, size = floor(1 / 3 * n), replace = F)

# original variables

# tune mtry parameter for ranger random forest
mtry_grid_ <- floor(seq(floor(ncol(geno_df_) / 10), floor(ncol(geno_df_) / 3),
  length.out = 4
))
list_opt_mtry_ <- tune_mtry_ranger_rf(
  X = geno_df_[idx_tune_hyperpara, ],
  Y = pheno_df[idx_tune_hyperpara, trait_],
  mtry_grid_,
  num_trees_ = 1000,
  pkgs_to_export_
)
plot(as.numeric(names(list_opt_mtry_$vect_acc_)),
  as.numeric(list_opt_mtry_$vect_acc_),
  type = "l"
)
opt_mtry <- list_opt_mtry_$opt_mtry

# tune gaussian svr on original variables

# A correct value for c_par according to Cherkassy and Ma (2004).
# Neural networks 17, 113-126
c_par <- max(
  abs(mean(pheno_df[idx_tune_hyperpara, trait_]) +
    3 * sd(pheno_df[idx_tune_hyperpara, trait_])),
  abs(mean(pheno_df[idx_tune_hyperpara, trait_])
  - 3 * sd(pheno_df[idx_tune_hyperpara, trait_]))
)
tune_gaussian_svr_model <- tune_eps_ksvm_reg_parallel(
  X = apply(geno_df_[idx_tune_hyperpara, ], 2, as.numeric),
  Y = pheno_df[idx_tune_hyperpara, trait_],
  kpar_ = "automatic", type_ = "eps-svr",
  kernel_ = "rbfdot", c_par_ = c_par,
  epsilon_grid_ = seq(1e-1, 1, length.out = 5),
  n_folds_ = 5, pkgs_to_export_
)
plot(tune_gaussian_svr_model$epsilon_grid_,
  tune_gaussian_svr_model$expected_loss_grid_,
  type = "l"
)
opt_eps <- tune_gaussian_svr_model$optimal_eps_

# tune gaussian krmm on original variables
tune_gaussian_krmm_model <- tune_krmm(
  Y = pheno_df[idx_tune_hyperpara, trait_],
  Matrix_covariates = geno_df_[idx_tune_hyperpara, ],
  rate_decay_grid = seq(1e-1, 1, length.out = 5),
  nb_folds = 5, method = "RKHS", kernel = "Gaussian"
)
plot(tune_gaussian_krmm_model$rate_decay_grid,
  tune_gaussian_krmm_model$mean_loss_grid,
  type = "l"
)
opt_h <- tune_gaussian_krmm_model$optimal_h

with_progress({
# detect number of cores and make cluster
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# compute parallelized genomic predictions and results for n_shuff_
df_result_ <- foreach(
  shuff_ = 1:n_shuff_, .packages = pkgs_to_export_,
  .combine = rbind
) %dopar% {
  set.seed(shuff_ + shift_seed_by_)
  print(paste0("shuff_ : ", shuff_))
  idx_train <- sample(1:n, size = floor(2 / 3 * n), replace = FALSE)

  # initialize vector of results for all tested models
  result <- c(
    "RF" = NA, "SVR" = NA, "RKHS" = NA, "GBLUP" = NA,
    "SVR_support_vectors" = NA
  )

  # genomic prediction based on original variables

  # train random forest on original variables
  rf_model <- ranger(
    y = pheno_df[idx_train, trait_],
    x = geno_df_[idx_train, ],
    mtry = opt_mtry,
    num.trees = 1000
  )
  # make prediction for test data and save results
  f_hat_test_rf <- predict(
    rf_model,
    geno_df_[-idx_train, ]
  )
  result[["RF"]] <- cor(
    f_hat_test_rf$predictions,
    pheno_df[-idx_train, trait_]
  )

  # train SVR with non linear gaussian kernel on original variables

  # A correct value for c_par according to Cherkassy and Ma (2004).
  # Neural networks 17, 113-126
  c_par <- max(
    abs(mean(pheno_df[idx_train, trait_]) + 3 * sd(pheno_df[idx_train, trait_])),
    abs(mean(pheno_df[idx_train, trait_]) - 3 * sd(pheno_df[idx_train, trait_]))
  )
  gaussian_svr_model <- ksvm(
    x = apply(geno_df_[idx_train, ], 2, as.numeric),
    y = pheno_df[idx_train, trait_],
    scaled = F, type = "eps-svr",
    kernel = "rbfdot",
    kpar = "automatic", C = c_par, epsilon = opt_eps
  )
  # get support vectors sv_ for gaussian svr model
  idx_sv_ <- attributes(gaussian_svr_model)$SVindex
  sv_ <- pheno_df[idx_train, "Genotype"][idx_sv_]

  # make prediction for test data and save results
  f_hat_test_gaussian_svr <- predict(
    gaussian_svr_model,
    apply(geno_df_[-idx_train, ], 2, as.numeric)
  )
  result[["SVR"]] <- cor(
    f_hat_test_gaussian_svr,
    pheno_df[-idx_train, trait_]
  )
  result[["SVR_support_vectors"]] <- paste0(sv_, collapse = ", ")

  # train krmm with linear kernel (i.e. dot product, GBLUP or RR-BLUP method)
  # on original variables
  linear_krmm_model <- krmm(
    Y = pheno_df[idx_train, trait_],
    Matrix_covariates = geno_df_[idx_train, ],
    method = "GBLUP"
  )
  # make prediction for test data and save results
  f_hat_test_linear_krmm <- predict_krmm(linear_krmm_model,
    Matrix_covariates = geno_df_[-idx_train, ],
    add_flxed_effects = T
  )
  result[["GBLUP"]] <- cor(
    f_hat_test_linear_krmm,
    pheno_df[-idx_train, trait_]
  )

  # train krmm with non linear gaussian kernel (i.e. RKHS method)
  # on original variables
  gaussian_krmm_model <- krmm(
    Y = pheno_df[idx_train, trait_],
    Matrix_covariates = geno_df_[idx_train, ],
    method = "RKHS", kernel = "Gaussian",
    rate_decay_kernel = opt_h
  )
  # make prediction for test data and save results
  f_hat_test_gaussian_krmm_model <- predict_krmm(gaussian_krmm_model,
    Matrix_covariates = geno_df_[-idx_train, ],
    add_flxed_effects = T
  )
  result[["RKHS"]] <- cor(
    f_hat_test_gaussian_krmm_model,
    pheno_df[-idx_train, trait_]
  )

  return(result)
}

# stop cluster
stopCluster(cl)
})

# create directory for trait_ graphics if it does not exist
if( !dir.exists(paste0(output_pred_graphics_path, trait_, '/')) ){
  dir.create(paste0(output_pred_graphics_path, trait_, '/'))
}

# convert to data frame format
df_result_ <- as.data.frame(df_result_)
df_ <- df_result_[, -match(
  "SVR_support_vectors",
  colnames(df_result_)
)]

# get methods names
df_ <- as.data.frame(apply(df_, 2, as.numeric))
method_names <- colnames(df_)

# initialize plot_ly boxplot graphic
boxplots_rpa_ <- plot_ly()

# add boxplots
for (method_ in method_names) {
  boxplots_rpa_ <- add_boxplot(
    boxplots_rpa_,
    y = df_[[method_]],
    name = method_,
    boxpoints = "all",
    jitter = 0.3,
    pointpos = -1.8
  )
}
# add layout
boxplots_rpa_ <- boxplots_rpa_ %>%
  layout(
    title = paste0(
      "Genomic prediction RPA distributions of methods for ",
      trait_, ", based on ", snp_sample_size_, " SNPs across ",
      n_shuff_, " shuffling scenarios"
    ),
    yaxis = list(title = "Relative Prediction Accuracy (RPA)", range = c(0, 1)),
    legend = list(title = list(text = "Prediction method"))
  )

# save boxplots_rpa_ graphics
saveWidget(boxplots_rpa_, file = paste0(
  output_pred_graphics_path, trait_, "/relative_prediction_accuracy_",
  trait_, "_", snp_sample_size_, "_SNP", ".html"
))

# save relative prediction accuracy results
rownames(df_) <- paste0("rpa_shuff_", 1:nrow(df_))
df_stat <- as.data.frame(rbind(apply(df_, 2, mean), apply(df_, 2, sd)))
rownames(df_stat) <- c("rpa_mean", "rpa_sd")
df_ <- rbind(df_, df_stat)
fwrite(df_,
  file = paste0(
    output_pred_results_path,
    "genomic_pred_results_", snp_sample_size_, "_SNP_",
    trait_, ".csv"
  ), row.names = T
)

# get genotypes used as support vectors
sv_geno <- unlist(str_split(df_result_$SVR_support_vectors,
  pattern = ", "
))

# create a frequency table
freq_table <- table(sv_geno)

# convert the frequency table to a dataframe
df_sv_ <- data.frame(
  Genotype = names(freq_table),
  Count = as.numeric(freq_table)
)
# merge genotype with family and origin info
df_sv_ <- merge(df_sv_, geno_fam_orig_df, by = "Genotype", all = T)
df_sv_$Family <- factor(df_sv_$Family, levels = unique(df_sv_$Family))
df_sv_$Origin <- factor(df_sv_$Origin, levels = unique(df_sv_$Origin))

# calculate the sum of counts for each origin
origins_sum <- df_sv_ %>%
  group_by(Origin) %>%
  summarize(total_count_origin = sum(Count))

# calculate the sum of counts for each family
families_sum <- df_sv_ %>%
  group_by(Family) %>%
  summarize(total_count_family = sum(Count))

# calculate the total Count for percentage by origin and by family
total_count <- sum(df_sv_$Count)

# join the original data with the sums calculated by origin and by family
df_sv_ <- df_sv_ %>%
  left_join(origins_sum, by = "Origin") %>%
  left_join(families_sum, by = "Family") %>%
  # calculate the percentage by origin
  mutate(Percentage_Origin = (total_count_origin / total_count) * 100) %>%
  # calculate the percentage by family
  mutate(Percentage_Family = (total_count_family / total_count) * 100)

# create barplot_sv_orig with colored bars based on origins
barplot_sv_orig <- plot_ly(df_sv_,
  x = ~Genotype, y = ~Count, type = "bar",
  color = ~Origin,
  colors = color_palette_origin
) %>%
  layout(
    title = paste0(
      "Counts of genotypes identified as support vectors across ",
      n_shuff_, " Gaussian SVR models, \n used for genomic prediction of ",
      trait_, " across ", n_shuff_, " shuffling scenarios"
    ),
    xaxis = list(
      categoryorder = "total descending", title = "Genotype",
      tickangle = 300
    ),
    yaxis = list(title = "Count"),
    legend = list(title = list(text = "Origin"))
  )
# save barplot_sv_orig graphics
saveWidget(barplot_sv_orig, file = paste0(
  output_pred_graphics_path, trait_, "/support_vector_genotypes_for_",
  trait_, "_", snp_sample_size_, "_SNP_origin_as_label", ".html"
))

# create barplot_sv_percent_orig based on percentages of origin
df_percent_origin <- unique(df_sv_[, c("Origin", "Percentage_Origin")])
barplot_sv_percent_orig <- plot_ly(df_percent_origin,
  x = ~Origin, y = ~Percentage_Origin, type = "bar",
  color = ~Origin,
  colors = color_palette_origin
) %>%
  layout(
    title = paste0(
      "Percentage breakdown for origins of genotypes, identified as support vectors across ",
      n_shuff_, " Gaussian SVR models, \n used for genomic prediction of ", trait_,
      " across ", n_shuff_, " shuffling scenarios"
    ),
    xaxis = list(
      categoryorder = "total descending", title = "Origin",
      tickangle = 300
    ),
    yaxis = list(title = "Percentage"),
    legend = list(title = list(text = "Origin"))
  )
# save barplot_sv_percent_orig graphics
saveWidget(barplot_sv_percent_orig, file = paste0(
  output_pred_graphics_path, trait_, 
  "/percentage_breakdown_for_origins_of_support_vectors_",
  trait_, "_", snp_sample_size_, "_SNP", ".html"
))

# create barplot_sv_fam with colored bars based on families
barplot_sv_fam <- plot_ly(df_sv_,
  x = ~Genotype, y = ~Count, type = "bar",
  color = ~Family,
  colors = color_palette_family
) %>%
  layout(
    title = paste0(
      "Counts of genotypes identified as support vectors across ",
      n_shuff_, " Gaussian SVR models, \n used for genomic prediction of ",
      trait_, " across ", n_shuff_, " shuffling scenarios"
    ),
    xaxis = list(
      categoryorder = "total descending", title = "Genotype",
      tickangle = 300
    ),
    yaxis = list(title = "Count"),
    legend = list(title = list(text = "Family (except accession)"))
  )
# save barplot_sv_fam graphics
saveWidget(barplot_sv_fam, file = paste0(
  output_pred_graphics_path, trait_, "/support_vector_genotypes_for_",
  trait_, "_", snp_sample_size_, "_SNP_family_as_label", ".html"
))

# create barplot_sv_percent_fam based on percentages of families
df_percent_family <- unique(df_sv_[, c("Family", "Percentage_Family")])
barplot_sv_percent_fam <- plot_ly(df_percent_family,
  x = ~Family, y = ~Percentage_Family,
  type = "bar",
  color = ~Family,
  colors = color_palette_family
) %>%
  layout(
    title = paste0(
      "Percentage breakdown for families of genotypes, identified as support vectors across ",
      n_shuff_, " Gaussian SVR models, \n used for genomic prediction of ", trait_,
      " across ", n_shuff_, " shuffling scenarios"
    ),
    xaxis = list(
      categoryorder = "total descending", title = "Family",
      tickangle = 300
    ),
    yaxis = list(title = "Percentage"),
    legend = list(title = list(text = "Family"))
  )
# save barplot_sv_percent_fam graphics
saveWidget(barplot_sv_percent_fam, file = paste0(
  output_pred_graphics_path, trait_, 
  "/percentage_breakdown_for_families_of_support_vectors_",
  trait_, "_", snp_sample_size_, "_SNP", ".html"
))
