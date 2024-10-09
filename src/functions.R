# functions which corrects for trait spatial heterogeneity for a site
spat_hetero_env_correct_trait <- function(trait_, envir_, df_,
                                          min_obs_ = 5) {
  tryCatch(
    {
      df_envir_ <- df_ %>% subset(Envir == envir_)
      df_envir_ <- drop_na(df_envir_)
      if (nrow(df_envir_) > min_obs_) {
        df_envir_$R <- as.factor(df_envir_$Row)
        df_envir_$P <- as.factor(df_envir_$Position)
        df_envir_$Genotype <- as.factor(df_envir_$Genotype)
        spats_formula <- as.formula(~ PSANOVA(Position, Row))
        spats_model <- SpATS(
          response = trait_,
          spatial = spats_formula,
          genotype = "Genotype", genotype.as.random = TRUE, fixed = NULL,
          random = ~ R + P, data = df_envir_,
          control = controlSpATS()
        )
        spats_geno_pred <- predict(spats_model, which = "Genotype")[
          , c("Genotype", "predicted.values")
        ]
        colnames(spats_geno_pred) <- c("Genotype", "spats_geno_pred_values")
        df_envir_$spats_residuals <- spats_model$residuals
        df_envir_ <- merge(df_envir_, spats_geno_pred, by = "Genotype", all = TRUE)
        df_envir_$spats_adj_pheno <- df_envir_$spats_geno_pred_values +
          df_envir_$spats_residuals
        return(list(
          "df_envir_" = df_envir_,
          "message_" = "no error"
        ))
      } else {
        return(list(
          "df_envir_" = data.frame(),
          "message_" = "no data available for this environment"
        ))
      }
    },
    error = function(e) {
      return(list(
        "df_envir_" = data.frame(),
        "message_" = paste0(
          "Error with spat_hetero_env_correct_trait, ",
          "here is the possible issue with data or computation : ",
          conditionMessage(e)
        )
      ))
    }
  )
}

# function which re-orders column names according to a specified order
reordered_cols <- function(col_names, prefix_order_patterns) {
  tryCatch(
    {
      # extract prefixes
      prefixes_ <- sub("_[0-9]{4}$", "", col_names)

      # order according to specified order
      ordered_countries <- factor(prefixes_, levels = prefix_order_patterns)

      # order the names
      ordered_cols <- col_names[order(ordered_countries, col_names)]

      return(ordered_cols)
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function which detects dates (does not work for all dates format, be careful)
contains_date <- function(string) {
  tryCatch(
    {
      !is.na(anydate(string))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function which computes individual location clonal mean heritability
compute_indiv_location_clonal_mean_h2 <- function(lmer_mod_, nr_bar_) {
  tryCatch(
    {
      var_cov <- VarCorr(lmer_mod_)
      sigma2G <- as.numeric(attr(var_cov$Genotype, "stddev")^2)
      sigma2E <- as.numeric(sigma(lmer_mod_)^2)
      h2 <- sigma2G / (sigma2G + sigma2E / nr_bar_)
      return(h2)
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function which computes individual location clonal mean heritability
# for test purposes
compute_indiv_location_clonal_mean_h2_var_comp <- function(lmer_mod_, nr_bar_) {
  tryCatch(
    {
      var_cov <- VarCorr(lmer_mod_)
      sigma2G <- as.numeric(attr(var_cov$Genotype, "stddev")^2)
      sigma2E <- as.numeric(sigma(lmer_mod_)^2)
      h2 <- sigma2G / (sigma2G + sigma2E / nr_bar_)
      return(list(
        "h2" = h2,
        "sigma2G" = sigma2G,
        "sigma2E" = sigma2E,
        "nr_bar_" = nr_bar_
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function which computes multi-location clonal mean heritability
compute_multi_location_clonal_mean_h2 <- function(lmer_mod_, nr_bar_, nl) {
  tryCatch(
    {
      var_cov <- VarCorr(lmer_mod_)
      sigma2G <- as.numeric(attr(var_cov$Genotype, "stddev")^2)
      sigma2Gl <- as.numeric(attr(var_cov$`Genotype:Envir`, "stddev")^2)
      sigma2E <- as.numeric(sigma(lmer_mod_)^2)
      h2 <- sigma2G / (sigma2G + sigma2Gl / nl + sigma2E / (nr_bar_ * nl))
      return(h2)
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function which converts character columns to factors
convert_col_char_to_fact <- function(col) {
  if (is.character(col)) {
    return(as.factor(col))
  } else {
    return(col)
  }
}

# functions which performs imputation using column means
impute_mean <- function(x) {
  # treat NaN as NA
  x[is.nan(x)] <- NA
  mean_value <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mean_value
  return(x)
}

# function which computes the number of necessary components to reach at
# least percent_explained_variance_
n_comp_required_for_percent_explained_var <- function(facto_mine_pca_model_,
                                                      percent_explained_var) {
  return(as.numeric(which(facto_mine_pca_model_$eig[, 3] >
    percent_explained_var))[1])
}


# function which divides a vector of indices into a list of n pieces
divide_indices <- function(indices, n) {
  if (n <= 0) {
    stop("The number of pieces should be superior to n")
  }

  # compute the length of each pieces
  piece_size <- length(indices) %/% n
  remainder <- length(indices) %% n

  # create the pieces
  pieces <- vector("list", n)

  start <- 1
  for (i in 1:n) {
    # compute the end of each piece
    end <- start + piece_size - 1 + min(1, remainder)

    # save the actual piece
    pieces[[i]] <- indices[start:end]

    # update the start of the next piece
    start <- end + 1
    remainder <- max(0, remainder - 1)
  }

  return(pieces)
}

# function which tests if a matrix is a squared one
is_square <- function(mat) {
  nrow(mat) == ncol(mat)
}

# function which recodes values for biallelic markers with codes A and B
recode_values <- function(value, row_index, df_new_code) {
  if (is.na(value)) {
    return(value)
  } else {
    if ((substr(value, 1, 1) == "A") && (substr(value, 2, 2) == "A")) {
      first_elem_ <- substr(df_new_code$SNP_A[row_index], 1, 1)
      second_elem_ <- substr(df_new_code$SNP_A[row_index], 1, 1)
    }
    if ((substr(value, 1, 1) == "A") && (substr(value, 2, 2) == "B")) {
      first_elem_ <- substr(df_new_code$SNP_A[row_index], 1, 1)
      second_elem_ <- substr(df_new_code$SNP_B[row_index], 1, 1)
    }
    if ((substr(value, 1, 1) == "B") && (substr(value, 2, 2) == "A")) {
      first_elem_ <- substr(df_new_code$SNP_B[row_index], 1, 1)
      second_elem_ <- substr(df_new_code$SNP_A[row_index], 1, 1)
    }
    if ((substr(value, 1, 1) == "B") && (substr(value, 2, 2) == "B")) {
      first_elem_ <- substr(df_new_code$SNP_B[row_index], 1, 1)
      second_elem_ <- substr(df_new_code$SNP_B[row_index], 1, 1)
    }
    return(paste0(first_elem_, second_elem_))
  }
}

# function which recodes values of a matrix for biallelic markers codes as A and B
recode_matrix_values <- function(df_, df_new_code) {
  # get row indices for replacement
  row_indices <- match(rownames(df_), df_new_code[, 1])

  # create a matrix copy to save results
  col_chr_idx <- match("chr", colnames(df_))
  if (!is.na(col_chr_idx) && length(col_chr_idx) >= 1) {
    chr_list_ <- df_$chr
    df_ <- df_[, -col_chr_idx]
    result_matrix <- as.matrix(df_)
  } else {
    result_matrix <- as.matrix(df_)
  }

  # recode matrix elements
  result_matrix <- mapply(recode_values, result_matrix, row_indices,
    MoreArgs = list(df_new_code = df_new_code)
  )

  # convert matrix to data frame
  result_df <- as.data.frame(matrix(result_matrix, ncol = ncol(df_)))
  rownames(result_df) <- rownames(df_)
  colnames(result_df) <- colnames(df_)

  # add chr column
  if (!is.na(col_chr_idx) && length(col_chr_idx) >= 1) {
    result_df$chr <- chr_list_
  }

  return(result_df)
}

# function which uses median absolute deviation (MAD) to get outliers and their indices
get_outliers_vect_mad <- function(vect_, h2_mad_value_factor) {
  # apply mad to detect outlier(s)
  h2_mad_value <- mad(vect_, constant = 1)
  h2_median <- median(vect_)
  h2_threshold <- h2_mad_value_factor * h2_mad_value

  # detect environment(s) to be excluded from adjusted pheno computed h2
  idx_outliers <- which(abs(vect_ - h2_median) > h2_threshold & vect_ < h2_median)
  names_outliers <- names(vect_)[idx_outliers]

  # remove the environments to be excluded from the lists of heritabilities
  if (length(idx_outliers) > 0) {
    vect_ <- vect_[-match(names_outliers, names(vect_))]
  }

  return(list(
    "names_outliers" = names_outliers,
    "idx_outliers" = idx_outliers,
    "data_no_outliers" = vect_
  ))
}

# function which uses median absolute deviation (MAD) to get outliers and their indices
get_outliers_list_mad <- function(list_, h2_mad_value_factor) {
  # apply mad to detect outlier(s)
  mad_values <- lapply(list_, mad, constant = 1)
  medians <- sapply(list_, median)

  h2_thresholds <- lapply(1:length(list_), function(i) {
    h2_mad_value <- mad_values[[i]]
    h2_median <- medians[i]
    h2_mad_value_factor * h2_mad_value
  })

  idx_outliers <- lapply(1:length(list_), function(i) {
    abs_vals <- abs(list_[[i]] - medians[i])
    which(abs_vals > h2_thresholds[[i]] & list_[[i]] < medians[i])
  })

  names_outliers <- lapply(1:length(list_), function(i) {
    names(list_[[i]])[idx_outliers[[i]]]
  })

  # remove the environments to be excluded from the lists of heritabilities
  list_no_outliers <- lapply(1:length(list_), function(i) {
    if (length(idx_outliers[[i]]) > 0) {
      list_[[i]] <- list_[[i]][-idx_outliers[[i]]]
    }
    return(list_[[i]])
  })

  return(list(
    "names_outliers" = names_outliers,
    "idx_outliers" = idx_outliers,
    "data_no_outliers" = list_no_outliers
  ))
}

# function which replaces missing values of parents by genotype values
replace_missing_parent_by_genotype <- function(df_) {
  for (i in 1:nrow(df_)) {
    if (is.na(df_$P1[i]) || df_$P1[i] == "-") {
      df_$P1[i] <- df_$Genotype[i]
    }
    if (is.na(df_$P2[i]) || df_$P2[i] == "-") {
      df_$P2[i] <- df_$Genotype[i]
    }
  }
  return(df_)
}

# function which creates pedigree incidence matrix
create_pedig_incid_mat <- function(df_) {
  # get unique genotypes, parents, families and origin
  genotypes <- unique(df_$Genotype)
  parents <- unique(c(df_$P1, df_$P2))

  # initialize incidence data with an empty data frame
  incid_df_ <- data.frame(Genotype = genotypes)

  # add columns for parent
  for (parent in parents) {
    incid_df_[[parent]] <- 0
  }

  # fill the incidence matrix
  for (i in 1:nrow(df_)) {
    genotype <- df_$Genotype[i]
    parent_1 <- df_$P1[i]
    parent_2 <- df_$P2[i]

    # assign 1 to indicate relatedness, family membership and origin
    incid_df_[incid_df_$Genotype == genotype, parent_1] <- 1
    incid_df_[incid_df_$Genotype == genotype, parent_2] <- 1
  }

  return(incid_df_)
}

# function which splits values of a column
split_column <- function(col) {
  # split values of the form "0|0", "0|1", etc.
  split_values <- strsplit(col, "\\|")
  # create two new columns
  allele1 <- sapply(split_values, `[`, 1)
  allele2 <- sapply(split_values, `[`, 2)
  # return the two new columns
  return(data.frame(
    allele1 = as.integer(allele1),
    allele2 = as.integer(allele2)
  ))
}

# function which tunes mtry for ranger random forest
tune_mtry_ranger_rf <- function(X, Y,
                                mtry_grid_,
                                num_trees_ = 500,
                                pkgs_to_export_) {
  # initialize the cluster
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  tryCatch(
    {
      vect_mse_ <- foreach(
        mtry_ = mtry_grid_, .combine = c,
        .packages = pkgs_to_export_
      ) %dopar% {
        # build model on train data
        rf_model <- ranger(
          y = Y,
          x = X,
          mtry = mtry_,
          num.trees = num_trees_
        )
        # correlate out-of-bag (OOB) predictions (almost asymptotically
        # equivalent to LOOCV on large samples) with observed values
        mse_ <- mean((rf_model$predictions - Y)^2)
        names(mse_) <- as.character(mtry_)
        mse_
      }
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
  # stop the cluster
  stopCluster(cl)
  return(list(
    "vect_mse_" = vect_mse_,
    "opt_mtry" = as.numeric(names(which.min(vect_mse_)))
  ))
}

# function which tunes the epsilon hyperparameter for support vector regression
tune_eps_ksvm_reg_parallel <- function(X, Y, kpar_, type_, kernel_, c_par_,
                                       epsilon_grid_, n_folds_, pkgs_to_export_) {
  expected_loss_grid_ <- rep(Inf, length(epsilon_grid_))
  n <- length(Y)
  Folds <- cvFolds(n, n_folds_, type = "consecutive")

  # initialize the cluster
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  tryCatch(
    {
      expected_loss_grid_ <- foreach(
        eps_ = epsilon_grid_, .packages = pkgs_to_export_,
        .combine = "c"
      ) %dopar% {
        vect_loss_folds <- foreach(fold_ = 1:n_folds_, .combine = "c") %dopar% {
          idx_fold_ <- which(Folds$which == fold_)

          # valid set
          y_val_ <- Y[idx_fold_]
          x_val_ <- X[idx_fold_, ]

          # train set
          y_train_ <- Y[-idx_fold_]
          x_train_ <- X[-idx_fold_, ]

          # build ksvm model on train set
          ksvm_model <- ksvm(
            x = x_train_, y = y_train_,
            scaled = F, type = type_,
            kernel = kernel_,
            kpar = kpar_, C = c_par_, epsilon = eps_
          )

          f_hat_val_ <- predict(ksvm_model, x_val_)

          # loss for fold_
          sum((y_val_ - f_hat_val_)^2)
        }
        mean(vect_loss_folds)
      }
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )

  # stop the cluster
  stopCluster(cl)

  # get the optimal epsilon
  optimal_eps_ <- epsilon_grid_[which.min(expected_loss_grid_)]

  # train the "optimal" model
  tuned_model <- ksvm(
    x = X, y = Y,
    scaled = F,
    kpar = kpar_, type = type_,
    kernel = kernel_, C = c_par_,
    epsilon = optimal_eps_
  )

  return(list(
    "tuned_ksvm" = tuned_model, "optimal_eps_" = optimal_eps_,
    "expected_loss_grid_" = expected_loss_grid_,
    "epsilon_grid_" = epsilon_grid_
  ))
}

# function which removes monomorphic markers
remove_monomorphic_markers <- function(geno_df) {
  if ("Genotype" %in% colnames(geno_df)) {
    # remove genotype column
    geno_df_ <- geno_df[, -match("Genotype", colnames(geno_df))]
  } else {
    geno_df_ <- geno_df
  }
  # identify monomorphic markers
  monomorphic_markers <- apply(
    geno_df_,
    2, function(col) length(unique(col)) == 1
  )
  # get the names of the monomorphic markers
  monomorphic_marker_names <- colnames(geno_df_)[
    monomorphic_markers
  ]

  if (length(monomorphic_markers) > 0) {
    # filter the monomorphic markers
    geno_df_filtered <- geno_df_[, !monomorphic_markers]

    if ("Genotype" %in% colnames(geno_df)) {
      # add genotype column
      geno_df_filtered <- cbind(geno_df$Genotype, geno_df_filtered)
      colnames(geno_df_filtered)[1] <- "Genotype"
    }

    # return the filtered data frame and the list of monomorphic markers
    return(list(
      "filtered_df" = geno_df_filtered,
      "monomorphic_markers" = monomorphic_marker_names
    ))
  } else {
    # return the filtered data frame and the list of monomorphic markers
    return(list(
      "filtered_df" = geno_df,
      "monomorphic_markers" = NULL
    ))
  }
}

# function which removes columns with variance below a threshold
remove_low_variance_columns <- function(geno_df, threshold) {
  if ("Genotype" %in% colnames(geno_df)) {
    # remove genotype column
    geno_df_ <- geno_df[, -match("Genotype", colnames(geno_df))]
  }

  # calculate variance for each column
  variances <- apply(geno_df_, 2, var)

  # identify columns with variance below the threshold
  low_variance_columns <- variances < threshold

  # get the names of the low variance columns
  low_variance_column_names <- colnames(geno_df_)[low_variance_columns]

  # filter out low variance columns
  geno_df_filtered <- geno_df_[, !low_variance_columns]
  if ("Genotype" %in% colnames(geno_df)) {
    # add genotype column
    geno_df_filtered <- cbind(geno_df$Genotype, geno_df_filtered)
    colnames(geno_df_filtered)[1] <- "Genotype"
  }
  # return the filtered data frame and the list of low variance column names
  return(list(
    "filtered_df" = geno_df_filtered,
    "low_variance_columns" = low_variance_column_names
  ))
}

# function to calculate the weighted sum of genotype columns without weights
sum_weighted_genotypes <- function(train_genotypes, test_genotype) {
  # calculate the Hamming distance between the test genotype and
  # every genotype in the training data
  hamming_distances <- apply(train_genotypes, 1, function(row) {
    hamming.distance(test_genotype, row)
  })

  # use Hamming distance as inverse weighting
  weights <- 1 / (hamming_distances + 1)

  # weighted sum of genotypes
  weighted_sums <- colSums(weights * train_genotypes)

  # return the weighted sums of genotypes for the test genotype
  return(weighted_sums)
}

# function to reduce genotype matrix by chunks of size n
reduce_genotype_matrix <- function(genotype_matrix, n) {
  # number of columns in the genotype matrix
  n_cols <- ncol(genotype_matrix)

  # number of chunks
  n_chunks <- n_cols %/% n

  # initialize a matrix to store the reduced genotypes
  reduced_matrix <- matrix(0, nrow = nrow(genotype_matrix), ncol = n_chunks)

  # loop over each chunk
  for (i in 1:n_chunks) {
    # define the start and end columns for the current chunk
    start_col <- (i - 1) * n + 1
    end_col <- i * n

    # Reduce the chunk by summing up columns
    reduced_matrix[, i] <- rowSums(genotype_matrix[, start_col:end_col])
  }

  # return the reduced genotype matrix
  return(reduced_matrix)
}

# function which regroups and unlists results
regroup_unlist_results <- function(miss_data_singular_model_h2_out_list_) {
  # initialize a list to save results
  result_list <- list()

  tryCatch(
    {
      # iterate over trait names
      for (trait_name in names(miss_data_singular_model_h2_out_list_)) {
        # get error messages for each trait
        error_messages <- miss_data_singular_model_h2_out_list_[[trait_name]]

        # add name to each error message
        trait_error_messages <- paste(trait_name, ": ", error_messages)

        # combine results for each trait
        result_list <- c(result_list, trait_error_messages)
      }
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )

  # return unlist result_list
  return(unlist(result_list))
}

# function which remove columns with na
remove_na_columns <- function(df) {
  columns_without_na <- colSums(is.na(df)) == 0
  df_filtered <- df[, columns_without_na]
  return(df_filtered)
}

# function which compute mean of each cell in a data.frame
compute_cell_mean <- function(cell_value) {
  # unlist cell if necessary
  values <- as.numeric(unlist(cell_value))

  # compute and return mean
  mean_ <- mean(values, na.rm = TRUE)

  return(mean_)
}

# function which detects if x contains only numeric values or NA/NaN
is_numeric_with_na <- function(x) {
  all(is.na(x) | is.nan(x) | is.numeric(x))
}

# function to remove columns with a rate of NA/NaN above a threshold
remove_col_with_na_thresh <- function(df_, threshold) {
  # calculate the proportion of NA/NaN for each column
  na_ratios <- sapply(df_, function(col) {
    if (is_numeric_with_na(col)) {
      mean(is.na(col) | is.nan(col))
    } else {
      NA
    }
  })

  # identify columns to remove (NA/NaN proportion greater than the threshold)
  cols_to_remove <- names(na_ratios)[na_ratios > threshold]

  # remove these columns from the dataframe
  filtered_df <- df_[, !(colnames(df_) %in% cols_to_remove)]

  # return the filtered dataframe and the list of removed column names
  return(list(
    "filtered_df" = filtered_df,
    "removed_columns" = cols_to_remove
  ))
}

# function which computes the mode of a continuous density
compute_density_mode <- function(data) {
  if (is.numeric(data) && all(data == floor(data))) {
    mode_value <- as.numeric(names(sort(table(data), decreasing = TRUE)[1]))
  } else {
    density_estimate <- density(data)
    mode_value <- density_estimate$x[which.max(density_estimate$y)]
  }
  return(mode_value)
}

# function which computes genomic h2
compute_genomic_h2 <- function(sigma2K, sigma2E) {
  h2 <- sigma2K / (sigma2K + sigma2E)
  return(h2)
}

# function which computes the mode of a vector
compute_vect_mode <- function(vect) {
  freq <- table(vect)
  mode <- names(freq)[freq == max(freq)][1]
  return(mode)
}

# function which returns columns with at least two unique values
find_columns_with_multiple_unique_values <- function(df_) {
  cols_with_multiple_values <- colnames(df_)[apply(
    df_, 2, function(col) length(unique(col)) > 1
  )]
  return(cols_with_multiple_values)
}

# reduce dataset based on genotype counts
reduce_dataset_based_on_genotypes <- function(df_, nrow_approx_lim = 10e3,
                                              min_samples = 30) {
  K <- nrow(df_) / nrow_approx_lim
  df_reduced <- df_ %>%
    group_by(Genotype) %>%
    group_modify(~ {
      n_genotype <- nrow(.x)
      if (n_genotype < min_samples) {
        samples_to_take <- n_genotype
      } else {
        samples_to_take <- ceiling(n_genotype / K)
      }
      sample_n(.x, samples_to_take)
    }) %>%
    ungroup()

  return(df_reduced)
}

# function which reduce dataset based on selected whitening method
reduce_dataset_based_on_selected_whitening <- function(
    whitening_method,
    raw_pheno_df,
    nrow_approx_lim_raw_dataset_zca_cor,
    nrow_approx_lim_raw_dataset_pca_cor,
    nrow_approx_lim_raw_dataset_chol) {
  set.seed(123)
  if (whitening_method == "ZCA-cor") {
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_approx_lim = nrow_approx_lim_raw_dataset_zca_cor
      )
    )
  } else if (whitening_method == "PCA-cor") {
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_approx_lim = nrow_approx_lim_raw_dataset_pca_cor
      )
    )
  } else {
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_approx_lim = nrow_approx_lim_raw_dataset_chol
      )
    )
  }
  return(raw_pheno_df)
}

# function which get raw phenotypes and marker based on common genotypes
raw_pheno_and_marker_based_on_trait_common_genotypes <- function(
    raw_pheno_df,
    omic_df,
    trait_,
    fixed_effects_vars,
    random_effects_vars) {
  # remove all rows with na w.r.t to trait_
  raw_pheno_df[, trait_] <- as.numeric(raw_pheno_df[, trait_])
  raw_pheno_df <- raw_pheno_df %>% drop_na(all_of(trait_))

  # define variables of interest
  sel_vars_ <- c(fixed_effects_vars, random_effects_vars, trait_)

  # get only variables of interest from raw_pheno_df
  raw_pheno_df <- raw_pheno_df[, sel_vars_]
  raw_pheno_df <- na.omit(raw_pheno_df)

  # droplevels in order to remove levels which don't exist anymore
  raw_pheno_df <- droplevels(raw_pheno_df)

  # get phenotypes and marker data based on common genotypes
  raw_pheno_df <- match_indices(raw_pheno_df, omic_df)
  omic_df <- omic_df[rownames(omic_df) %in% unique(raw_pheno_df$Genotype), ]

  return(list(
    "raw_pheno_df" = raw_pheno_df,
    "omic_df" = omic_df
  ))
}

# function which computes fixed effect vars as factors for those declared as
compute_fixed_effect_vars_declared_as_factors <- function(
    raw_pheno_df,
    fixed_effects_vars,
    fixed_effects_vars_computed_as_factor,
    site_var,
    fixed_effects_vars_computed_as_factor_by_site) {
  # convert fixed effects variables to factors for those declared
  for (fix_eff_var_ in fixed_effects_vars) {
    if ((length(fixed_effects_vars_computed_as_factor) > 0) &&
      (fix_eff_var_ %in% fixed_effects_vars_computed_as_factor)
    ) {
      # if fix_eff_var_ is equal to management, remove buffer if exists
      if ("buffer" %in% tolower(raw_pheno_df[, fix_eff_var_])) {
        raw_pheno_df <- raw_pheno_df[
          tolower(raw_pheno_df[, fix_eff_var_]) != "buffer",
        ]
      }
      # compute fix_eff_var_ as factor by site, for those declared as,
      # otherwise compute fix_eff_var_ as a factor only
      if ((length(site_var) > 0 &&
        length(fixed_effects_vars_computed_as_factor_by_site) > 0) &&
        (fix_eff_var_ %in% fixed_effects_vars_computed_as_factor_by_site)
      ) {
        site_var_fix_eff_var_ <- paste0(
          raw_pheno_df[, site_var], "_", fix_eff_var_, "_",
          raw_pheno_df[, fix_eff_var_]
        )
        raw_pheno_df[, fix_eff_var_] <- as.factor(site_var_fix_eff_var_)
      } else {
        raw_pheno_df[, fix_eff_var_] <- as.factor(
          raw_pheno_df[, fix_eff_var_]
        )
      }
    }
  }
  return(droplevels(raw_pheno_df))
}

# function which computes the Gram matrix (i.e. genomic covariance matrix)
compute_gram_matrix <- function(omic_df, kernel_type) {
  geno_names <- rownames(omic_df)
  omic_df <- apply(omic_df, 2, as.numeric)
  if (kernel_type == "linear") {
    k_mat <- tcrossprod(scale(apply(omic_df, 2, as.numeric),
      center = T, scale = F
    ))
  } else {
    # kernel identity is not recommended due to constrained hypothesis about
    # genotypes independence which may lead to low precision
    k_mat <- as.matrix(diag(nrow(omic_df)))
  }
  # test positive definiteness and force it if necessary
  if (!is.positive.definite(k_mat, tol = 1e-8)) {
    k_mat <- as.matrix(nearPD(k_mat)$mat)
  }
  # assign genotype rownames and colnames to k_mat
  colnames(k_mat) <- rownames(k_mat) <- geno_names
  return(k_mat)
}

# function which computes incidence matrices for fixed and random effects
# NB. column of ones is added for intercept, and is associated to first
# fixed effect (which must be factor or numeric) during construction
compute_incidence_matrices_fixed_and_random_effects <- function(
    fixed_effects_vars,
    fixed_effects_vars_computed_as_factor,
    random_effects_vars,
    raw_pheno_df) {
  # define list of incidence matrices for fixed effects
  list_x_mat <- vector("list", length(fixed_effects_vars))
  names(list_x_mat) <- fixed_effects_vars

  # add incidence matrix for first fixed effect to list of matrices
  fix_eff_var_ <- fixed_effects_vars[1]
  if ((length(fixed_effects_vars_computed_as_factor) > 0) &&
    (fix_eff_var_ %in% fixed_effects_vars_computed_as_factor)
  ) {
    list_x_mat[[fix_eff_var_]] <- model.matrix(
      as.formula(paste0("~", fix_eff_var_)),
      data = raw_pheno_df
    )
    colnames(list_x_mat[[fix_eff_var_]]) <- str_replace_all(
      colnames(list_x_mat[[fix_eff_var_]]),
      pattern = fix_eff_var_, replacement = paste0(fix_eff_var_, "_")
    )
  } else {
    list_x_mat[[fix_eff_var_]] <- cbind(
      rep(1, nrow(raw_pheno_df)),
      raw_pheno_df[, fix_eff_var_]
    )
    colnames(list_x_mat[[fix_eff_var_]]) <- c("Intercept", fix_eff_var_)
  }
  # add incidence matrices (without intercept) for other fixed effects to list
  for (fix_eff_var_ in fixed_effects_vars[-1]) {
    if ((length(fixed_effects_vars_computed_as_factor) > 0) &&
      (fix_eff_var_ %in% fixed_effects_vars_computed_as_factor)
    ) {
      list_x_mat[[fix_eff_var_]] <- model.matrix(
        as.formula(paste0("~", fix_eff_var_, " - 1")),
        data = raw_pheno_df
      )
      colnames(list_x_mat[[fix_eff_var_]]) <- str_replace_all(
        colnames(list_x_mat[[fix_eff_var_]]),
        pattern = fix_eff_var_, replacement = paste0(fix_eff_var_, "_")
      )
    } else {
      list_x_mat[[fix_eff_var_]] <- raw_pheno_df[, fix_eff_var_]
      names(list_x_mat[[fix_eff_var_]]) <- fix_eff_var_
    }
  }
  x_mat <- do.call(cbind, list_x_mat)
  x_mat <- apply(x_mat, 2, as.numeric)

  # define list of incidence matrices for random effects
  list_z_mat <- vector("list", length(random_effects_vars))
  names(list_z_mat) <- random_effects_vars

  # add incidence matrices for random effects to list
  for (rand_eff_var in random_effects_vars) {
    # make sure effect is indeed a factor
    raw_pheno_df[, rand_eff_var] <- as.factor(raw_pheno_df[, rand_eff_var])
    # build incidence matrix for random effect rand_eff_var
    list_z_mat[[rand_eff_var]] <- model.matrix(
      as.formula(paste0("~", rand_eff_var, " - 1")),
      data = raw_pheno_df
    )
    colnames(list_z_mat[[rand_eff_var]]) <- str_replace_all(
      colnames(list_z_mat[[rand_eff_var]]),
      pattern = rand_eff_var, replacement = paste0(rand_eff_var, "_")
    )
  }
  z_mat <- do.call(cbind, list_z_mat)
  z_mat <- apply(z_mat, 2, as.numeric)
  return(
    list(
      "x_mat" = x_mat,
      "z_mat" = z_mat
    )
  )
}

# function which computes the whitening matrix for sig_mat_ based on the
# selected whitening method
compute_whitening_matrix_for_sig_mat_ <- function(whitening_method,
                                                  regularization_method,
                                                  parallelized_cholesky,
                                                  sig_mat_, alpha_) {
  # regularize covariance matrix, by adding a strictly positive value to the
  # diagonal of Σ, to ensure its positive definiteness
  if (regularization_method == "frobenius_norm_regularization") {
    sig_mat_ <- frobenius_norm_regularization(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "frobenius_norm_shrinkage") {
    sig_mat_ <- frobenius_norm_shrinkage(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "frobenius_norm_partial_shrinkage") {
    sig_mat_ <- frobenius_norm_partial_shrinkage(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "trace_sample_variance_regularization") {
    sig_mat_ <- trace_sample_variance_regularization(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "trace_sample_variance_shrinkage") {
    sig_mat_ <- trace_sample_variance_shrinkage(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "trace_sample_variance_partial_shrinkage") {
    sig_mat_ <- trace_sample_variance_partial_shrinkage(
      sig_mat_, alpha_
    )
  }

  # compute whitening matrix, either from ZCA-cor or cholesky
  # decomposition (i.e. Σ = LL' )
  if (whitening_method == "ZCA-cor") {
    # compute w_mat from ZCA-cor
    w_mat <- whiteningMatrix(sig_mat_, method = "ZCA-cor")
  } else if (whitening_method == "PCA-cor") {
    # compute w_mat from ZCA-cor
    w_mat <- whiteningMatrix(sig_mat_, method = "PCA-cor")
  } else {
    # compute w_mat = L^−1 from Cholesky decomposition
    if (parallelized_cholesky) {
      L <- t(cholesky(sig_mat_, parallel = T))
    } else {
      L <- t(cholesky(sig_mat_, parallel = F))
    }
    w_mat <- forwardsolve(L, diag(nrow(L)))
  }
  return(list(
    "sig_mat_" = sig_mat_,
    "w_mat" = w_mat
  ))
}

# function which match indices (not only the first one)
match_indices <- function(raw_pheno_df, omic_df) {
  common_indices <- unlist(lapply(rownames(omic_df), function(g) {
    which(raw_pheno_df$Genotype == g)
  }))
  return(raw_pheno_df[common_indices, ])
}

# function which normalizes the columns of a matrix
normalize_matrix_columns <- function(matrix_) {
  column_norms <- apply(matrix_, 2, function(col) sqrt(sum(col^2)))
  matrix_ <- sweep(matrix_, 2, column_norms, FUN = "/")
  return(matrix_)
}

# function which computes the log of determinant
log_det <- function(Sigma) {
  return(2 * sum(log(diag(cholesky(Sigma, parallel = T)))))
}

# function which computes KL divergence in a more stable way numerically
kl_divergence_ <- function(Sigma_1, Sigma_2) {
  Sigma_1 <- regularize_covariance_mean_eigenvalues(Sigma_1)
  Sigma_2 <- regularize_covariance_mean_eigenvalues(Sigma_2)

  inv_Sigma_2 <- ginv(Sigma_2)
  term_1 <- sum(diag(inv_Sigma_2 %*% Sigma_1))
  term_2 <- log_det(Sigma_2) - log_det(Sigma_1)

  k <- nrow(Sigma_1)
  kl_div <- 0.5 * (term_1 - k + term_2)

  return(kl_div)
}

# function which computes the trace of a matrix
trace_mat <- function(matrix_) {
  return(sum(diag(matrix_)))
}

# function which computes the l2 norm of a matrix
frobenius_norm <- function(matrix_) {
  return(sqrt(sum(matrix_^2)))
}

# function which makes a covariance matrix positive definite by adding
# a strictly positive diagonal matrix based on the Frobenius norm
frobenius_norm_regularization <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)
  cov_mat_ <- cov_mat_ + alpha_ * frobenius_norm(cov_mat_) * diag(n)
  return(cov_mat_)
}

# function which makes a covariance matrix positive definite by using
# a convex shrinkage estimator which adds a strictly positive diagonal matrix
# based on the Frobenius norm
frobenius_norm_shrinkage <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)
  cov_mat_ <- (1 - alpha_) * cov_mat_ + alpha_ * frobenius_norm(cov_mat_) * diag(n)
  return(cov_mat_)
}

# function which applies partial shrinkage to the diagonal elements
# of a covariance matrix by adding a positive diagonal matrix based on
# the Frobenius norm. This does not modify the off-diagonal elements.
frobenius_norm_partial_shrinkage <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)

  # create a new matrix where only the diagonal is modified
  cov_mat_reg <- cov_mat_
  diag(cov_mat_reg) <- (1 - alpha_) * diag(cov_mat_) + alpha_ * frobenius_norm(cov_mat_)

  return(cov_mat_reg)
}

# function which makes a covariance matrix positive definite by adding
# a positive diagonal matrix based on the sum of sample variances
# (i.e. trace of sample covariance matrix)
trace_sample_variance_regularization <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)
  cov_mat_ <- cov_mat_ + alpha_ * trace_mat(cov_mat_) * diag(n)
  return(cov_mat_)
}

# function which makes a covariance matrix positive definite by using
# a convex shrinkage estimator which adds a positive diagonal matrix based on
# the sum of sample variances (i.e. trace of sample covariance matrix)
trace_sample_variance_shrinkage <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)
  cov_mat_ <- (1 - alpha_) * cov_mat_ + alpha_ * trace_mat(cov_mat_) * diag(n)
  return(cov_mat_)
}

# function which applies partial shrinkage to the diagonal elements
# of a covariance matrix by adding a positive diagonal matrix based on
# the sum of sample variances (i.e. trace of sample covariance matrix).
# This does not modify the off-diagonal elements.
trace_sample_variance_partial_shrinkage <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)

  # create a new matrix where only the diagonal is modified
  cov_mat_reg <- cov_mat_
  diag(cov_mat_reg) <- (1 - alpha_) * diag(cov_mat_) + alpha_ * trace_mat(cov_mat_)

  return(cov_mat_reg)
}

# function to simulate phenotype data
simulate_y <- function(x_mat, z_mat, beta_hat, sigma2_u, sigma2_e, k_mat) {
  # get incidence matrices dimensions
  n <- nrow(x_mat)
  q <- ncol(z_mat)

  # simulate u and eps
  u <- mvrnorm(1, mu = rep(0, q), Sigma = sigma2_u * k_mat)
  eps <- rnorm(n, mean = 0, sd = sqrt(sigma2_e))

  # compute simulated y
  y_sim <- x_mat %*% beta_hat + z_mat %*% u + eps
  return(y_sim)
}

# function to compute distance between observed and simulated phenotypes
squared_l2_norm <- function(y, y_sim) {
  return(sum((y - y_sim)^2))
}

# function to simulate phenotypes and compute distance between simulated and
# observed values
simulate_and_compute_squared_l2_norm <- function(y, x_mat, z_mat, k_mat, beta_hat,
                                                 prior_sigma2_u, prior_sigma2_e) {
  # sample random values for variance components for prior ranges
  sigma2_u <- runif(1, prior_sigma2_u[1], prior_sigma2_u[2])
  sigma2_e <- runif(1, prior_sigma2_e[1], prior_sigma2_e[2])

  # simulate phenotypes
  y_sim <- simulate_y(x_mat, z_mat, beta_hat, sigma2_u, sigma2_e, k_mat)

  # compute distances
  dist_y_y_sim <- squared_l2_norm(y, y_sim)

  return(c(sigma2_u, sigma2_e, dist_y_y_sim))
}

# abc function to compute variance components
abc_variance_component_estimation <- function(y, x_mat, z_mat, k_mat, beta_hat,
                                              prior_sigma2_u, prior_sigma2_e,
                                              n_sim_abc, seed_abc,
                                              quantile_threshold_abc) {
  df_results <- future.apply::future_lapply(
    1:n_sim_abc,
    future.seed = T,
    function(sim_num) {
      set.seed(sim_num * seed_abc)
      simulate_and_compute_squared_l2_norm(
        y, x_mat, z_mat, k_mat, beta_hat, prior_sigma2_u, prior_sigma2_e
      )
    },
    future.packages = c("MASS")
  )
  df_results <- do.call(rbind, df_results)

  # assign colnames to df_results
  df_results <- as.data.frame(df_results)
  colnames(df_results) <- c("sigma2_u_hat", "sigma2_e_hat", "distance")

  # extract df_results
  vect_distances <- as.numeric(df_results[, "distance"])

  # get rejection threshold based on define quantile_threshold_abc
  reject_thresh <- quantile(vect_distances, quantile_threshold_abc)

  # get accepted variance components parameters for rejection threshold
  accepted_params <- as.data.frame(
    df_results[vect_distances <= reject_thresh, ]
  )

  # compute the average of the accepted parameters
  sigma2_u_hat_mean <- mean(accepted_params[, "sigma2_u_hat"])
  sigma2_e_hat_mean <- mean(accepted_params[, "sigma2_e_hat"])

  return(list(
    "complete_results" = df_results,
    "sigma2_u_hat_mean" = sigma2_u_hat_mean,
    "sigma2_e_hat_mean" = sigma2_e_hat_mean,
    "accepted_params" = accepted_params,
    "rejection_threshold" = reject_thresh
  ))
}

# function which computes transformed fixed variables and least squares
compute_transformed_vars_and_ols_estimates <- function(
    omic_df, raw_pheno_df, trait_,
    fixed_effects_vars,
    fixed_effects_vars_computed_as_factor,
    site_var,
    fixed_effects_vars_computed_as_factor_by_site,
    random_effects_vars,
    sigma2_u, sigma2_e, kernel_type,
    whitening_method,
    regularization_method,
    alpha_,
    parallelized_cholesky,
    reduce_raw_dataset_size_,
    nrow_approx_lim_raw_dataset_zca_cor,
    nrow_approx_lim_raw_dataset_pca_cor,
    nrow_approx_lim_raw_dataset_chol) {
  tryCatch(
    {
      # get raw phenotypes and omic data based on common genotypes
      raw_data_obj <-
        raw_pheno_and_marker_based_on_trait_common_genotypes(
          raw_pheno_df,
          omic_df,
          trait_,
          fixed_effects_vars,
          random_effects_vars
        )
      raw_pheno_df <- raw_data_obj$raw_pheno_df
      omic_df <- raw_data_obj$omic_df

      # should raw dataset size be reduced wrt to selected whitening method ?
      if (reduce_raw_dataset_size_ && (
        (whitening_method == "ZCA-cor" &&
          nrow(raw_pheno_df) > nrow_approx_lim_raw_dataset_zca_cor) ||
          (whitening_method == "PCA-cor" &&
            nrow(raw_pheno_df) > nrow_approx_lim_raw_dataset_pca_cor) ||
          (whitening_method == "Cholesky" &&
            nrow(raw_pheno_df) > nrow_approx_lim_raw_dataset_chol)
      )
      ) {
        raw_pheno_df <- reduce_dataset_based_on_selected_whitening(
          whitening_method,
          raw_pheno_df,
          nrow_approx_lim_raw_dataset_zca_cor,
          nrow_approx_lim_raw_dataset_pca_cor,
          nrow_approx_lim_raw_dataset_chol
        )
      }

      # computes fixed effect vars as factors for those declared as
      raw_pheno_df <- compute_fixed_effect_vars_declared_as_factors(
        raw_pheno_df,
        fixed_effects_vars,
        fixed_effects_vars_computed_as_factor,
        site_var,
        fixed_effects_vars_computed_as_factor_by_site
      )

      # get omic data associated to common genotypes
      omic_df <- omic_df[rownames(omic_df) %in% unique(raw_pheno_df$Genotype), ]

      # compute Gram matrix (e.g. genomic covariance matrix)
      k_mat <- compute_gram_matrix(omic_df, kernel_type)

      # retain fixed effects with variance or non unique level for factors
      if (!is.null(ncol(raw_pheno_df[, fixed_effects_vars])) &&
        ncol(raw_pheno_df[, fixed_effects_vars]) > 1) {
        fixed_effects_vars <- find_columns_with_multiple_unique_values(
          raw_pheno_df[, fixed_effects_vars]
        )
      }

      # get incidence matrices for fixed and random effects
      # NB. column of ones is added for intercept associated to fixed effects
      incid_obj <- compute_incidence_matrices_fixed_and_random_effects(
        fixed_effects_vars,
        fixed_effects_vars_computed_as_factor,
        random_effects_vars,
        raw_pheno_df
      )
      x_mat <- incid_obj$x_mat
      z_mat <- incid_obj$z_mat

      # compute Σu, i.e. sig_mat_ here
      sig_mat_ <- sigma2_u * crossprod(t(z_mat), tcrossprod(k_mat, z_mat))

      # compute the whitening matrix for Σu based on the selected
      # whitening method
      white_obj <- compute_whitening_matrix_for_sig_mat_(
        whitening_method,
        regularization_method,
        parallelized_cholesky,
        sig_mat_, alpha_
      )
      w_mat <- white_obj$w_mat
      sig_mat_ <- white_obj$sig_mat_

      # whiten x_mat using w_mat
      # NB. intercept is already present in x_mat and x_mat_tilde
      x_mat_tilde <- w_mat %*% x_mat

      # get raw phenotypes associated to common genotypes
      y <- as.numeric(raw_pheno_df[, trait_])

      # get ols estimates for fixed effects and xi
      beta_hat <- ginv(t(x_mat_tilde) %*% x_mat_tilde) %*% t(x_mat_tilde) %*% y
      y_hat <- x_mat_tilde %*% beta_hat
      xi_hat <- y - y_hat

      return(list(
        "omic_df" = omic_df,
        "sig_mat_u" = sig_mat_,
        "w_mat" = w_mat,
        "x_mat" = x_mat,
        "x_mat_tilde" = x_mat_tilde,
        "z_mat" = z_mat,
        "k_mat" = k_mat,
        "beta_hat" = beta_hat,
        "y_hat" = y_hat,
        "xi_hat" = xi_hat,
        "y" = y
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}

# function which computes phenotypes approximating genetic values using whitening
estimate_wiser_phenotype <- function(omic_df, raw_pheno_df, trait_,
                                     fixed_effects_vars = c(
                                       "Envir", "Country", "Year",
                                       "Row", "Position", "Management"
                                     ),
                                     fixed_effects_vars_computed_as_factor = c(
                                       "Envir", "Country", "Year",
                                       "Row", "Position", "Management"
                                     ),
                                     site_var = "Country",
                                     fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
                                     random_effects_vars = "Genotype",
                                     init_sigma2_u = 1,
                                     init_sigma2_e = 1,
                                     prior_sigma2_lower_bound = 1e-2,
                                     n_sim_abc = 100,
                                     seed_abc = 123,
                                     quantile_threshold_abc = 0.05,
                                     nb_iter_abc = 1,
                                     kernel_type = "linear",
                                     whitening_method = "ZCA-cor",
                                     regularization_method = "frobenius_norm_regularization",
                                     alpha_ = 0.01,
                                     parallelized_cholesky = T,
                                     reduce_raw_dataset_size_ = T,
                                     nrow_approx_lim_raw_dataset_zca_cor = 20e3,
                                     nrow_approx_lim_raw_dataset_pca_cor = 20e3,
                                     nrow_approx_lim_raw_dataset_chol = 40e3) {
  tryCatch(
    {
      # compute transformed variables associated to fixed effects and least-squares
      # to estimate these
      transform_and_ls_obj <- compute_transformed_vars_and_ols_estimates(
        omic_df, raw_pheno_df, trait_,
        fixed_effects_vars,
        fixed_effects_vars_computed_as_factor,
        site_var,
        fixed_effects_vars_computed_as_factor_by_site,
        random_effects_vars,
        sigma2_u = init_sigma2_u,
        sigma2_e = init_sigma2_e,
        kernel_type,
        whitening_method,
        regularization_method,
        alpha_,
        parallelized_cholesky,
        reduce_raw_dataset_size_,
        nrow_approx_lim_raw_dataset_zca_cor,
        nrow_approx_lim_raw_dataset_pca_cor,
        nrow_approx_lim_raw_dataset_chol
      )

      # get an upper bound for sigma2_u et sigma2_e priors
      prior_sigma2_upper_bound <- var(
        transform_and_ls_obj$y
      )

      for (iter_ in 1:nb_iter_abc) {
        # compute variance components using abc
        var_comp_abc_obj <- abc_variance_component_estimation(
          y = transform_and_ls_obj$y,
          x_mat = transform_and_ls_obj$x_mat,
          z_mat = transform_and_ls_obj$z_mat,
          k_mat = transform_and_ls_obj$k_mat,
          beta_hat = transform_and_ls_obj$beta_hat,
          prior_sigma2_u = c(
            prior_sigma2_lower_bound,
            prior_sigma2_upper_bound
          ),
          prior_sigma2_e = c(
            prior_sigma2_lower_bound,
            prior_sigma2_upper_bound
          ),
          n_sim_abc, seed_abc,
          quantile_threshold_abc
        )
        # compute variance components again with abc using new estimates
        transform_and_ls_obj <- compute_transformed_vars_and_ols_estimates(
          omic_df, raw_pheno_df, trait_,
          fixed_effects_vars,
          fixed_effects_vars_computed_as_factor,
          site_var,
          fixed_effects_vars_computed_as_factor_by_site,
          random_effects_vars,
          sigma2_u = var_comp_abc_obj$sigma2_u_hat_mean,
          sigma2_e = var_comp_abc_obj$sigma2_e_hat_mean,
          kernel_type,
          whitening_method,
          regularization_method,
          alpha_,
          parallelized_cholesky,
          reduce_raw_dataset_size_,
          nrow_approx_lim_raw_dataset_zca_cor,
          nrow_approx_lim_raw_dataset_pca_cor,
          nrow_approx_lim_raw_dataset_chol
        )
      }

      # get estimated and/or modified components after abc

      # extract marker data in case of modification
      omic_df <- transform_and_ls_obj$omic_df

      # extract estimated fixed effects
      beta_hat <- transform_and_ls_obj$beta_hat

      # compute phenotypic values using ols
      v_hat <- ginv(t(transform_and_ls_obj$z_mat) %*% transform_and_ls_obj$z_mat) %*%
        t(transform_and_ls_obj$z_mat) %*% transform_and_ls_obj$xi_hat

      # save wiser fixed effects estimates in a data frame
      wiser_pheno_df <- data.frame(
        "Genotype" = str_replace_all(
          colnames(transform_and_ls_obj$z_mat),
          pattern = "Genotype_",
          replacement = ""
        ),
        "v_hat" = v_hat
      )

      # save wiser phenotypes in a data frame
      wiser_fix_eff_df <- data.frame(
        "fixed_effect_var" = colnames(transform_and_ls_obj$x_mat),
        "beta_hat_var" = beta_hat
      )

      return(list(
        "wiser_omic_data" = omic_df,
        "sig_mat_u" = transform_and_ls_obj$sig_mat_,
        "w_mat" = transform_and_ls_obj$w_mat,
        "wiser_fixed_effect_estimates" = wiser_fix_eff_df,
        "wiser_abc_variance_component_estimates" = var_comp_abc_obj,
        "wiser_phenotypes" = wiser_pheno_df,
        "wiser_z_mat" = transform_and_ls_obj$z_mat,
        "wiser_x_mat" = transform_and_ls_obj$x_mat,
        "wiser_x_mat_tilde" = transform_and_ls_obj$x_mat_tilde,
        "wiser_xi_hat" = transform_and_ls_obj$xi_hat,
        "wiser_y_hat" = transform_and_ls_obj$y_hat,
        "wiser_y" = transform_and_ls_obj$y
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}


# function which performs parallelized k-folds cv using several prediction methods
perform_kfold_cv_wiser <- function(omic_df, raw_pheno_df, trait_,
                                   fixed_effects_vars,
                                   fixed_effects_vars_computed_as_factor,
                                   site_var,
                                   fixed_effects_vars_computed_as_factor_by_site,
                                   random_effects_vars,
                                   whitening_method,
                                   reg_method, alpha_,
                                   pred_method, k_folds,
                                   wiser_obj_local) {
  # extract the local wiser object
  omic_df <- wiser_obj_local$wiser_omic_data
  v_hat <- wiser_obj_local$wiser_phenotypes$v_hat

  # set seed for reproducibility and get set of indices
  set.seed(123)
  idx <- 1:nrow(omic_df)

  # create folds for k-folds cv
  folds <- cvFolds(nrow(omic_df), K = k_folds, type = "consecutive")

  # use future_lapply for folds
  results <- future_lapply(1:k_folds,
    future.seed = T,
    function(fold) {
      idx_train <- idx[folds$which != fold]
      idx_val <- idx[folds$which == fold]

      # train and predict with random forest (using ranger package)
      if (pred_method == "rf") {
        rf_model <- ranger(
          y = v_hat[idx_train],
          x = omic_df[idx_train, ],
          mtry = ncol(omic_df) / 3,
          num.trees = 1000
        )
        f_hat_val_rf <- predict(rf_model, omic_df[idx_val, ])
        pa_ <- cor(
          f_hat_val_rf$predictions,
          v_hat[idx_val]
        )
        # train and predict with non-linear svr (using kernlab package)
      } else if (pred_method == "svr") {
        c_par <- max(
          abs(mean(v_hat[idx_train])
          + 3 * sd(v_hat[idx_train])),
          abs(mean(v_hat[idx_train])
          - 3 * sd(v_hat[idx_train]))
        )
        gaussian_svr_model <- ksvm(
          x = as.matrix(omic_df[idx_train, ]),
          y = v_hat[idx_train],
          scaled = FALSE, type = "eps-svr",
          kernel = "rbfdot",
          kpar = "automatic", C = c_par, epsilon = 0.1
        )
        f_hat_val_gaussian_svr <- predict(
          gaussian_svr_model,
          as.matrix(omic_df[idx_val, ])
        )
        pa_ <- cor(
          f_hat_val_gaussian_svr,
          v_hat[idx_val]
        )
        # train and predict with gblup (using KRMM package)
      } else if (pred_method == "gblup") {
        linear_krmm_model <- krmm(
          Y = v_hat[idx_train],
          Matrix_covariates = omic_df[idx_train, ],
          method = "GBLUP"
        )
        f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
          Matrix_covariates = omic_df[idx_val, ],
          add_fixed_effects = T
        )
        pa_ <- cor(
          f_hat_val_linear_krmm,
          v_hat[idx_val]
        )
        # train and predict with rkhs (using KRMM package)
      } else if (pred_method == "rkhs") {
        gaussian_krmm_model <- krmm(
          Y = v_hat[idx_train],
          Matrix_covariates = omic_df[idx_train, ],
          method = "RKHS", kernel = "Gaussian",
          rate_decay_kernel = 0.1
        )
        f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
          Matrix_covariates = omic_df[idx_val, ],
          add_fixed_effects = T
        )
        pa_ <- cor(
          f_hat_val_gaussian_krmm,
          v_hat[idx_val]
        )
        # train and predict with lasso (using glmnet package)
      } else {
        cv_fit_lasso_model <- cv.glmnet(
          intercept = TRUE, y = v_hat[idx_train],
          x = as.matrix(omic_df[idx_train, ]),
          type.measure = "mse", alpha = 1.0, nfold = 10,
          parallel = TRUE
        )
        f_hat_val_lasso <- predict(cv_fit_lasso_model,
          newx = as.matrix(omic_df[idx_val, ]),
          s = "lambda.min"
        )
        pa_ <- suppressWarnings(cor(
          f_hat_val_lasso,
          v_hat[idx_val]
        ))
      }
      data.frame(pa = pa_)
    }, future.packages = c(
      "ranger", "KRMM", "kernlab",
      "glmnet", "cvTools", "dplyr",
      "stringr", "matrixcalc",
      "Matrix", "whitening", "mixOmics"
    )
  )

  df_results <- do.call(rbind, results)
  mean_pa <- mean(df_results$pa, na.rm = T)
  return(mean_pa)
}

# function which finds the optimal whitening method and regularization
optimize_whitening_and_regularization <- function(
    omic_df, raw_pheno_df, trait_,
    fixed_effects_vars = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    fixed_effects_vars_computed_as_factor = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    site_var = "Country",
    fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
    random_effects_vars = "Genotype",
    prediction_method = c("rf", "svr", "gblup", "rkhs", "lasso"),
    whitening_method_grid = c("ZCA-cor", "PCA-cor", "Cholesky"),
    regularization_method_ = "frobenius_norm_regularization",
    alpha_grid = c(0.01, 0.1),
    reduce_raw_dataset_size_ = T,
    nrow_approx_lim_raw_dataset_ = 5e3,
    parallelized_cholesky = T,
    k_folds_ = 5,
    nb_cores_ = 12) {
  # remove all rows with na w.r.t to trait_
  raw_pheno_df <- raw_pheno_df %>% drop_na(all_of(trait_))

  # get raw phenotypes and marker data based on common genotypes
  raw_pheno_df <- match_indices(raw_pheno_df, omic_df)
  omic_df <- omic_df[rownames(omic_df) %in% raw_pheno_df$Genotype, ]

  # define variables of interest
  sel_vars_ <- c(fixed_effects_vars, random_effects_vars, trait_)

  # get only variables of interest from raw_pheno_df
  raw_pheno_df <- raw_pheno_df[, sel_vars_]
  raw_pheno_df <- na.omit(raw_pheno_df)

  # downsize dataset for computation time optimization purpose
  if (reduce_raw_dataset_size_ &&
    (nrow(raw_pheno_df) > nrow_approx_lim_raw_dataset_)) {
    set.seed(123)
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_approx_lim = nrow_approx_lim_raw_dataset_
      )
    )
  }

  # create a grid combining whitening methods, regularization parameter and
  # prediction methods
  grid_ <- expand.grid(
    whitening_method = whitening_method_grid,
    alpha_ = alpha_grid,
    pred_method = prediction_method
  )

  # pre-compute unique wiser object for unique combinations
  wiser_cache <- list()
  unique_combinations <- unique(grid_[, c("whitening_method", "alpha_")])

  for (j in 1:nrow(unique_combinations)) {
    method <- unique_combinations$whitening_method[j]
    alpha <- unique_combinations$alpha_[j]

    wiser_obj <- estimate_wiser_phenotype(
      omic_df, raw_pheno_df, trait_,
      fixed_effects_vars,
      fixed_effects_vars_computed_as_factor,
      site_var,
      fixed_effects_vars_computed_as_factor_by_site,
      random_effects_vars,
      whitening_method = method,
      regularization_method = regularization_method_,
      alpha_ = alpha,
      reduce_raw_dataset_size_ = FALSE
    )

    cache_key <- paste(method, alpha, sep = "_")
    wiser_cache[[cache_key]] <- wiser_obj
  }

  # configure parallelization
  plan(multisession, workers = nb_cores_)

  df_results <- future_lapply(
    1:nrow(grid_),
    future.seed = T,
    function(i) {
      # retrieve the precomputed wiser_obj for this combination
      cache_key <- paste(grid_$whitening_method[i], grid_$alpha_[i], sep = "_")
      wiser_obj_local <- wiser_cache[[cache_key]]

      mean_pa <- tryCatch(
        {
          perform_kfold_cv_wiser(
            omic_df, raw_pheno_df, trait_,
            fixed_effects_vars,
            fixed_effects_vars_computed_as_factor,
            site_var,
            fixed_effects_vars_computed_as_factor_by_site,
            random_effects_vars,
            whitening_method = grid_$whitening_method[i],
            reg_method = regularization_method_,
            alpha_ = grid_$alpha_[i],
            pred_method = grid_$pred_method[i],
            k_folds = k_folds_,
            wiser_obj_local = wiser_obj_local
          )
        },
        error = function(e) {
          cat("Error during iteration", i, ":", e$message, "\n")
          return(NA)
        }
      )
      data.frame(
        "whitening_method" = grid_$whitening_method[i],
        "alpha_" = grid_$alpha_[i],
        "prediction_method" = grid_$pred_method[i],
        "mean_pa" = mean_pa
      )
    }, future.packages = c(
      "ranger", "KRMM", "kernlab",
      "glmnet", "cvTools", "dplyr",
      "stringr", "matrixcalc",
      "Matrix", "whitening", "mixOmics"
    )
  )
  df_results <- na.omit(do.call(rbind, df_results))

  # get optimal whitening method based on mean pa for each prediction method
  df_opt_ <- data.frame()
  for (method_ in df_results$prediction_method) {
    df_res_method_ <- df_results[
      df_results$prediction_method == method_,
    ]
    df_res_method_$white_reg_combination <- paste0(
      df_res_method_$whitening_method,
      "/", df_res_method_$alpha_
    )
    df_opt_ <- rbind(
      df_opt_,
      unique(df_res_method_[
        which.max(df_res_method_$mean_pa)[1],
      ])
    )
  }
  opt_mode_ <- compute_vect_mode(df_opt_$white_reg_combination)
  opt_mode_ <- unlist(str_split(opt_mode_, pattern = "/"))
  opt_whitening_method <- opt_mode_[1]
  opt_alpha_ <- as.numeric(opt_mode_[2])

  # stop parallelization
  plan(sequential)

  return(list(
    "opt_results" = unique(df_opt_),
    "opt_alpha_" = opt_alpha_,
    "opt_whitening_method" = opt_whitening_method
  ))
}

# function which make a scatter plot for ls-mean and wiser with a linear fit
create_scatter_plot_with_linear_fit <- function(df_, trait_, env_) {
  # fit simple linear regression
  fit <- lm(v_hat ~ get(trait_), data = df_)
  intercept <- coef(fit)[1]
  slope <- coef(fit)[2]
  r_squared <- summary(fit)$r.squared

  # format data
  intercept_formatted <- round(intercept, 2)
  slope_formatted <- round(slope, 2)
  r_squared_formatted <- round(r_squared, 2)

  # create equation for text
  eqn_text <- paste0(
    "Y = ", slope_formatted, "X ",
    ifelse(intercept_formatted >= 0, "+ ", "- "),
    abs(intercept_formatted),
    "\nR² = ", r_squared_formatted,
    "\nPearson r = ", signif(cor(
      df_[[trait_]],
      df_$v_hat
    ), 2)
  )

  # set annotation positions
  x_annotation <- min(df_[[trait_]]) +
    (max(df_[[trait_]]) - min(df_[[trait_]])) * 0.05
  y_annotation <- max(df_$v_hat) -
    (max(df_$v_hat) - min(df_$v_hat)) * 0.1

  # make scatter plot with linear fit
  fig_x_y <- plot_ly(
    data = df_
  ) %>%
    add_trace(
      x = ~ get(trait_),
      y = ~v_hat,
      type = "scatter",
      mode = "markers",
      marker = list(color = "blue"),
      name = "Data points"
    ) %>%
    add_lines(
      x = ~ get(trait_),
      y = fitted(fit),
      line = list(color = "red"),
      name = "Linear fit"
    ) %>%
    layout(
      title = list(
        text = paste0(
          trait_, " LS-means versus WISER phenotypes for ",
          env_
        ),
        x = 0.5
      ),
      xaxis = list(
        title = "Genotype LS-means phenotype"
      ),
      yaxis = list(
        title = "Genotype WISER phenotype"
      ),
      annotations = list(
        x = x_annotation,
        y = y_annotation,
        text = eqn_text,
        showarrow = FALSE,
        xanchor = "left",
        yanchor = "top",
        font = list(
          size = 12,
          color = "black"
        ),
        align = "left",
        bgcolor = "rgba(255, 255, 255, 0.8)",
        bordercolor = "black",
        borderwidth = 1
      )
    )
  return(fig_x_y)
}
