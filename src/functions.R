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

# function which compute average of each cell in a data.frame
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

# function which normalizes the columns of a matrix
normalize_matrix_columns <- function(matrix_) {
  column_norms <- apply(matrix_, 2, function(col) sqrt(sum(col^2)))
  matrix_ <- sweep(matrix_, 2, column_norms, FUN = "/")
  return(matrix_)
}

# function which computes trace of a matrix
trace_mat <- function(matrix_) {
  return(sum(diag(matrix_)))
}

# function which regularizes a covariance matrix by adding a small
# positive value delta to the diagonal
regularize_covariance <- function(cov_mat_, alpha_ = 1e-2) {
  n <- nrow(cov_mat_)
  delta_ <- alpha_ * trace_mat(cov_mat_)
  cov_mat_ <- cov_mat_ + delta_ * diag(n)
  return(cov_mat_)
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
  # register parallel backend
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)

  # compute simulated phenotypes and distances
  df_results <- foreach(
    sim_num = 1:n_sim_abc,
    .export = c(
      "y", "x_mat", "z_mat", "k_mat", "beta_hat", "prior_sigma2_u", "prior_sigma2_e",
      "simulate_y", "squared_l2_norm", "simulate_and_compute_squared_l2_norm"
    ),
    .packages = c("MASS"),
    .combine = rbind
  ) %dopar% {
    set.seed(sim_num * seed_abc)
    simulate_and_compute_squared_l2_norm(
      y, x_mat, z_mat, k_mat, beta_hat, prior_sigma2_u, prior_sigma2_e
    )
  }
  # stop the parallel backend
  stopCluster(cl)
  registerDoSEQ()

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
    geno_df, raw_pheno_df, fixed_effects_vars, random_effect_vars, trait_,
    sigma2_u, sigma2_e, kernel_type,
    rate_decay_kernel,
    whitening_method) {
  # sigma2_u <- 1
  # sigma2_e <- 1
  # kernel_type <- "gaussian"
  # rate_decay_kernel <- 0.1
  # random_effect_vars <- 'Genotype'
  # fixed_effects_vars <- c("Envir", "Country",
  # "Year", "Row", "Position","Management")
  tryCatch(
    {
      # remove all rows with na w.r.t to trait_
      raw_pheno_df <- raw_pheno_df %>% drop_na(all_of(trait_))

      # define variables of interest
      sel_vars_ <- c(fixed_effects_vars, random_effect_vars, trait_)

      # get only variables of interest from raw_pheno_df
      raw_pheno_df <- raw_pheno_df[, sel_vars_]
      raw_pheno_df <- na.omit(raw_pheno_df)

      # compute Gram matrix (i.e. genomic covariance matrix)
      geno_names <- rownames(geno_df)
      geno_df <- apply(geno_df, 2, as.numeric)

      if (kernel_type == "linear") {
        k_mat <- tcrossprod(scale(apply(geno_df, 2, as.numeric),
          center = T, scale = F
        ))
      } else if (kernel_type == "gaussian") {
        kernel_function <- rbfdot(sigma = (1 / ncol(geno_df)) * rate_decay_kernel)
        k_mat <- kernelMatrix(kernel_function, geno_df)
      } else {
        # kernel identity is not recommended due to constrained hypothesis about
        # genotypes independence which may lead to low precision
        k_mat <- as.matrix(diag(nrow(geno_df)))
      }

      # test positive definiteness and force it if necessary
      if (!is.positive.definite(k_mat, tol = 1e-8)) {
        k_mat <- as.matrix(nearPD(k_mat)$mat)
      }

      # assign genotype rownames and colnames to k_mat
      colnames(k_mat) <- rownames(k_mat) <- geno_names

      # get common geontype between raw_pheno_df and geno_df
      raw_pheno_df <- raw_pheno_df[
        raw_pheno_df$Genotype %in% rownames(k_mat),
      ]

      # convert fixed effects variables to factors, and remove
      # buffer for management if exists
      for (fix_eff_var_ in fixed_effects_vars) {
        raw_pheno_df[, fix_eff_var_] <- as.factor(raw_pheno_df[, fix_eff_var_])
        if ("BUFFER" %in% raw_pheno_df[, fix_eff_var_]) {
          raw_pheno_df <- raw_pheno_df[
            raw_pheno_df[, fix_eff_var_] != "BUFFER",
          ]
        }
      }
      # droplevels in order to remove levels which don't exist anymore
      raw_pheno_df <- droplevels(raw_pheno_df)

      # get raw phenotypes associated to common genotypes
      y <- raw_pheno_df[, trait_]

      # get incidence matrices for fixed and random effects
      # NB. column of ones is added for intercept associated to fixed effects

      # define list of incidence matrices for fixed effects
      list_x_mat <- vector("list", length(fixed_effects_vars))
      names(list_x_mat) <- fixed_effects_vars

      # add incidence matrix for first fixed effect to list of matrices
      fix_eff_var_ <- fixed_effects_vars[1]
      list_x_mat[[fix_eff_var_]] <- model.matrix(
        as.formula(paste0("~", fix_eff_var_)),
        data = raw_pheno_df
      )
      colnames(list_x_mat[[fix_eff_var_]]) <- str_replace_all(
        colnames(list_x_mat[[fix_eff_var_]]),
        pattern = fix_eff_var_, replacement = paste0(fix_eff_var_, "_")
      )

      # add incidence matrices (without intercept) for other fixed effects to list
      for (fix_eff_var_ in fixed_effects_vars[-1]) {
        list_x_mat[[fix_eff_var_]] <- model.matrix(
          as.formula(paste0("~", fix_eff_var_, " - 1")),
          data = raw_pheno_df
        )
        colnames(list_x_mat[[fix_eff_var_]]) <- str_replace_all(
          colnames(list_x_mat[[fix_eff_var_]]),
          pattern = fix_eff_var_, replacement = paste0(fix_eff_var_, "_")
        )
      }
      x_mat <- do.call(cbind, list_x_mat)
      x_mat <- apply(x_mat, 2, as.numeric)

      # define list of incidence matrices for random effects
      list_z_mat <- vector("list", length(random_effect_vars))
      names(list_z_mat) <- random_effect_vars

      # add incidence matrices for random effects to list
      for (rand_eff_var in random_effect_vars) {
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

      # compute Σ (sig_mat_) and its Cholesky decomposition, i.e. Σ = LL'
      sig_mat_ <- sigma2_u * crossprod(t(z_mat), tcrossprod(k_mat, z_mat))
      sig_mat_ <- regularize_covariance(sig_mat_, alpha_ = 0.01)
      # regularize_covariance() adds α * trace(Σ) * I_n to the diagonal of Σ
      # to ensure its positive semi definiteness (PSD)

      if (whitening_method == "Cholesky") {
        # compute w_mat = L^−1 from Cholesky decomposition
        w_mat <- Matrix::solve(t(cholesky(sig_mat_, parallel = T)))
      } else {
        # compute w_mat from ZCA-cor
        w_mat <- whiteningMatrix(sig_mat_, method = "ZCA-cor")
      }

      # NB. intercept is already present in x_mat and x_mat_tilde
      x_mat_tilde <- w_mat %*% x_mat

      # get ols estimates for fixed effects and xi
      beta_hat <- ginv(t(x_mat_tilde) %*% x_mat_tilde) %*% t(x_mat_tilde) %*% y
      xi_hat <- y - x_mat_tilde %*% beta_hat

      return(list(
        "x_mat" = x_mat,
        "z_mat" = z_mat,
        "k_mat" = k_mat,
        "beta_hat" = beta_hat,
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
estimate_wiser_phenotype <- function(geno_df, raw_pheno_df, trait_,
                                    fixed_effects_vars = c(
                                      "Envir", "Country", "Year",
                                      "Row", "Position", "Management"
                                    ),
                                    random_effect_vars = "Genotype",
                                    init_sigma2_u = 1,
                                    init_sigma2_e = 1,
                                    n_sim_abc = 100,
                                    seed_abc = 123,
                                    quantile_threshold_abc = 0.05,
                                    nb_iter_abc = 1,
                                    kernel_type = "linear",
                                    rate_decay_kernel = 0.1,
                                    whitening_method = "Cholesky") {
  # init_sigma2_u = 1
  # init_sigma2_e = 1
  # n_sim_abc = 100
  # seed_abc = 123
  # quantile_threshold_abc = 0.05
  # nb_iter_abc = 1
  # kernel_type = "gaussian"
  # rate_decay_kernel = 0.1
  tryCatch(
    {
      # compute transformed variables associated to fixed effects and least-squares
      # to estimate these
      transform_and_ls_obj <- compute_transformed_vars_and_ols_estimates(
        geno_df, raw_pheno_df, fixed_effects_vars, random_effect_vars, trait_,
        sigma2_u = init_sigma2_u,
        sigma2_e = init_sigma2_e,
        kernel_type, rate_decay_kernel,
        whitening_method
      )

      # get an upper bound for sigma2_u et sigma2_e priors
      prior_sigma2_upper_bound <- var(
        transform_and_ls_obj$y
      )

      for (iter_ in 1:nb_iter_abc) {
        # print(paste0('iter : ', iter_))
        # compute variance components using abc
        var_comp_abc_obj <- abc_variance_component_estimation(
          y = transform_and_ls_obj$y,
          x_mat = transform_and_ls_obj$x_mat,
          z_mat = transform_and_ls_obj$z_mat,
          k_mat = transform_and_ls_obj$k_mat,
          beta_hat = transform_and_ls_obj$beta_hat,
          prior_sigma2_u = c(1e-2, prior_sigma2_upper_bound),
          prior_sigma2_e = c(1e-2, prior_sigma2_upper_bound),
          n_sim_abc, seed_abc,
          quantile_threshold_abc
        )
        # compute variance components again with abc using new estimates
        transform_and_ls_obj <- compute_transformed_vars_and_ols_estimates(
          geno_df, raw_pheno_df, fixed_effects_vars, random_effect_vars, trait_,
          sigma2_u = var_comp_abc_obj$sigma2_u_hat_mean,
          sigma2_e = var_comp_abc_obj$sigma2_e_hat_mean,
          kernel_type, rate_decay_kernel,
          whitening_method
        )
      }

      # get estimated components after abc

      # extract estimated fixed effects
      beta_hat <- transform_and_ls_obj$beta_hat

      # compute phenotypic values using ols
      v_hat <- Matrix::solve(t(transform_and_ls_obj$z_mat) %*% transform_and_ls_obj$z_mat) %*%
        t(transform_and_ls_obj$z_mat) %*% transform_and_ls_obj$xi_hat

      return(list(
        "var_comp_abc_obj" = var_comp_abc_obj,
        "beta_hat" = beta_hat,
        "v_hat" = v_hat
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}
