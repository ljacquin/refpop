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
        return(df_envir_)
      } else {
        return(data.frame())
      }
    },
    error = function(e) {
      return(data.frame())
      cat(
        "Error with spat_hetero_env_correct_trait,
      here is the possible issue with data and/or computation : ",
        conditionMessage(e), "\n"
      )
    }
  )
}


# function to re-order column names according to a specified order
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

# function to detect dates (does not work for all dates format, be careful)
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

# function which multi-location clonal mean heritability (WIP)
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

# (DEPRECATED) functions which corrects for trait spatial heterogeneity in raw_data
correct_trait_spatial_heterogeneity_raw_data <- function(
    trait_, geno_df_, df_year_, raw_data_site_names_, year_) {
  tryCatch(
    {
      for (site_name_ in raw_data_site_names_) {
        df_site_ <- df_year_ %>% subset(Country == site_name_)
        df_site_ <- df_site_[, -match("Country", colnames(df_site_))]
        df_site_ <- drop_na(df_site_)
        df_site_$R <- as.factor(df_site_$Row)
        df_site_$P <- as.factor(df_site_$Position)
        df_site_$Genotype <- as.factor(df_site_$Genotype)
        spats_formula <- as.formula(~ PSANOVA(Position, Row))
        spats_model <- SpATS(
          response = trait_,
          spatial = spats_formula,
          genotype = "Genotype", genotype.as.random = TRUE, fixed = NULL,
          random = ~ R + P, data = df_site_,
          control = controlSpATS()
        )
        spats_pred <- predict(spats_model, which = "Genotype")
        spats_pred <- spats_pred[, c("Genotype", "predicted.values")]
        colnames(spats_pred)[2] <- paste0(site_name_, "_", year_)
        spats_pred$Genotype <- as.character(spats_pred$Genotype)
        geno_df_ <- merge(geno_df_, spats_pred, by = "Genotype", all = TRUE)
      }
      return(geno_df_)
    },
    error = function(e) {
      cat(
        "Error with correct_trait_spatial_heterogeneity_raw_data,
      here is the possible issue with data and/or computation : ",
        conditionMessage(e), "\n"
      )
    }
  )
}

# (DEPRECATED) functions which corrects for trait spatial heterogeneity in raw_data
correct_trait_spatial_heterogeneity_munq_data <- function(
    trait_, geno_df_, list_df_year_, data_site_names_year_) {
  tryCatch(
    {
      for (site_name_ in data_site_names_year_) {
        print(site_name_)
        df_site_ <- list_df_year_[[site_name_]]
        df_site_[is.na(df_site_$MUNQ), "MUNQ"] <- df_site_[is.na(df_site_$MUNQ), "Genotype"]
        df_site_ <- df_site_[, -match("Genotype", colnames(df_site_))]
        df_site_ <- drop_na(df_site_)
        df_site_$R <- as.factor(df_site_$Row)
        df_site_$P <- as.factor(df_site_$Position)
        df_site_$MUNQ <- as.factor(df_site_$MUNQ)
        spats_formula <- as.formula(~ PSANOVA(Position, Row))
        spats_model <- SpATS(
          response = trait_,
          spatial = spats_formula,
          genotype = "MUNQ", genotype.as.random = TRUE, fixed = NULL,
          random = ~ R + P, data = df_site_,
          control = controlSpATS()
        )
        spats_pred <- predict(spats_model, which = "MUNQ")
        spats_pred <- spats_pred[, c("MUNQ", "predicted.values")]
        spats_pred$MUNQ <- as.character(spats_pred$MUNQ)
        colnames(spats_pred) <- c("Genotype", site_name_)
        geno_df_ <- merge(geno_df_, spats_pred, by = "Genotype", all = TRUE)
      }
      return(geno_df_)
    },
    error = function(e) {
      cat(
        "Error with correct_trait_spatial_heterogeneity_raw_data,
      here is the possible issue with data and/or computation : ",
        conditionMessage(e), "\n"
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
  mean_value <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mean_value
  return(x)
}

# function which computes the number of necessary components to reach at least percent_explained_variance_
n_comp_required_for_percent_explained_var <- function(facto_mine_pca_model_,
                                                      percent_explained_var) {
  return(as.numeric(which(facto_mine_pca_model_$eig[, 3] >
    percent_explained_var))[1])
}


# function which divide a vector of indices into a list ofn pieces
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

# function which test if a matrix is a squared one
is_square <- function(mat) {
  nrow(mat) == ncol(mat)
}

# function which recode values for biallelic markers codes as A and B
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

# function which recode values of a matrix for biallelic markers codes as A and B
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

# function which replace missing values of parents by genotype values
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

# function which creates pedigree incidence matrix
create_pedig_incid_mat_2 <- function(df_) {
  # get unique genotypes, parents, families and origin
  genotypes <- unique(df_$Genotype)
  parents <- unique(c(df_$P1, df_$P2))
  families <- unique(df_$Family)
  origins <- unique(df_$Origin)

  # initialize incidence data with an empty data frame
  incid_df_ <- data.frame(Genotype = genotypes)

  # add columns for parent, family and origin
  for (parent in parents) {
    incid_df_[[parent]] <- 0
  }
  for (family in families) {
    incid_df_[[family]] <- 0
  }
  for (origin in origins) {
    incid_df_[[origin]] <- 0
  }

  # fill the incidence matrix
  for (i in 1:nrow(df_)) {
    genotype <- df_$Genotype[i]
    parent_1 <- df_$P1[i]
    parent_2 <- df_$P2[i]
    family <- df_$Family[i]
    origin <- df_$Origin[i]

    # assign 1 to indicate relatedness, family membership and origin
    incid_df_[incid_df_$Genotype == genotype, parent_1] <- 1
    incid_df_[incid_df_$Genotype == genotype, parent_2] <- 1
    incid_df_[incid_df_$Genotype == genotype, family] <- 1
    incid_df_[incid_df_$Genotype == genotype, origin] <- 1
  }
  return(incid_df_)
}

# function to split values of a column
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

# function which tune mtry for ranger random foest
tune_mtry_ranger_rf <- function(X, Y, 
                                mtry_grid_,
                                num_trees_ = 500,
                                pkgs_to_export_
                                ) {
  # initialize the cluster
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  tryCatch(
    {
      vect_acc_ <- foreach(mtry_ = mtry_grid_, .combine = c,
                           .packages = pkgs_to_export_) %dopar% {
        # build model on train data
        rf_model <- ranger(
          y = Y,
          x = X,
          mtry = mtry_,
          num.trees = num_trees_
        )
        # correlate out-of-bag (OOB) predictions (asymptotically equivalent to 
        # LOOCV on large samples) with observed values
        acc_ <- cor(rf_model$predictions, Y)
        names(acc_) <- as.character(mtry_)
        acc_
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
    "vect_acc_" = vect_acc_,
    "opt_mtry" = as.numeric(names(which.max(vect_acc_)))
  ))
}

# function which tunes the epsilon hyperparameter for support vector regression
tune_eps_ksvm_reg <- function(X, Y, kpar_, type_, kernel_, c_par_,
                              epsilon_grid_, n_folds_) {
  expected_loss_grid_ <- rep(Inf, length(epsilon_grid_))
  n <- length(Y)
  Folds <- cvFolds(n, n_folds_, type = "consecutive")

  tryCatch(
    {
      l <- 1
      for (eps_ in epsilon_grid_)
      {
        vect_loss_folds <- rep(0, n_folds_)

        for (fold_ in 1:n_folds_)
        {
          idx_fold_ <- which(Folds$which == fold_)

          # valid set
          y_val_ <- Y[idx_fold_]
          x_val_ <- X[idx_fold_, ]

          # train set
          y_train_ <- Y[-idx_fold_]
          x_train_ <- X[-idx_fold_, ]

          # build ksvm model on train set
          ksvm_model <- ksvm(
            x = X_minus_fold_, y = Y_minus_fold_, kpar = kpar_,
            type = type_, kernel = kernel_,
            C = c_par_, epsilon = eps_
          )

          f_hat_val_ <- predict(ksvm_model, x_val_)

          # loss for fold_
          vect_loss_folds[fold_] <- sum((y_val_ - f_hat_val_)^2)
        }

        expected_loss_grid_[l] <- mean(vect_loss_folds)
        l <- l + 1
      }
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )

  # get the optimal epsilon
  optimal_eps_ <- epsilon_grid_[which.min(expected_loss_grid_)]

  # train the "optimal" model
  tuned_model <- ksvm(
    x = X, y = Y, kpar = kpar_, type = type_,
    kernel = kernel_, C = c_par_,
    epsilon = optimal_eps_
  )

  return(list(
    "tuned_ksvm" = tuned_model, "optimal_eps_" = optimal_eps_,
    "expected_loss_grid_" = expected_loss_grid_,
    "epsilon_grid_" = epsilon_grid_
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

# function which remove monomorphic markers
remove_monomorphic_markers <- function(geno_df) {
  # identify monomorphic markers
  monomorphic_markers <- apply(
    geno_df[, -match("Genotype", colnames(geno_df))],
    2, function(col) length(unique(col)) == 1
  )

  # filter the monomorphic markers
  geno_df_filtered <- geno_df[, !monomorphic_markers]

  # get the names of the monomorphic markers
  monomorphic_marker_names <- colnames(geno_df)[monomorphic_markers]

  # return the filtered data frame and the list of monomorphic markers
  return(list(
    "filtered_df" = geno_df_filtered,
    "monomorphic_markers" = monomorphic_marker_names
  ))
}


# function to remove columns with variance below a threshold
remove_low_variance_columns <- function(geno_df, threshold) {
  # calculate variance for each column
  variances <- apply(geno_df[, -match("Genotype", colnames(geno_df))], 2, var)

  # identify columns with variance below the threshold
  low_variance_columns <- variances < threshold

  # filter out low variance columns
  geno_df_filtered <- geno_df[, !low_variance_columns]

  # get the names of the low variance columns
  low_variance_column_names <- colnames(geno_df)[low_variance_columns]

  # return the filtered data frame and the list of low variance column names
  return(list(
    "filtered_df" = geno_df_filtered,
    "low_variance_columns" = low_variance_column_names
  ))
}

# function to calculate the weighted sum of genotype columns without weights
sum_weighted_genotypes <- function(train_genotypes, test_genotype) {
  # calculate the Hamming distance between the test genotype and every genotype in the training data
  hamming_distances <- apply(train_genotypes, 1, function(row) hamming.distance(test_genotype, row))

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
    # Define the start and end columns for the current chunk
    start_col <- (i - 1) * n + 1
    end_col <- i * n

    # Reduce the chunk by summing up columns
    reduced_matrix[, i] <- rowSums(genotype_matrix[, start_col:end_col])
  }

  # return the reduced genotype matrix
  return(reduced_matrix)
}

# # example usage
# train_genotypes <- matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, byrow = TRUE)
# test_genotype <- c(0, 1, 0)
#
# # calculate the weighted sum of genotype columns for the test genotype
# result <- sum_weighted_genotypes(train_genotypes, test_genotype)
# print(result)
#
# # example usage for reducing genotype matrix
# reduced_matrix <- reduce_genotype_matrix(train_genotypes, 2)
# print(reduced_matrix)
