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
  x[is.na(x)] <- (mean_value + rnorm(1, 0, 10))
  return(x)
}

# function which computes the number of necessary components to reach at least percent_explained_variance_
n_comp_required_for_percent_explained_var <- function(facto_mine_pca_model_,
                                                      percent_explained_var) {
  return(as.numeric(which(facto_mine_pca_model_$eig[, 3] >=
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

# a modified version of the mahalanobis distance function from the robustbase
# package
mahalanobis_dist_ <- function(x, center, cov, inverted = FALSE, ...) {
  x <- if (is.vector(x)) {
    matrix(x, ncol = length(x))
  } else {
    as.matrix(x)
  }
  if (!isFALSE(center)) {
    x <- sweep(x, 2L, center)
  }
  if (!inverted) {
    cov <- as.matrix(nearPD(cov)$mat)
  }
  cov <- solve(cov, ...)
  setNames(rowSums(x %*% cov * x), rownames(x))
}

# a modified version of the minimum covariance determinant (MCD) function
# from the robustbase package
covMcd_ <- function(x, cor = FALSE, raw.only = FALSE, alpha = control$alpha,
                    nsamp = control$nsamp, nmini = control$nmini, kmini = control$kmini,
                    scalefn = control$scalefn, maxcsteps = control$maxcsteps,
                    initHsets = NULL, save.hsets = FALSE, names = TRUE, seed = control$seed,
                    tolSolve = control$tolSolve, trace = control$trace, use.correction = control$use.correction,
                    wgtFUN = control$wgtFUN, control = rrcov.control()) {
  logdet.Lrg <- 50
  if (length(seed) > 0) {
    if (length(seed) < 3L || seed[1L] < 100L) {
      stop("invalid 'seed'. Must be compatible with .Random.seed !")
    }
    if (!is.null(seed.keep <- get0(".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    ))) {
      on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
    assign(".Random.seed", seed, envir = .GlobalEnv)
  }
  defCtrl <- if (missing(control)) {
    control
  } else {
    rrcov.control()
  }
  if (missing(wgtFUN)) {
    getDefCtrl("wgtFUN", defCtrl)
  }
  if (is.null(nmini)) {
    getDefCtrl("nmini", defCtrl)
  }
  if (is.numeric(nsamp) && nsamp <= 0) {
    stop("Invalid number of trials nsamp = ", nsamp, "!")
  }
  if (is.data.frame(x)) {
    x <- data.matrix(x, rownames.force = FALSE)
  } else if (!is.matrix(x)) {
    x <- matrix(x, length(x), 1, dimnames = if (names) {
      list(names(x), deparse(substitute(x)))
    })
  }
  if (!names) {
    dimnames(x) <- NULL
  }
  ok <- is.finite(x %*% rep.int(1, ncol(x)))
  x <- x[ok, , drop = FALSE]
  if (!length(dx <- dim(x))) {
    stop("All observations have missing values!")
  }
  n <- dx[1]
  p <- dx[2]
  if (names) {
    dimn <- dimnames(x)
  }
  h <- h.alpha.n(alpha, n, p)
  if (n <= p + 1) {
    stop(if (n <= p) {
      "n <= p -- you can't be serious!"
    } else {
      "n == p+1  is too small sample size for MCD"
    })
  }
  if (n < 2 * p) {
    warning("n < 2 * p, i.e., possibly too small sample size")
  }
  if (h > n) {
    stop("Sample size n  <  h(alpha; n,p) := size of \"good\" subsample")
  } else if (2 * h < n) {
    warning("subsample size\t h < n/2  may be too small")
  }
  if (is.character(wgtFUN)) {
    if (is.function(mkWfun <- .wgtFUN.covMcd[[wgtFUN]])) {
      wgtFUN <- mkWfun(p = p, n = n, control)
    }
  }
  if (!is.function(wgtFUN)) {
    stop(
      gettextf(
        "'wgtFUN' must be a function or one of the strings %s.",
        pasteK(paste0("\"", names(.wgtFUN.covMcd), "\""))
      ),
      domain = NA
    )
  }
  raw.cnp2 <- cnp2 <- c(1, 1)
  ans <- list(call = match.call(), nsamp = nsamp, method = sprintf(
    "MCD(alpha=%g ==> h=%d)",
    alpha, h
  ))
  if (h == n) {
    mcd <- cov(x)
    loc <- as.vector(colMeans(x))
    obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
    if (-obj / p > logdet.Lrg) {
      ans$cov <- mcd
      if (names) {
        dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
      }
      if (cor) {
        ans$cor <- cov2cor(ans$cov)
      }
      ans$center <- loc
      if (names && length(dimn[[2]])) {
        names(ans$center) <- dimn[[2]]
      }
      ans$n.obs <- n
      ans$singularity <- list(kind = "classical")
      weights <- 1
    } else {
      mah <- mahalanobis_dist_(x, loc, mcd, tol = tolSolve)
      weights <- wgtFUN(mah)
      sum.w <- sum(weights)
      ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
      if (sum.w != n) {
        cnp2[1] <- .MCDcons(p, sum.w / n)
        ans$cov <- ans$cov * cnp2[1]
      }
      obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
      if (-obj / p > logdet.Lrg) {
        ans$singularity <- list(kind = "reweighted.MCD")
      } else {
        mah <- mahalanobis_dist_(x, ans$center, ans$cov, tol = tolSolve)
        weights <- wgtFUN(mah)
      }
    }
    ans$alpha <- alpha
    ans$quan <- h
    ans$raw.cov <- mcd
    ans$raw.center <- loc
    if (names && !is.null(nms <- dimn[[2]])) {
      names(ans$raw.center) <- nms
      dimnames(ans$raw.cov) <- list(nms, nms)
    }
    ans$crit <- obj
    ans$method <- paste(
      ans$method, "\nalpha = 1: The minimum covariance determinant estimates based on",
      n, "observations \nare equal to the classical estimates."
    )
    ans$mcd.wt <- rep.int(NA, length(ok))
    ans$mcd.wt[ok] <- weights
    if (names && length(dimn[[1]])) {
      names(ans$mcd.wt) <- dimn[[1]]
    }
    ans$wt <- NULL
    ans$X <- x
    if (names) {
      if (length(dimn[[1]])) {
        dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
      } else {
        dimnames(ans$X) <- list(
          seq(along = ok)[ok],
          NULL
        )
      }
    }
    if (trace) {
      cat(ans$method, "\n")
    }
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "mcd"
    return(ans)
  }
  mcd <- if (nsamp == "deterministic") {
    ans$method <- paste("Deterministic", ans$method)
    .detmcd(x, h,
      hsets.init = initHsets, save.hsets = save.hsets,
      scalefn = scalefn, maxcsteps = maxcsteps, trace = as.integer(trace),
      names = names
    )
  } else {
    ans$method <- paste0(
      "Fast ", ans$method, "; nsamp = ",
      nsamp, "; (n,k)mini = (", nmini, ",", kmini, ")"
    )
    .fastmcd(x, h, nsamp, nmini, kmini, trace = as.integer(trace))
  }
  calpha <- .MCDcons(p, h / n)
  correct <- if (use.correction) {
    .MCDcnp2(p, n, alpha)
  } else {
    1
  }
  raw.cnp2 <- c(calpha, correct)
  if (p == 1) {
    ans$method <- paste("Univariate", ans$method)
    scale <- sqrt(calpha * correct) * as.double(mcd$initcovariance)
    center <- as.double(mcd$initmean)
    if (abs(scale - 0) < 1e-07) {
      ans$singularity <- list(kind = "identicalObs", q = h)
      ans$raw.cov <- ans$cov <- matrix(0)
      ans$raw.center <- ans$center <- center
      ans$n.obs <- n
      ans$alpha <- alpha
      ans$quan <- h
      if (names && !is.null(nms <- dimn[[2]][1])) {
        names(ans$raw.center) <- names(ans$center) <- nms
        dimnames(ans$raw.cov) <- dimnames(ans$cov) <- list(
          nms,
          nms
        )
      }
      ans$crit <- -Inf
      weights <- as.numeric(abs(x - center) < 1e-07)
    } else {
      weights <- wgtFUN(((x - center) / scale)^2)
      sum.w <- sum(weights)
      ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
      if (sum.w != n) {
        cdelta.rew <- .MCDcons(p, 0.975)
        correct.rew <- if (use.correction) {
          .MCDcnp2.rew(p, n, alpha)
        } else {
          1
        }
        cnp2 <- c(cdelta.rew, correct.rew)
        ans$cov <- cdelta.rew * correct.rew * ans$cov
      }
      ans$alpha <- alpha
      ans$quan <- h
      ans$raw.cov <- as.matrix(scale^2)
      ans$raw.center <- as.vector(center)
      if (names && !is.null(nms <- dimn[[2]][1])) {
        dimnames(ans$raw.cov) <- list(nms, nms)
        names(ans$raw.center) <- nms
      }
      ans$crit <- log(sum(sort((x - as.double(mcd$initmean))^2,
        partial = h
      )[1:h]) / max(1, h - 1))
      center <- ans$center
      scale <- as.vector(sqrt(ans$cov))
      weights <- wgtFUN(((x - center) / scale)^2)
    }
  } else {
    mcd$initcovariance <- matrix(
      calpha * correct * mcd$initcovariance,
      p, p
    )
    if (raw.only || mcd$exactfit != 0) {
      if (!is.null(mcd$coeff)) {
        dim(mcd$coeff) <- c(5, p)
      }
      ans$cov <- ans$raw.cov <- mcd$initcovariance
      ans$center <- ans$raw.center <- as.vector(mcd$initmean)
      if (names && !is.null(nms <- dimn[[2]])) {
        dimnames(ans$cov) <- list(nms, nms)
        names(ans$center) <- nms
      }
      ans$n.obs <- n
      if (mcd$exactfit != 0) {
        if (!(mcd$exactfit %in% c(1, 2, 3))) {
          stop(
            "Unexpected 'exactfit' code ", mcd$exactfit,
            ". Please report!"
          )
        }
        ans$singularity <- list(
          kind = "on.hyperplane",
          exactCode = mcd$exactfit, p = p, h = h, count = mcd$kount,
          coeff = mcd$coeff[1, ]
        )
        ans$crit <- -Inf
        weights <- mcd$weights
      } else {
        ans$raw.only <- TRUE
        ans$crit <- mcd$mcdestimate
        weights <- mcd$weights
        if (is.null(mcd$weights)) {
          mah <- mahalanobis_dist_(x, mcd$initmean, mcd$initcovariance,
            tol = tolSolve
          )
          weights <- wgtFUN(mah)
        }
      }
      ans$alpha <- alpha
      ans$quan <- h
      if (names && !is.null(nms <- dimn[[2]])) {
        names(ans$raw.center) <- nms
        dimnames(ans$raw.cov) <- list(nms, nms)
      }
    } else {
      mah <- mahalanobis_dist_(x, mcd$initmean, mcd$initcovariance,
        tol = tolSolve
      )
      weights <- wgtFUN(mah)
      sum.w <- sum(weights)
      ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
      sing.rewt <- any(apply(ans$cov == 0, 2, all))
      if (!sing.rewt && sum.w != n) {
        cdelta.rew <- .MCDcons(p, 0.975)
        correct.rew <- if (use.correction) {
          .MCDcnp2.rew(p, n, alpha)
        } else {
          1
        }
        cnp2 <- c(cdelta.rew, correct.rew)
        ans$cov <- cdelta.rew * correct.rew * ans$cov
      }
      ans$best <- sort(as.vector(mcd$best))
      ans$alpha <- alpha
      ans$quan <- h
      ans$raw.cov <- mcd$initcovariance
      ans$raw.center <- as.vector(mcd$initmean)
      if (names && !is.null(nms <- dimn[[2]])) {
        names(ans$raw.center) <- nms
        dimnames(ans$raw.cov) <- list(nms, nms)
      }
      ans$raw.weights <- weights
      ans$crit <- mcd$mcdestimate
      ans$raw.mah <- mah
      if (sing.rewt || -determinant(ans$cov, logarithm = TRUE)$modulus[1] / p >
        logdet.Lrg) {
        ans$singularity <- list(kind = paste0(
          "reweighted.MCD",
          if (sing.rewt) "(zero col.)"
        ))
        ans$mah <- mah
      } else {
        mah <- mahalanobis_dist_(x, ans$center, ans$cov, tol = tolSolve)
        ans$mah <- mah
        weights <- wgtFUN(mah)
      }
    }
  }
  ans$mcd.wt <- rep.int(NA, length(ok))
  ans$mcd.wt[ok] <- weights
  if (names) {
    if (length(dimn[[1]])) {
      names(ans$mcd.wt) <- dimn[[1]]
    }
    if (length(dimn[[1]])) {
      dimnames(x)[[1]] <- names(ans$mcd.wt)[ok]
    } else {
      dimnames(x) <- list(seq(along = ok)[ok], NULL)
    }
  }
  ans$X <- x
  ans$wt <- NULL
  if (trace) {
    cat(ans$method, "\n")
  }
  ans$raw.cnp2 <- raw.cnp2
  ans$cnp2 <- cnp2
  if (nsamp == "deterministic") {
    ans <- c(ans, mcd[c("iBest", "n.csteps", if (save.hsets) "initHsets")])
  }
  class(ans) <- "mcd"
  if (is.list(ans$singularity)) {
    warning(paste(strwrap(.MCDsingularityMsg(
      ans$singularity,
      ans$n.obs
    )), collapse = "\n"), domain = NA)
  }
  ans
}
