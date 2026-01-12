# Title: MGM–Guided Multiple Imputation by Chained Equations (MICE) Simulation

# 1) Simulation setup -----------------------------------
library(MASS)
library(mgm)
library(mice)
library(dplyr)
library(caret)
library(qgraph)

all_fpr_tpr_data <- data.frame()
seeds <- 1:10

# 2) Simulation runs ---------------------------------
for (seednum in seeds) {
  # 2-1) Data generation and Type Specification-------
  mgmsampler_user <- function (factors, interactions, thresholds, sds, type, level,
                               N, nIter = 250, pbar = TRUE, divWarning = 10^3, returnChains = FALSE) {

    p <- length(type)
    n_order <- length(lapply(interactions, function(x) 1))

    if (any(level[type != "c"] != 1))
      stop("The levels of all non-categorical variables must be specified as c(1).")

    for (ord in 1:n_order) {
      rows_ord <- nrow(factors[[ord]])
      if (!is.null(rows_ord)) {
        for (row in 1:rows_ord) {
          if (!all(level[factors[[ord]][row, ]] == dim(interactions[[ord]][[row]]))) {
            stop(paste0("Incorrect dimensions for interaction ",
                        paste0(factors[[ord]][row, ], collapse = "-"),
                        ". Please correct the 'interactions' argument."))
          }
        }
      }
    }

    for (i in 1:n_order) {
      n_ints <- length(interactions[[i]])
      if (n_ints > 0) {
        for (row in 1:n_ints) {
          if (!inherits(interactions[[i]][[row]], "array") &&
              !inherits(interactions[[i]][[row]], "matrix")) {
            stop("Parameters of each interaction must be provided as a k-dimensional array.")
          }
        }
      }
    }

    for (i in 1:n_order) if (!is.null(factors[[i]]))
      if (ncol(factors[[i]]) != i + 1)
        stop(paste0("Matrix specifying ", i + 1, "-way interactions has incorrect dimensions."))

    for (i in 1:n_order) {
      n_ints <- nrow(factors[[i]])
      if (!is.null(n_ints)) {
        for (row in 1:n_ints) {
          check <- all.equal(level[factors[[i]][row, ]], dim(interactions[[i]][[row]]))
          if (any(check == FALSE))
            stop(paste0("Incorrect parameter array dimensions for interaction: ",
                        paste(level[factors[[i]][row, ]], collapse = " ")))
        }
      }
    }

    FlagSymmetric <- function(x) {
      vec_sim <- rep(NA, nrow(x))
      ind_ord <- ncol(x)
      counter <- 1
      for (i in 1:nrow(x)) {
        if (is.na(vec_sim[i])) {
          vec_sim[i] <- counter
          for (j in (i + 1):nrow(x)) {
            if ((i + 1) > nrow(x)) next
            ind <- x[j, ] %in% x[i, ]
            if (sum(ind) == ind_ord) vec_sim[j] <- counter
          }
          counter <- counter + 1
        }
      }
      return(vec_sim)
    }

    for (i in 1:n_order) {
      if (!is.null(factors[[i]])) {
        symflag <- FlagSymmetric(factors[[i]])
        dup <- duplicated(symflag)
        if (any(dup))
          stop(paste0("Interaction ",
                      paste(factors[[i]][min(which(dup)), ], collapse = " "),
                      " is specified twice."))
      }
    }

    if (!("g" %in% type)) {
      if (missing(sds)) sds <- NULL
    } else {
      if (any(is.na(sds[type == "g"]))) stop("Missing values in standard deviations.")
      if (any(sds[type == "g"] < 0)) stop("Standard deviations must be positive.")
      if (any(!is.finite(sds[type == "g"]))) stop("Standard deviations must be finite.")
    }

    if (length(type) != length(level))
      stop("Length of 'type' does not match length of 'level'.")
    if (length(thresholds) != length(type))
      stop("Length of 'thresholds' list does not match length of 'type'.")
    if (!all.equal(unlist(lapply(thresholds, length)), level))
      stop("Specified thresholds do not match the number of categories in 'level'.")

    for (i in 1:n_order) {
      if (!is.null(factors[[i]])) {
        factors[[i]] <- as.matrix(factors[[i]], ncol = ncol(factors[[i]]), nrow = nrow(factors[[i]]))
      }
    }

    mgmsamp_obj <- list(call = NULL, data = NULL, chains = NULL)
    mgmsamp_obj$call <- list(factors = factors, interactions = interactions,
                             thresholds = thresholds, sds = sds, type = type, N = N,
                             level = level, nIter = nIter, pbar = pbar, divWarning = divWarning)

    if (pbar == TRUE)
      pb <- txtProgressBar(min = 0, max = N, initial = 0, char = "-", style = 3)
    if (returnChains)
      chains <- array(NA, dim = c(N, nIter, p))

    data <- matrix(NA, nrow = N, ncol = p)

    for (case in 1:N) {
      sampling_seq <- matrix(NA, nrow = nIter, ncol = p)

      for (v in 1:p) {
        if (type[v] == "g") sampling_seq[1, v] <- rnorm(1)
        if (type[v] == "p") sampling_seq[1, v] <- rpois(1, lambda = 0.5)
        if (type[v] == "c") sampling_seq[1, v] <- sample(1:level[v], size = 1)
      }

      for (iter in 2:nIter) {
        for (v in 1:p) {
          if (type[v] != "c") {
            l_natPar <- list()
            l_natPar[[1]] <- thresholds[[v]]

            for (ord in 1:n_order) {
              n_rows <- nrow(factors[[ord]])
              l_row_terms <- list()
              row_counter <- 1
              if (!is.null(n_rows)) {
                for (row in 1:n_rows) {
                  f_ind <- factors[[ord]][row, ]
                  if (v %in% f_ind) {
                    fill_in <- rep(NA, ord + 1)
                    get_cont <- rep(NA, ord + 1)
                    k_counter <- 1
                    for (k in f_ind) {
                      gibbs_shift <- if (k >= v) (iter - 1) else iter
                      if (level[k] == 1) fill_in[k_counter] <- 1 else fill_in[k_counter] <- sampling_seq[gibbs_shift, k]
                      if (level[k] == 1) get_cont[k_counter] <- sampling_seq[gibbs_shift, k] else get_cont[k_counter] <- 1
                      k_counter <- k_counter + 1
                    }
                    if (any(abs(sampling_seq[!is.na(sampling_seq)]) > divWarning))
                      warning("Gibbs sampler diverged; adapt parameterization of continuous variables.")
                    m_fill_in <- matrix(fill_in, ncol = length(fill_in))
                    row_int <- interactions[[ord]][[row]][m_fill_in]
                    l_row_terms[[row_counter]] <- prod(get_cont[-which(f_ind == v)]) * row_int
                    row_counter <- row_counter + 1
                  }
                }
              }
              l_natPar[[ord + 1]] <- unlist(l_row_terms)
            }

            natPar <- sum(unlist(l_natPar))

            if (type[v] == "g")
              sampling_seq[iter, v] <- rnorm(1, mean = natPar, sd = sds[v])
            if (type[v] == "p") {
              if (natPar <= 0)
                stop(paste0("Lambda parameter of Poisson variable ", v, " is nonpositive. Sampler returns NA."))
              sampling_seq[iter, v] <- rpois(1, lambda = natPar)
            }
          }

          if (type[v] == "c") {
            v_Potential <- rep(NA, level[v])
            for (cat in 1:level[v]) {
              l_natPar <- list()
              l_natPar[[1]] <- thresholds[[v]][cat]
              for (ord in 1:n_order) {
                n_rows <- nrow(factors[[ord]])
                l_row_terms <- list()
                row_counter <- 1
                if (!is.null(n_rows)) {
                  for (row in 1:n_rows) {
                    f_ind <- factors[[ord]][row, ]
                    if (v %in% f_ind) {
                      fill_in <- rep(NA, ord + 1)
                      get_cont <- rep(NA, ord + 1)
                      k_counter <- 1
                      for (k in f_ind) {
                        gibbs_shift <- if (k >= v) (iter - 1) else iter
                        if (level[k] == 1) fill_in[k_counter] <- 1 else fill_in[k_counter] <- sampling_seq[gibbs_shift, k]
                        if (level[k] == 1) get_cont[k_counter] <- sampling_seq[gibbs_shift, k] else get_cont[k_counter] <- 1
                        k_counter <- k_counter + 1
                      }
                      if (any(abs(sampling_seq[!is.na(sampling_seq)]) > divWarning))
                        warning("Gibbs sampler diverged; adapt parameterization of continuous variables.")
                      fill_in[which(v == f_ind)] <- cat
                      m_fill_in <- matrix(fill_in, ncol = length(fill_in))
                      row_int <- interactions[[ord]][[row]][m_fill_in]
                      l_row_terms[[row_counter]] <- prod(get_cont[-which(f_ind == v)]) * row_int
                      row_counter <- row_counter + 1
                    }
                  }
                }
                l_natPar[[ord + 1]] <- unlist(l_row_terms)
              }
              v_Potential[cat] <- sum(unlist(l_natPar))
            }
            v_probabilities <- exp(v_Potential) / sum(exp(v_Potential))
            sampling_seq[iter, v] <- sample(1:level[v], size = 1, prob = v_probabilities)
          }
        }
      }

      if (returnChains) chains[case, , ] <- sampling_seq
      data[case, ] <- sampling_seq[nIter, ]

      if (pbar == TRUE) setTxtProgressBar(pb, case)
    }

    mgmsamp_obj$data <- data
    if (returnChains) mgmsamp_obj$chains <- chains
    return(mgmsamp_obj)
  }

  generate_mixed_data <- function(adj, nsample, seednum) {
    # This function generates mixed data with 
    # p/2 continuous types and p/2 binary types.
    #
    # adj    : The dependency structure is copied from "adj".
    # nsample: The number of samples to generate
    # seednum: To specify seeds

    p <- nrow(adj)

    idx_pair <- which(adj != 0, arr.ind = TRUE)
    idx_pair <- idx_pair[idx_pair[, 1] < idx_pair[, 2], ]

    type <- rep(c("g", "c"), each = p / 2)
    level <- rep(c(1, 2), each = p / 2)

    # 1. Specify Model
    # (a) Edge set
    factors <- list()
    factors[[1]] <- as.array(idx_pair)

    # (b) Edge weights by type combination
    interactions <- list()
    interactions[[1]] <- list()
    for (i in 1:nrow(factors[[1]])) {
      type_pair <- paste0(sort(type[factors[[1]][i, ]]), collapse = "-")
      switch(
        type_pair,
        "g-g" = {
          interactions[[1]][[i]] <-
            array(0.5, dim = c(level[factors[[1]][i, 1]], level[factors[[1]][i, 2]]))
        },
        "c-g" = {
          interactions[[1]][[i]] <-
            array(c(1, 0), dim = c(level[factors[[1]][i, 1]], level[factors[[1]][i, 2]]))
        },
        "c-c" = {
          interactions[[1]][[i]] <-
            array(c(1, 1, 0, 0),
                  dim = c(level[factors[[1]][i, 1]], level[factors[[1]][i, 2]]))
        }
      )
    }

    # (c) Define Thresholds
    thresholds <- list()
    for (i in seq_along(type)) {
      switch(
        type[i],
        "g" = { thresholds[[i]] <- 0 },
        "c" = { thresholds[[i]] <- rep(0, level[i]) }
      )
    }

    # (d) Define Variances
    sds <- rep(.1, length(type))

    # 2. Sample cases
    set.seed(seednum)
    data <- mgmsampler_user(factors = factors,
                            interactions = interactions,
                            thresholds = thresholds,
                            sds = sds,
                            type = type,
                            level = level,
                            N = nsample,
                            nIter = 100,
                            pbar = FALSE)

    return(list(data = data))
  }

  p <- 20
  Omega <- diag(1, p)
  for (i in 1:(p - 1)) {
    Omega[i, i + 1] <- 0.5
    Omega[i + 1, i] <- 0.5
  }

  # Generate mixed data
  res_mgmsample <- generate_mixed_data(Omega, nsample = 100, seednum = seednum * 100)

  # Convert selected columns to categorical (factors)
  convert_to_categorical <- function(data, categorical_cols) {
    if (is.matrix(data)) data <- as.data.frame(data)
    for (col in categorical_cols) {
      data[[col]] <- as.factor(data[[col]])
    }
    return(data)
  }
  completed_data <- convert_to_categorical(res_mgmsample$data$data, 11:20)
  head(completed_data)

  # Check column classes
  check_column_types <- function(data) {
    column_types <- sapply(data, class)
    return(column_types)
  }
  check_column_types(completed_data)

  # 2-2) Inject MCAR missingness at 10% ---------
  set.seed(seednum * 100)
      
  create_MCAR <- function(data, perc = 0.1) {
    data_na <- data
    total_cells <- nrow(data_na) * ncol(data_na)
    num_na <- round(perc * total_cells)
    na_indices <- sample(seq_len(total_cells), size = num_na)
    for (index in na_indices) {
      row <- ceiling(index / ncol(data_na))
      col <- (index - 1) %% ncol(data_na) + 1
      data_na[row, col] <- NA
    }
    return(data_na)
  }
  data_na <- create_MCAR(completed_data)
  head(data_na)

  # 2-3) Initial mean imputation --------
  initial_imputation <- function(data) {
    data_imputed <- data
    for (col in 1:ncol(data_imputed)) {
      if (is.numeric(data_imputed[, col])) {
        mean_value <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[is.na(data_imputed[, col]), col] <- mean_value
      } else {
        mode_value <- names(sort(table(data_imputed[, col]), decreasing = TRUE))[1]
        data_imputed[is.na(data_imputed[, col]), col] <- mode_value
      }
    }
    return(data_imputed)
  }
  initial_imputed_data <- initial_imputation(data_na)
  head(initial_imputed_data)

  # 2-4) MICE-by-mgm iterative imputation--------
  data <- initial_imputed_data  # working copy

  # Automatically determine mgm 'type' and 'level' from data frame
  determine_type_and_level <- function(data) {
    types <- sapply(data, function(col) {
      if (is.numeric(col)) {
        return("g")
      } else if (is.factor(col)) {
        return("c")
      } else {
        stop("Unsupported data type encountered.")
      }
    })
    levels <- sapply(data, function(col) {
      if (is.factor(col)) {
        return(length(levels(col)))
      } else {
        return(1)
      }
    })
    return(list(type = types, level = levels))
  }

  # Convert factors to numeric for mgm input
  convert_factors_to_numeric <- function(data) {
    data_numeric <- data
    for (col in 1:ncol(data_numeric)) {
      if (is.factor(data_numeric[, col])) {
        data_numeric[, col] <- as.numeric(data_numeric[, col])
      }
    }
    return(data_numeric)
  }

  # Configure MICE methods by column type
  mice_methods <- make.method(data)
  for (col in names(data)) {
    if (is.factor(data_na[[col]])) {
      mice_methods[col] <- "polyreg"  # categorical
    } else {
      mice_methods[col] <- "pmm"      # numeric
    }
  }

  # Missingness pattern (where to impute)
  missing_pattern <- is.na(data_na)

  # Generate 5 pseudo-completed datasets using the proposed iterative scheme
  for (data_set_number in 1:5) {
    for (iter in 1:5) {
      # Determine types/levels
      type_level <- determine_type_and_level(data)
      types <- type_level$type
      levels <- type_level$level

      # Numeric conversion for mgm
      numeric_imputed_data <- convert_factors_to_numeric(data)
      imputed_matrix <- as.matrix(numeric_imputed_data)

      # Fit mgm model
      sink(tempfile())
      mgm_model <- mgm(
        data = imputed_matrix, type = types, level = levels,
        lambdaSel = "CV", lambdaFolds = 10, ruleReg = "AND"
      )
      sink()
      adj_matrix <- mgm_model$pairwise$wadj

      # Build predictor matrix from adjacency
      non_zero_elements <- abs(adj_matrix) > 0
      pred_matrix <- ifelse(non_zero_elements, 1, 0)
      diag(pred_matrix) <- 0

      # One-step MICE using that graph as predictors
      mice_imputation <- mice(
        data,
        predictorMatrix = pred_matrix,
        method = mice_methods,
        where = missing_pattern,
        m = 1,
        maxit = 1,
        printFlag = FALSE
      )
      completed_data <- complete(mice_imputation)

      # Update data for next iteration
      data <- completed_data
    }
    # Save each pseudo-completed dataset to global env (to match original behavior)
    assign(paste("pseudo_completed", data_set_number, sep = ""), data, envir = .GlobalEnv)
  }

  # 2-5) Final adjacency matrices ---------

  final_adj_matrix <- function(data_list) {
    adj_matrices <- list()
    for (i in 1:length(data_list)) {
      data <- data_list[[i]]
      type_level <- determine_type_and_level(data)
      types <- type_level$type
      levels <- type_level$level

      numeric_imputed_data <- convert_factors_to_numeric(data)
      imputed_matrix <- as.matrix(numeric_imputed_data)

      sink(tempfile())
      mgm_model <- mgm(
        data = imputed_matrix, type = types, level = levels,
        lambdaSel = "CV", lambdaFolds = 10, ruleReg = "AND"
      )
      sink()
      adj_matrix <- mgm_model$pairwise$wadj

      # binarize adjacency
      binary_adj_matrix <- ifelse(adj_matrix > 0, 1, 0)
      adj_matrices[[i]] <- binary_adj_matrix
    }
    return(list(adj_matrices = adj_matrices))
  }

  # Majority vote across a list of binary adjacency matrices
  calculate_majority_vote <- function(results) {
    sum_matrix <- Reduce("+", results)
    majority_vote_matrix <- ifelse(sum_matrix > length(results) / 2, 1, 0)
    return(majority_vote_matrix)
  }

  # (i) Completed-data (single run)
  result_completed_data <- final_adj_matrix(list(completed_data))
  adj_matrix_completed_data <- result_completed_data$adj_matrices[[1]]

  # (ii) Proposed method (majority over 5 pseudo-completed sets)
  result_proposed <- final_adj_matrix(
    list(pseudo_completed1, pseudo_completed2, pseudo_completed3, pseudo_completed4, pseudo_completed5)
  )
  adj_matrix_proposed <- calculate_majority_vote(result_proposed$adj_matrices)

  # (iii) Plain MICE (m=5, majority vote)
  mice_data <- mice(data_na, m = 5, printFlag = FALSE)
  for (i in 1:5) {
    assign(paste0("mice_data", i), complete(mice_data, action = i))
  }
  result_mice <- final_adj_matrix(list(mice_data1, mice_data2, mice_data3, mice_data4, mice_data5))
  adj_matrix_mice <- calculate_majority_vote(result_mice$adj_matrices)

  # (iv) Ground-truth binary adjacency from Omega
  Omega_bin <- ifelse(abs(Omega) > 0, 1, 0)
  diag(Omega_bin) <- 0

  # 2-6) Evaluation of FPR/TPR across Omega values --------

  calculate_fpr_tpr_points <- function(predicted_adj, true_adj) {
    predicted <- as.vector(predicted_adj)
    true <- as.vector(true_adj)
    cm <- confusionMatrix(factor(predicted, levels = c(0, 1)),
                          factor(true, levels = c(0, 1)))
    tpr <- cm$byClass["Sensitivity"]
    fpr <- 1 - cm$byClass["Specificity"]
    return(data.frame(FPR = fpr, TPR = tpr, Model = "Model"))
  }

  fpr_tpr_proposed  <- calculate_fpr_tpr_points(adj_matrix_proposed,        Omega_bin); fpr_tpr_proposed$Model  <- "Proposed"
  fpr_tpr_completed <- calculate_fpr_tpr_points(adj_matrix_completed_data,  Omega_bin); fpr_tpr_completed$Model <- "Completed"
  fpr_tpr_mice      <- calculate_fpr_tpr_points(adj_matrix_mice,            Omega_bin); fpr_tpr_mice$Model      <- "Mice"

  fpr_tpr_data <- rbind(fpr_tpr_proposed, fpr_tpr_completed, fpr_tpr_mice)
  fpr_tpr_data$Seed <- seednum
  all_fpr_tpr_data <- rbind(all_fpr_tpr_data, fpr_tpr_data)
}

# 3) Final results ---------------------------

# 3-1) Print all FPR/TPR pairs
print(all_fpr_tpr_data)

# 3-2) Model-wise means
mean_fpr_tpr_data <- all_fpr_tpr_data %>%
  group_by(Model) %>%
  summarise(FPR = mean(FPR), TPR = mean(TPR)) %>%
  mutate(Seed = "Mean")
print(mean_fpr_tpr_data)

# 3-3) ROC scatter (requires ggplot2)
# If ggplot2 is not already attached in your R session, run: install.packages("ggplot2"); library(ggplot2)
ggplot(all_fpr_tpr_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_point(size = 3) +
  geom_text(aes(label = Seed), hjust = -1.5, vjust = 0.3, size = 3) +
  geom_point(data = mean_fpr_tpr_data, shape = 4, size = 3, stroke = 1.5) +
  geom_text(data = mean_fpr_tpr_data, aes(label = Seed), hjust = -0.5, vjust = 0.2, size = 3, color = "red") +
  labs(title = "ROC Plot by run",
       x = "False Positive Rate (FPR)",
       y = "True Positive Rate (TPR)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
