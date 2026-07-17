###############################################################################
# Full workflow for UAV-based wheat yield prediction
#
# Input file format:
#   - One Excel workbook with six sheets.
#   - The six sheets represent:
#       1. Jointing
#       2. Booting
#       3. Heading
#       4. Flowering
#       5. Grain filling
#       6. Maturity
#   - In each sheet:
#       column 1: grain yield
#       columns 2 to n-1: vegetation indices
#       last column: NRCT
#
# Main workflow:
#   1. Fit each vegetation index or NRCT to grain yield by linear regression.
#   2. Use the fitted yield values to perform grey relational analysis (GRA).
#   3. Add predictors sequentially according to GRA ranking and select the
#      subset with the minimum training error.
#   4. Compare VI-only and VI+NRCT single-stage ENR models using 50 repeated
#      10-fold cross-validation.
#   5. Compare multistage feature-level fusion and entropy weight fusion (EWF)
#      using the same 50 repeated 10-fold cross-validation.
###############################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

packages <- c("readxl", "glmnet", "caret", "writexl")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(readxl)
library(glmnet)
library(caret)
library(writexl)

###############################################################################
# User settings
###############################################################################

# Replace this example path with the path used on your own computer.
input_file <- "E:/Wheat_Yield_Project/data/Wheat_VI_NRCT_by_growth_stage.xlsx"
output_dir <- "E:/Wheat_Yield_Project/results/full_workflow_outputs"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stage_names <- paste0("T", 1:6)
stage_labels <- c(
  T1 = "Jointing",
  T2 = "Booting",
  T3 = "Heading",
  T4 = "Flowering",
  T5 = "Grain filling",
  T6 = "Maturity"
)

fusion_combos <- list(
  C1 = c("T1", "T2", "T3"),
  C2 = c("T1", "T2", "T3", "T4"),
  C3 = c("T1", "T2", "T3", "T4", "T5"),
  C4 = c("T1", "T2", "T3", "T4", "T5", "T6"),
  C5 = c("T2", "T3", "T4"),
  C6 = c("T3", "T4", "T5"),
  C7 = c("T4", "T5", "T6"),
  C8 = c("T3", "T4", "T6")
)

set.seed(2024)
n_repeats <- 50
n_folds <- 10
alpha_value <- 0.5
lambda_rule <- "lambda.min"
gra_rho <- 0.5
selection_metric <- "RMSE"  # available options: "RMSE" or "MAE"

###############################################################################
# Helper functions
###############################################################################

to_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

normalize01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || abs(rng[2] - rng[1]) < .Machine$double.eps) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

r2_lm_style <- function(obs, pred) {
  obs <- as.numeric(obs)
  pred <- as.numeric(pred)
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (length(obs) < 3 || sd(obs) == 0 || sd(pred) == 0) return(NA_real_)
  summary(lm(obs ~ pred))$r.squared
}

calc_metrics <- function(obs, pred) {
  obs <- as.numeric(obs)
  pred <- as.numeric(pred)
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  c(
    R2 = r2_lm_style(obs, pred),
    RMSE = rmse,
    RRMSE = rmse / mean(obs) * 100,
    MAE = mae
  )
}

linear_fitted_yield <- function(y, x) {
  df <- data.frame(y = as.numeric(y), x = as.numeric(x))
  fit <- lm(y ~ x, data = df)
  as.numeric(predict(fit, newdata = df))
}

gra_rank <- function(y, X, rho = 0.5) {
  fitted_mat <- sapply(seq_len(ncol(X)), function(j) {
    linear_fitted_yield(y, X[, j])
  })
  if (is.null(dim(fitted_mat))) fitted_mat <- matrix(fitted_mat, ncol = 1)

  y_norm <- normalize01(y)
  fit_norm <- apply(fitted_mat, 2, normalize01)
  if (is.null(dim(fit_norm))) fit_norm <- matrix(fit_norm, ncol = 1)

  diff_mat <- abs(sweep(fit_norm, 1, y_norm, "-"))
  d_min <- min(diff_mat, na.rm = TRUE)
  d_max <- max(diff_mat, na.rm = TRUE)
  coeff <- (d_min + rho * d_max) / (diff_mat + rho * d_max)
  grd <- colMeans(coeff, na.rm = TRUE)

  data.frame(
    Rank = seq_along(grd),
    Feature = colnames(X),
    GRD = as.numeric(grd),
    stringsAsFactors = FALSE
  )[order(-grd), ]
}

make_foldid <- function(y, k = 10, seed_value = 1) {
  set.seed(seed_value)
  folds <- caret::createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  foldid <- integer(length(y))
  for (i in seq_along(folds)) foldid[folds[[i]]] <- i
  foldid
}

preprocess_xy <- function(X_train, X_test) {
  pp <- caret::preProcess(as.data.frame(X_train), method = c("center", "scale"))
  list(
    X_train = as.matrix(predict(pp, as.data.frame(X_train))),
    X_test = as.matrix(predict(pp, as.data.frame(X_test)))
  )
}

fit_predict_enr <- function(X_train, y_train, X_test, foldid_inner,
                            alpha = 0.5, lambda_rule = "lambda.min") {
  pp <- preprocess_xy(X_train, X_test)

  if (ncol(pp$X_train) == 1) {
    fit <- cv.glmnet(
      x = pp$X_train,
      y = y_train,
      alpha = alpha,
      foldid = foldid_inner,
      standardize = FALSE
    )
  } else {
    fit <- cv.glmnet(
      x = pp$X_train,
      y = y_train,
      alpha = alpha,
      foldid = foldid_inner,
      standardize = FALSE
    )
  }

  selected_lambda <- if (lambda_rule == "lambda.1se") fit$lambda.1se else fit$lambda.min

  list(
    pred_train = as.numeric(predict(fit, newx = pp$X_train, s = selected_lambda)),
    pred_test = as.numeric(predict(fit, newx = pp$X_test, s = selected_lambda)),
    lambda_min = fit$lambda.min,
    lambda_1se = fit$lambda.1se,
    selected_lambda = selected_lambda,
    n_predictors = ncol(pp$X_train)
  )
}

entropy_weights <- function(metric_matrix) {
  M <- as.matrix(metric_matrix)
  M <- apply(M, 2, normalize01)
  if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  M <- M + 1e-12
  P <- sweep(M, 2, colSums(M), "/")
  n <- nrow(P)
  e <- -colSums(P * log(P)) / log(n)
  d <- 1 - e
  if (sum(d) == 0 || any(!is.finite(d))) {
    rep(1 / ncol(M), ncol(M))
  } else {
    as.numeric(d / sum(d))
  }
}

###############################################################################
# Read data
###############################################################################

sheet_names <- excel_sheets(input_file)
if (length(sheet_names) < 6) {
  stop("The input workbook must contain at least six sheets.")
}

raw_data <- list()
for (i in seq_along(stage_names)) {
  dat <- as.data.frame(read_excel(input_file, sheet = sheet_names[i]))
  names(dat) <- make.names(names(dat), unique = TRUE)
  dat[] <- lapply(dat, to_numeric)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  raw_data[[stage_names[i]]] <- dat
}

y <- raw_data[[1]][[1]]
if (any(sapply(raw_data, nrow) != length(y))) {
  stop("All six sheets must contain the same number of observations after removing incomplete rows.")
}

X_stage <- list()
vi_names <- list()
nrct_name <- list()

for (s in stage_names) {
  dat <- raw_data[[s]]
  X <- dat[, -1, drop = FALSE]
  nrct_col <- names(X)[ncol(X)]
  vi_cols <- names(X)[seq_len(ncol(X) - 1)]
  X_stage[[s]] <- X
  vi_names[[s]] <- vi_cols
  nrct_name[[s]] <- nrct_col
}

###############################################################################
# GRA ranking and feature selection
###############################################################################

select_features <- function(stage, include_nrct = FALSE) {
  X <- X_stage[[stage]]
  vi_cols <- vi_names[[stage]]
  nrct_col <- nrct_name[[stage]]

  if (include_nrct) {
    X_for_gra <- X[, c(vi_cols, nrct_col), drop = FALSE]
  } else {
    X_for_gra <- X[, vi_cols, drop = FALSE]
  }

  ranking <- gra_rank(y, X_for_gra, rho = gra_rho)
  ranking$Stage <- stage
  ranking$StageLabel <- stage_labels[[stage]]
  ranking$InputSet <- if (include_nrct) "VI+NRCT" else "VI only"

  if (include_nrct) {
    vi_order <- ranking$Feature[ranking$Feature != nrct_col]
    subset_list <- lapply(seq_along(vi_order), function(k) {
      unique(c(nrct_col, vi_order[seq_len(k)]))
    })
  } else {
    vi_order <- ranking$Feature
    subset_list <- lapply(seq_along(vi_order), function(k) {
      vi_order[seq_len(k)]
    })
  }

  full_foldid <- make_foldid(y, k = n_folds, seed_value = 1000 + match(stage, stage_names))
  curve <- list()

  for (k in seq_along(subset_list)) {
    feats <- subset_list[[k]]
    pred <- fit_predict_enr(
      X_train = X[, feats, drop = FALSE],
      y_train = y,
      X_test = X[, feats, drop = FALSE],
      foldid_inner = full_foldid,
      alpha = alpha_value,
      lambda_rule = lambda_rule
    )
    met <- calc_metrics(y, pred$pred_train)
    curve[[k]] <- data.frame(
      Stage = stage,
      StageLabel = stage_labels[[stage]],
      InputSet = if (include_nrct) "VI+NRCT" else "VI only",
      Step = k,
      FeatureCount = length(feats),
      AddedFeature = if (include_nrct) vi_order[k] else feats[length(feats)],
      Features = paste(feats, collapse = ";"),
      R2 = met["R2"],
      RMSE = met["RMSE"],
      RRMSE = met["RRMSE"],
      MAE = met["MAE"],
      stringsAsFactors = FALSE
    )
  }

  curve_df <- do.call(rbind, curve)
  best_i <- if (selection_metric == "MAE") {
    which.min(curve_df$MAE)
  } else {
    which.min(curve_df$RMSE)
  }

  list(
    ranking = ranking,
    curve = curve_df,
    selected_features = strsplit(curve_df$Features[best_i], ";", fixed = TRUE)[[1]],
    best_step = curve_df$Step[best_i],
    best_metric_value = if (selection_metric == "MAE") curve_df$MAE[best_i] else curve_df$RMSE[best_i]
  )
}

rankings <- list()
curves <- list()
selected_vi <- list()
selected_vinrct <- list()
selected_records <- list()

for (s in stage_names) {
  sel_vi <- select_features(s, include_nrct = FALSE)
  sel_vinrct <- select_features(s, include_nrct = TRUE)

  rankings[[paste0(s, "_VI")]] <- sel_vi$ranking
  rankings[[paste0(s, "_VI_NRCT")]] <- sel_vinrct$ranking
  curves[[paste0(s, "_VI")]] <- sel_vi$curve
  curves[[paste0(s, "_VI_NRCT")]] <- sel_vinrct$curve

  selected_vi[[s]] <- sel_vi$selected_features
  selected_vinrct[[s]] <- sel_vinrct$selected_features

  selected_records[[length(selected_records) + 1]] <- data.frame(
    Stage = s,
    StageLabel = stage_labels[[s]],
    InputSet = "VI only",
    SelectedFeatures = paste(sel_vi$selected_features, collapse = ";"),
    FeatureCount = length(sel_vi$selected_features),
    BestStep = sel_vi$best_step,
    SelectionMetric = selection_metric,
    TrainingError = sel_vi$best_metric_value,
    stringsAsFactors = FALSE
  )

  selected_records[[length(selected_records) + 1]] <- data.frame(
    Stage = s,
    StageLabel = stage_labels[[s]],
    InputSet = "VI+NRCT",
    SelectedFeatures = paste(sel_vinrct$selected_features, collapse = ";"),
    FeatureCount = length(sel_vinrct$selected_features),
    BestStep = sel_vinrct$best_step,
    SelectionMetric = selection_metric,
    TrainingError = sel_vinrct$best_metric_value,
    stringsAsFactors = FALSE
  )
}

gra_df <- do.call(rbind, rankings)
curve_df <- do.call(rbind, curves)
selected_df <- do.call(rbind, selected_records)

###############################################################################
# 50 repeated 10-fold cross-validation
###############################################################################

single_stage_metrics <- list()
fusion_metrics <- list()
entropy_weight_records <- list()
hyper_records <- list()

record_hyper <- function(store, repeat_id, fold_id, model_name, stage, combo, pred_obj) {
  store[[length(store) + 1]] <- data.frame(
    Repeat = repeat_id,
    Fold = fold_id,
    Model = model_name,
    Stage = ifelse(is.null(stage), NA, stage),
    Combination = ifelse(is.null(combo), NA, combo),
    Alpha = alpha_value,
    LambdaMin = pred_obj$lambda_min,
    Lambda1SE = pred_obj$lambda_1se,
    SelectedLambda = pred_obj$selected_lambda,
    NumPredictors = pred_obj$n_predictors,
    stringsAsFactors = FALSE
  )
  store
}

for (r in seq_len(n_repeats)) {
  foldid_outer <- make_foldid(y, k = n_folds, seed_value = 10000 + r)

  for (fold in seq_len(n_folds)) {
    test_idx <- which(foldid_outer == fold)
    train_idx <- setdiff(seq_along(y), test_idx)
    y_train <- y[train_idx]
    y_test <- y[test_idx]
    inner_foldid <- make_foldid(y_train, k = n_folds, seed_value = 20000 + r * 100 + fold)

    pred_train_stage <- list()
    pred_test_stage <- list()

    for (s in stage_names) {
      pred_vi <- fit_predict_enr(
        X_train = X_stage[[s]][train_idx, selected_vi[[s]], drop = FALSE],
        y_train = y_train,
        X_test = X_stage[[s]][test_idx, selected_vi[[s]], drop = FALSE],
        foldid_inner = inner_foldid,
        alpha = alpha_value,
        lambda_rule = lambda_rule
      )
      met_vi <- calc_metrics(y_test, pred_vi$pred_test)
      hyper_records <- record_hyper(hyper_records, r, fold, paste0("VI_", s), s, NULL, pred_vi)

      single_stage_metrics[[length(single_stage_metrics) + 1]] <- data.frame(
        Repeat = r,
        Fold = fold,
        Stage = s,
        StageLabel = stage_labels[[s]],
        InputSet = "VI only",
        R2 = met_vi["R2"],
        RMSE = met_vi["RMSE"],
        RRMSE = met_vi["RRMSE"],
        MAE = met_vi["MAE"],
        stringsAsFactors = FALSE
      )

      pred_vinrct <- fit_predict_enr(
        X_train = X_stage[[s]][train_idx, selected_vinrct[[s]], drop = FALSE],
        y_train = y_train,
        X_test = X_stage[[s]][test_idx, selected_vinrct[[s]], drop = FALSE],
        foldid_inner = inner_foldid,
        alpha = alpha_value,
        lambda_rule = lambda_rule
      )
      met_vinrct <- calc_metrics(y_test, pred_vinrct$pred_test)
      hyper_records <- record_hyper(hyper_records, r, fold, paste0("VI_NRCT_", s), s, NULL, pred_vinrct)

      single_stage_metrics[[length(single_stage_metrics) + 1]] <- data.frame(
        Repeat = r,
        Fold = fold,
        Stage = s,
        StageLabel = stage_labels[[s]],
        InputSet = "VI+NRCT",
        R2 = met_vinrct["R2"],
        RMSE = met_vinrct["RMSE"],
        RRMSE = met_vinrct["RRMSE"],
        MAE = met_vinrct["MAE"],
        stringsAsFactors = FALSE
      )

      pred_train_stage[[s]] <- pred_vinrct$pred_train
      pred_test_stage[[s]] <- pred_vinrct$pred_test
    }

    for (combo_name in names(fusion_combos)) {
      stages <- fusion_combos[[combo_name]]

      Xtr_list <- list()
      Xte_list <- list()
      for (s in stages) {
        feats <- selected_vinrct[[s]]
        Xtr <- X_stage[[s]][train_idx, feats, drop = FALSE]
        Xte <- X_stage[[s]][test_idx, feats, drop = FALSE]
        colnames(Xtr) <- paste0(s, "_", colnames(Xtr))
        colnames(Xte) <- paste0(s, "_", colnames(Xte))
        Xtr_list[[s]] <- Xtr
        Xte_list[[s]] <- Xte
      }

      Xtr_fusion <- do.call(cbind, Xtr_list)
      Xte_fusion <- do.call(cbind, Xte_list)

      pred_feature <- fit_predict_enr(
        X_train = Xtr_fusion,
        y_train = y_train,
        X_test = Xte_fusion,
        foldid_inner = inner_foldid,
        alpha = alpha_value,
        lambda_rule = lambda_rule
      )
      met_feature <- calc_metrics(y_test, pred_feature$pred_test)
      hyper_records <- record_hyper(hyper_records, r, fold, "Feature_level_fusion", NULL, combo_name, pred_feature)

      fusion_metrics[[length(fusion_metrics) + 1]] <- data.frame(
        Repeat = r,
        Fold = fold,
        Combination = combo_name,
        GrowthStages = paste(stage_labels[stages], collapse = "+"),
        FusionMethod = "Feature-level fusion",
        R2 = met_feature["R2"],
        RMSE = met_feature["RMSE"],
        RRMSE = met_feature["RRMSE"],
        MAE = met_feature["MAE"],
        stringsAsFactors = FALSE
      )

      train_pred_matrix <- do.call(cbind, pred_train_stage[stages])
      test_pred_matrix <- do.call(cbind, pred_test_stage[stages])
      colnames(train_pred_matrix) <- stages
      colnames(test_pred_matrix) <- stages

      w <- entropy_weights(train_pred_matrix)
      ewf_pred <- as.numeric(test_pred_matrix %*% w)
      met_ewf <- calc_metrics(y_test, ewf_pred)

      fusion_metrics[[length(fusion_metrics) + 1]] <- data.frame(
        Repeat = r,
        Fold = fold,
        Combination = combo_name,
        GrowthStages = paste(stage_labels[stages], collapse = "+"),
        FusionMethod = "Entropy weight fusion",
        R2 = met_ewf["R2"],
        RMSE = met_ewf["RMSE"],
        RRMSE = met_ewf["RRMSE"],
        MAE = met_ewf["MAE"],
        stringsAsFactors = FALSE
      )

      entropy_weight_records[[length(entropy_weight_records) + 1]] <- data.frame(
        Repeat = r,
        Fold = fold,
        Combination = combo_name,
        Stage = stages,
        StageLabel = unname(stage_labels[stages]),
        Weight = w,
        stringsAsFactors = FALSE
      )
    }
  }
}

single_stage_df <- do.call(rbind, single_stage_metrics)
fusion_df <- do.call(rbind, fusion_metrics)
weights_df <- do.call(rbind, entropy_weight_records)
hyper_df <- do.call(rbind, hyper_records)

###############################################################################
# Summaries and output
###############################################################################

single_stage_summary <- aggregate(
  cbind(R2, RMSE, RRMSE, MAE) ~ Stage + StageLabel + InputSet,
  data = single_stage_df,
  FUN = mean
)

fusion_summary <- aggregate(
  cbind(R2, RMSE, RRMSE, MAE) ~ Combination + GrowthStages + FusionMethod,
  data = fusion_df,
  FUN = mean
)

hyper_summary <- aggregate(
  cbind(SelectedLambda, NumPredictors) ~ Model,
  data = hyper_df,
  FUN = function(z) c(Min = min(z, na.rm = TRUE), Median = median(z, na.rm = TRUE), Max = max(z, na.rm = TRUE))
)

write.csv(gra_df, file.path(output_dir, "GRA_rankings.csv"), row.names = FALSE)
write.csv(curve_df, file.path(output_dir, "Feature_selection_training_error_curves.csv"), row.names = FALSE)
write.csv(selected_df, file.path(output_dir, "Selected_features.csv"), row.names = FALSE)
write.csv(single_stage_df, file.path(output_dir, "Single_stage_50x10CV_metrics.csv"), row.names = FALSE)
write.csv(single_stage_summary, file.path(output_dir, "Single_stage_mean_metrics.csv"), row.names = FALSE)
write.csv(fusion_df, file.path(output_dir, "Fusion_50x10CV_metrics.csv"), row.names = FALSE)
write.csv(fusion_summary, file.path(output_dir, "Fusion_mean_metrics.csv"), row.names = FALSE)
write.csv(weights_df, file.path(output_dir, "Entropy_weight_records.csv"), row.names = FALSE)
write.csv(hyper_df, file.path(output_dir, "ENR_hyperparameter_records.csv"), row.names = FALSE)

write_xlsx(
  list(
    GRA_rankings = gra_df,
    Feature_selection_curves = curve_df,
    Selected_features = selected_df,
    Single_stage_50x10CV_metrics = single_stage_df,
    Single_stage_mean_metrics = single_stage_summary,
    Fusion_50x10CV_metrics = fusion_df,
    Fusion_mean_metrics = fusion_summary,
    Entropy_weight_records = weights_df,
    ENR_hyperparameter_records = hyper_df
  ),
  path = file.path(output_dir, "Full_workflow_results.xlsx")
)

cat("\nWorkflow finished.\n")
cat("Outputs saved to:\n")
cat(output_dir, "\n")
