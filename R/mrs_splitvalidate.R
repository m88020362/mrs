#' @title Perform split-validation on MRS model to select best lambda
#' @name mrs_splitvalidate
#' @description
#' Split-validation procedure to select the best-performing lambda based on validation performance.
#' For binary phenotypes, selection can be based on AUC, accuracy, or F1-score.
#' For continuous phenotypes, correlation ("cor") is used. The held-out test set is then used to
#' compute corresponding metrics.
#'
#' @param risk_mat A numeric matrix of risk scores. Rows are samples, columns are lambdas.
#' @param test_data A data.frame with a `Class` column (phenotype: 0/1 or continuous).
#' @param lambda A numeric vector of lambda values; must match `ncol(risk_mat)`.
#' @param seed Optional integer to fix random split for reproducibility.
#' @param select_by One of "auc", "accuracy", "f1", or "cor". Determines how the best lambda is selected.
#'
#' @return A list with elements: `best_lambda`, `best_lambda_index`, `test_auc`, `test_accuracy`,
#' `test_f1`, `test_correlation`, `validation_correlation`, `metric_matrix`, and `results_table`.
#'
#' @importFrom pROC roc auc
#' @importFrom stats cor
#' @export
mrs_splitvalidate <- function(risk_mat, test_data, lambda, seed = NULL, select_by = "auc") {
  stopifnot(nrow(risk_mat) == nrow(test_data))
  stopifnot(length(lambda) == ncol(risk_mat))
  stopifnot("Class" %in% colnames(test_data))

  if (!is.null(seed)) set.seed(seed)
  n <- nrow(test_data)
  split <- sample(1:2, n, replace = TRUE)

  pheno <- test_data$Class
  val_idx <- which(split == 1)
  test_idx <- which(split == 2)

  risk_val <- risk_mat[val_idx, , drop = FALSE]
  pheno_val <- pheno[val_idx]

  metric_to_use <- tolower(select_by)

  if (metric_to_use == "cor") {
    # For continuous phenotype
    cor_vec <- apply(risk_val, 2, function(x) cor(x, pheno_val, use = "complete.obs"))
    best_lambda_index <- which.max(abs(cor_vec))
    best_lambda <- lambda[best_lambda_index]
    risk_test <- risk_mat[test_idx, best_lambda_index]
    pheno_test <- pheno[test_idx]
    test_cor <- cor(risk_test, pheno_test, use = "complete.obs")

    result_df <- data.frame(
      ID = seq_len(n),
      pheno = pheno,
      best_pgs = NA_real_,
      split = split
    )
    result_df$best_pgs[test_idx] <- risk_test

    return(list(
      best_lambda = best_lambda,
      best_lambda_index = best_lambda_index,
      validation_correlation = cor_vec,
      test_correlation = test_cor,
      results_table = result_df
    ))
  } else {
    # For binary phenotype
    get_metrics <- function(score, truth) {
      pred_bin <- ifelse(score > 0, 1, 0)
      auc <- suppressWarnings(pROC::auc(pROC::roc(truth, score)))
      acc <- mean(pred_bin == truth)
      tp <- sum(pred_bin == 1 & truth == 1)
      fp <- sum(pred_bin == 1 & truth == 0)
      fn <- sum(pred_bin == 0 & truth == 1)
      f1 <- if ((2 * tp + fp + fn) == 0) NA else (2 * tp) / (2 * tp + fp + fn)
      return(c(auc = auc, accuracy = acc, f1 = f1))
    }

    metric_mat <- t(apply(risk_val, 2, get_metrics, truth = pheno_val))
    if (!(metric_to_use %in% colnames(metric_mat))) stop("Invalid select_by metric for binary outcome")

    best_lambda_index <- which.max(metric_mat[, metric_to_use])
    best_lambda <- lambda[best_lambda_index]

    risk_test <- risk_mat[test_idx, best_lambda_index]
    pheno_test <- pheno[test_idx]
    pred_test <- ifelse(risk_test > 0, 1, 0)

    auc_test <- suppressWarnings(pROC::auc(pROC::roc(pheno_test, risk_test)))
    acc_test <- mean(pred_test == pheno_test)
    tp <- sum(pred_test == 1 & pheno_test == 1)
    fp <- sum(pred_test == 1 & pheno_test == 0)
    fn <- sum(pred_test == 0 & pheno_test == 1)
    f1_test <- if ((2 * tp + fp + fn) == 0) NA else (2 * tp) / (2 * tp + fp + fn)

    result_df <- data.frame(
      ID = seq_len(n),
      pheno = pheno,
      best_pgs = NA_real_,
      split = split
    )
    result_df$best_pgs[test_idx] <- risk_test

    return(list(
      best_lambda = best_lambda,
      best_lambda_index = best_lambda_index,
      metric_matrix = metric_mat,
      selected_metric = metric_to_use,
      test_auc = auc_test,
      test_accuracy = acc_test,
      test_f1 = f1_test,
      results_table = result_df
    ))
  }
}



