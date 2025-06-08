#' @title Perform split-validation on MRS model to select best lambda
#' @name mrs_splitvalidate
#' @description
#' Split-validation procedure to select the best-performing lambda based on validation set correlation.
#' The held-out test set is then used to compute AUC, accuracy, and F1-score.
#'
#' @param risk_mat A numeric matrix of risk scores. Rows are samples, columns are lambdas.
#' @param test_data A data.frame with a `Class` column (phenotype: 0/1) and optional SNPs.
#' @param lambda A numeric vector of lambda values; must match `ncol(risk_mat)`.
#' @param seed Optional integer to fix random split for reproducibility.
#'
#' @return A list with elements: `best_lambda`, `test_auc`, `test_accuracy`, `test_f1`,
#' `test_correlation`, `validation_correlation`, `results_table`.
#'
#' @importFrom pROC roc auc
#' @importFrom stats cor
#' @export
mrs_splitvalidate <- function(risk_mat, test_data, lambda, seed = NULL) {
  stopifnot(nrow(risk_mat) == nrow(test_data))
  stopifnot(length(lambda) == ncol(risk_mat))
  stopifnot("Class" %in% colnames(test_data))

  if (!is.null(seed)) set.seed(seed)
  n <- nrow(test_data)
  split <- sample(1:2, n, replace = TRUE)

  pheno <- test_data$Class
  val_idx <- which(split == 1)
  test_idx <- which(split == 2)

  # Validation
  risk_val <- risk_mat[val_idx, , drop = FALSE]
  pheno_val <- pheno[val_idx]
  cors <- apply(risk_val, 2, function(x) cor(x, pheno_val, use = "complete.obs"))
  best_lambda_index <- which.max(cors)
  best_lambda <- lambda[best_lambda_index]

  # Testing
  risk_test <- risk_mat[test_idx, best_lambda_index]
  pheno_test <- pheno[test_idx]
  test_cor <- cor(risk_test, pheno_test, use = "complete.obs")

  pred_test <- ifelse(risk_test > 0, 1, 0)

  # Metrics
  test_auc <- suppressWarnings(pROC::auc(pROC::roc(pheno_test, risk_test)))
  test_acc <- mean(pred_test == pheno_test)
  tp <- sum(pred_test == 1 & pheno_test == 1)
  fp <- sum(pred_test == 1 & pheno_test == 0)
  fn <- sum(pred_test == 0 & pheno_test == 1)
  test_f1 <- if ((2 * tp + fp + fn) == 0) NA else (2 * tp) / (2 * tp + fp + fn)

  result_df <- data.frame(
    ID = seq_len(n),
    pheno = pheno,
    best_pgs = NA_real_,
    split = split
  )
  result_df$best_pgs[test_idx] <- risk_test

  return(list(
    best_lambda = best_lambda,
    test_correlation = test_cor,
    test_auc = test_auc,
    test_accuracy = test_acc,
    test_f1 = test_f1,
    validation_correlation = cors,
    results_table = result_df
  ))
}



