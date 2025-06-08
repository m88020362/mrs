#' @title Plot ROC curve on test set using best lambda
#' @description Draw ROC curve using predicted risk scores from `mrs_splitvalidate()` result.
#'
#' @param res Result list from `mrs_splitvalidate()`, must contain `results_table` and `best_lambda`.
#'
#' @return A ggplot ROC curve with AUC annotation.
#'
#' @import ggplot2
#' @importFrom graphics abline legend
#' @importFrom pROC roc auc
#' @export
plot_mrs_roc <- function(res) {
  df <- res$results_table
  test_idx <- df$split == 2
  true <- df$pheno[test_idx]
  score <- df$best_pgs[test_idx]

  roc_obj <- pROC::roc(response = true, predictor = score)
  auc_val <- pROC::auc(roc_obj)

  plot_df <- data.frame(
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities)
  )

  ggplot(plot_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(size = 1.2) +
    geom_abline(linetype = "dashed", color = "gray") +
    labs(
      title = paste0("ROC Curve (AUC = ", round(auc_val, 3), ")"),
      x = "1 - Specificity", y = "Sensitivity"
    ) +
    theme_minimal()
}


