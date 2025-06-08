#' @title Plot risk score distribution by phenotype (test set)
#' @name plot_mrs_score_dist
#' @description Violin + box plot showing score distribution across phenotypes.
#'
#' @param res Result list from mrs_splitvalidate() with results_table.
#'
#' @return A ggplot violin plot annotated with AUC and Accuracy.
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme_minimal
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' @importFrom graphics legend
#' @importFrom utils globalVariables
#' @export
utils::globalVariables(c("pheno", "best_pgs"))

plot_mrs_score_dist <- function(res) {
  df <- res$results_table
  test_df <- df[df$split == 2, ]

  auc <- round(res$test_auc, 3)
  acc <- round(res$test_accuracy, 3)

  ggplot(test_df, aes(x = as.factor(pheno), y = best_pgs)) +
    geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(
      title = paste0(
        "Risk Score Distribution at Best lambda = ", res$best_lambda,
        "\nAUC = ", auc, " | Accuracy = ", acc
      ),
      x = "Phenotype (0 = control, 1 = case)",
      y = "Risk Score"
    ) +
    theme_minimal()
}


