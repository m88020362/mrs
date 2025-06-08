#' @title Statistical test for group separation of MRS
#' @name mrs_group_test
#' @description Performs Wilcoxon rank-sum test on MRS scores between case and control groups.
#'
#' @param res Result object from `mrs_splitvalidate()`, containing `results_table` and phenotype.
#'
#' @return A named list with test method and p-value.
#'
#' @importFrom stats wilcox.test
#' @export
mrs_group_test <- function(res) {
  df <- res$results_table
  test_df <- df[df$split == 2, ]
  out <- wilcox.test(best_pgs ~ pheno, data = test_df)
  return(list(
    method = out$method,
    p_value = out$p.value
  ))
}

