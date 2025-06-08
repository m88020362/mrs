#' @title Full MRS Pipeline Using MB-MDR Output and Raw Model File
#' @description
#' One-step modeling process for multilocus risk score (MRS)
#' construction using MB-MDR results. Performs shrinkage, risk scoring, and validation.
#'
#' @param base_output A data frame with MB-MDR results. Must contain `pair` and `F_value` columns.
#' @param test_data A data frame with phenotype in column 1 and SNP genotypes in the remaining columns.
#' @param raw_model A matrix or data frame representing the MB-MDR `.model.txt` output.
#' @param lambda Numeric vector of penalty values for soft-thresholding.
#' @param seed Optional integer for reproducible split-validation.
#' @return A list containing:
#' \describe{
#'   \item{best_lambda}{Lambda with best validation performance.}
#'   \item{test_correlation}{Correlation between MRS and phenotype.}
#'   \item{test_auc}{AUC in test set.}
#'   \item{test_accuracy}{Accuracy in test set.}
#'   \item{test_f1}{F1-score in test set.}
#'   \item{validation_correlation}{Correlation for each lambda in base split.}
#'   \item{results_table}{Data frame with phenotype, split label, and MRS.}
#' }
#' @seealso \code{\link{mrs_indeplasso}}, \code{\link{parse_model_to_mdr}}, \code{\link{mrs.default}}, \code{\link{mrs_splitvalidate}}
#' @importFrom stats cor
#' @importFrom dplyr %>% mutate
#' @export
mrs.pipeline <- function(
    base_output,
    test_data,
    raw_model,
    lambda = exp(seq(log(0.001), log(0.1), length.out = 20)),
    seed = NULL
) {
  beta_obj <- mrs_indeplasso(base_output, lambda)

  mdr_model <- parse_model_to_mdr(raw_model)

  risk_mat <- mrs.default(
    snp_data = test_data[, -1],
    mdr_model = mdr_model,
    beta = beta_obj$beta,
    pair_names = beta_obj$pair_names
  )

  res <- mrs_splitvalidate(
    risk_mat = risk_mat,
    test_data = test_data,
    lambda = beta_obj$lambda,
    seed = seed
  )

  return(res)
}


