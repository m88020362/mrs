#' @title Flexible MRS Pipeline for Binary or Continuous Outcome
#' @description
#' Generalized pipeline to construct MRS for either binary or continuous phenotype.
#' Allows shrinkage, HLO-based scoring, and split-validation with metric selection.
#'
#' @param base_output MB-MDR output data frame with columns "pair" and "F_value".
#' @param test_data Data frame with phenotype in column `Class` and SNPs in others.
#' @param raw_model Parsed MB-MDR .model file.
#' @param lambda Numeric vector of penalty parameters.
#' @param base_data Optional base data to determine observed genotypes in parsing.
#' @param select_by Metric used for lambda selection. One of: "cor", "auc", "accuracy", "f1".
#' @param seed Optional seed for split reproducibility.
#' @param verbose If TRUE, prints details during scoring.
#' @return A list containing validation results and test performance.
#' @export
mrs.pipeline <- function(
    base_output,
    test_data,
    raw_model,
    lambda = exp(seq(log(0.001), log(0.1), length.out = 20)),
    base_data = NULL,
    select_by = "auc",
    seed = NULL,
    verbose = FALSE
) {
  # Step 1: shrinkage of F-values
  beta_obj <- mrs_indeplasso(base_output, lambda)

  # Step 2: parse HLO matrices
  if (is.null(base_data)) {
    base_data <- test_data[, !(colnames(test_data) %in% "Class"), drop = FALSE]
  }
  mdr_model <- parse_model_to_mdr(raw_model, base_data = base_data)

  # Step 3: risk matrix construction
  risk_mat <- compute_one_mrs_fast(
    snp_data = test_data[, !(colnames(test_data) %in% "Class"), drop = FALSE],
    mdr_model = mdr_model,
    beta = beta_obj$beta,
    pair_names = beta_obj$pair_names
  )

  # Step 4: split-validation using specified metric
  res <- mrs_splitvalidate(
    risk_mat = risk_mat,
    test_data = test_data,
    lambda = beta_obj$lambda,
    seed = seed,
    select_by = select_by
  )

  return(res)
}


