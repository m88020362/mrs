#' @name mrs_indeplasso
#' @rdname mrs_indeplasso
#' @title Shrink F-values for MRS modeling using LASSO-style soft-thresholding
#' @description Performs lambda-based shrinkage on MB-MDR derived F-values to obtain beta weights for each SNP pair.
#' @param output_df A data.frame with columns `pair` (e.g., SNP1_SNP2) and `F_value`.
#' @param lambda A numeric vector of lambda shrinkage values. Default is 20 log-scaled values from 0.001 to 0.1.
#' @return A list containing:
#' \describe{
#'   \item{lambda}{The lambda values used.}
#'   \item{beta}{A numeric matrix of shrunken F-values (rows = SNP pairs, columns = lambda values).}
#'   \item{pair_names}{Character vector of SNP pair names corresponding to rows in `beta`.}
#' }
#' @examples
#' out <- mrs_indeplasso(output_2D)
#' head(out$beta)
#' @export
mrs_indeplasso <- function(output_df, lambda = exp(seq(log(0.001), log(0.1), length.out = 20))) {
  stopifnot(all(c("pair", "F_value") %in% colnames(output_df)))

  F <- output_df$F_value

  beta <- outer(F, lambda, function(f, l) sign(f) * pmax(abs(f) - l, 0))

  return(list(
    lambda = lambda,
    beta = beta,
    pair_names = output_df$pair
  ))
}

