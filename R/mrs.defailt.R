#' @title Multilocus Risk Score calculation using MB-MDR outputs
#' @name mrs.default
#' @description Calculates multilocus risk scores (MRS) per individual using F-weighted HLO matrices from MB-MDR.
#'
#' @param snp_data A genotype matrix or data.frame (rows = individuals, columns = SNPs).
#' @param mdr_model A named list of MB-MDR model objects, each containing an HLO matrix.
#' @param beta A numeric matrix of F-value weights. Rows = SNP pairs, columns = lambda.
#' @param pair_names A character vector of pair names matching beta row order.
#' @param cluster Optional parallel backend created by `parallel::makeCluster`.
#'
#' @return A numeric matrix: rows = individuals, columns = lambda. Each cell = MRS score.
#'
#' @importFrom parallel parLapply
#' @export
mrs.default <- function(snp_data, mdr_model, beta, pair_names, cluster = NULL) {
  n_subj <- nrow(snp_data)
  n_lambda <- ncol(beta)

  compute_one <- function(subj_idx) {
    rs <- numeric(n_lambda)
    for (j in seq_along(pair_names)) {
      pair <- pair_names[j]
      snps <- strsplit(pair, "_")[[1]]
      snp1_val <- snp_data[[snps[1]]][subj_idx]
      snp2_val <- snp_data[[snps[2]]][subj_idx]

      hlo <- mdr_model[[pair]]$HLO
      if (!is.null(hlo)) {
        rs <- rs + hlo[snp1_val + 1, snp2_val + 1] * beta[j, ]
      }
    }
    return(rs)
  }

  if (!is.null(cluster)) {
    result <- do.call(rbind, parallel::parLapply(cluster, seq_len(n_subj), compute_one))
  } else {
    result <- t(sapply(seq_len(n_subj), compute_one))
  }

  return(result)
}

