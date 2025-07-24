#' @title Compute full MRS matrix across all individuals (fast version)
#' @description Efficiently computes MRS scores using cleaned genotype data and MB-MDR HLO matrices.
#' @param snp_data Data.frame or matrix. Genotype data (0/1/2) with SNPs as columns.
#' @param mdr_model Named list. Each element is a 3x3 HLO matrix.
#' @param beta Matrix of F-value weights. Rows = SNP pairs, columns = lambda.
#' @param pair_names Character vector. Names of SNP pairs in same order as beta.
#' @return A numeric matrix (n_individuals x n_lambda) of MRS scores.
#' @details Assumes clean input. Invalid values should be removed before use.
#' @export
compute_one_mrs_fast <- function(snp_data, mdr_model, beta, pair_names) {

  split_pair <- strsplit(pair_names, "_")
  snp1_vec <- sapply(split_pair, function(x) paste(x[1:(length(x) - 2)], collapse = "_"))
  snp2_vec <- sapply(split_pair, function(x) paste(x[(length(x) - 1):length(x)], collapse = "_"))

  rs_mat <- matrix(0, nrow = nrow(snp_data), ncol = ncol(beta))

  for (j in seq_along(pair_names)) {
    snp1 <- snp1_vec[j]
    snp2 <- snp2_vec[j]
    pair <- pair_names[j]

    if (!(snp1 %in% colnames(snp_data)) || !(snp2 %in% colnames(snp_data))) next
    if (!(pair %in% names(mdr_model))) next

    geno1 <- as.integer(snp_data[[snp1]])
    geno2 <- as.integer(snp_data[[snp2]])
    valid_idx <- which(!is.na(geno1) & !is.na(geno2) & geno1 %in% 0:2 & geno2 %in% 0:2)

    if (length(valid_idx) == 0) next

    hlo <- mdr_model[[pair]]
    if (!is.matrix(hlo) || any(dim(hlo) != c(3, 3))) next

    risk_contrib <- hlo[cbind(geno1[valid_idx] + 1, geno2[valid_idx] + 1)]
    rs_mat[valid_idx, ] <- rs_mat[valid_idx, ] + risk_contrib * rep(1, ncol(beta)) * rep(beta[j, ], each = length(valid_idx))
  }

  return(rs_mat)
}
