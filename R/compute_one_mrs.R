#' @title Compute one MRS per individual (safe version)
#' @description Safely compute a multilocus risk score (MRS) for a single individual using MB-MDR HLO and F-values.
#' @param subj_idx Integer. Row index for the individual.
#' @param snp_data Data.frame or matrix. Genotype data (0/1/2) with SNPs as columns.
#' @param mdr_model Named list. Each element is a 3x3 HLO matrix from MB-MDR.
#' @param beta Matrix of F-value weights. Rows = SNP pairs, columns = lambda.
#' @param pair_names Character vector. Names of SNP pairs in the same order as beta.
#' @param verbose Logical. If TRUE, prints debug info for first few iterations.
#' @return A numeric vector of MRS scores for one subject across lambdas.
#' @export
compute_one_mrs <- function(subj_idx, snp_data, mdr_model, beta, pair_names, verbose = FALSE) {
  rs <- numeric(ncol(beta))

  split_pair_from_pairname <- function(pairname) {
    parts <- strsplit(pairname, "_", fixed = TRUE)[[1]]
    len <- length(parts)
    if (len < 4) return(c(NA, NA))
    snp1 <- paste(parts[1:(len - 2)], collapse = "_")
    snp2 <- paste(parts[(len - 1):len], collapse = "_")
    return(c(snp1, snp2))
  }

  for (j in seq_along(pair_names)) {
    pair <- pair_names[j]
    snps <- tryCatch(split_pair_from_pairname(pair), error = function(e) return(c(NA, NA)))
    snp1 <- snps[1]; snp2 <- snps[2]

    if (is.na(snp1) || is.na(snp2)) next
    if (!(snp1 %in% colnames(snp_data))) next
    if (!(snp2 %in% colnames(snp_data))) next
    if (!(pair %in% names(mdr_model))) next

    snp1_val <- tryCatch(as.integer(snp_data[[snp1]][subj_idx]), error = function(e) NA)
    snp2_val <- tryCatch(as.integer(snp_data[[snp2]][subj_idx]), error = function(e) NA)
    if (is.na(snp1_val) || is.na(snp2_val)) next
    if (!(snp1_val %in% 0:2) || !(snp2_val %in% 0:2)) next

    hlo <- mdr_model[[pair]]
    if (!is.matrix(hlo) || any(dim(hlo) != c(3, 3))) next

    if (verbose && j <= 5 && subj_idx == 1) {
      cat("pair:", pair, "\n")
      cat("snp1:", snp1, "val:", snp1_val, "\n")
      cat("snp2:", snp2, "val:", snp2_val, "\n")
      cat("beta[j, ]:", beta[j, ], "\n")
      cat("HLO[snp1_val+1, snp2_val+1]:", hlo[snp1_val + 1, snp2_val + 1], "\n\n")
    }

    rs <- rs + hlo[snp1_val + 1, snp2_val + 1] * beta[j, ]
  }

  return(rs)
}
