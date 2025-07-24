#' @title Parse MB-MDR .model File to Structured MDR Model List
#' @description Parses a raw `.model.txt` file into a structured list containing 3×3 HLO matrices derived from genotype data.
#'
#' @param raw_model A `data.frame` read from a `.model.txt` file.
#' @param base_data A genotype matrix (samples × SNPs) used to determine actual genotype combinations per SNP pair.
#' @param filter_zero Logical; whether to exclude HLO matrices with all-zero entries. Default is TRUE.
#'
#' @return A named list where each element is a SNP pair and includes:
#' \describe{
#'   \item{HLO}{3×3 matrix of -1, 0, 1 encoding derived from MB-MDR result.}
#' }
#'
#' @examples
#' \dontrun{
#' raw_model <- fread("model_2D.txt", fill = TRUE, header = FALSE)
#' mdr_model <- parse_model_to_mdr(raw_model, base_data)
#' }
#'
#' @importFrom dplyr mutate across recode everything
#' @importFrom progress progress_bar
#' @export
parse_model_to_mdr <- function(raw_model, base_data, filter_zero = TRUE) {
  raw_model <- as.data.frame(raw_model)
  anchor_rows <- which(raw_model[[1]] == "HLO" & raw_model[[2]] == "matrix:")
  mdr_model <- list()

  total <- length(anchor_rows)
  pb <- progress::progress_bar$new(
    format = "  [:bar] :current/:total (:percent) eta: :eta",
    total = total
  )
  start_time <- Sys.time()

  for (i in anchor_rows) {
    pb$tick()

    hlo_matrix_raw <- raw_model[(i + 1):(i + 3), , drop = FALSE] |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
      as.data.frame()

    hlo_matrix <- as.data.frame(
      hlo_matrix_raw[
        apply(hlo_matrix_raw, 1, function(row) any(row %in% c("H", "L", "O"))),
        , drop = FALSE
      ]
    )

    nrow_hlo <- nrow(hlo_matrix)
    ncol_hlo <- max(apply(hlo_matrix, 2, function(col) any(col %in% c("H", "L", "O"))))

    snp_row <- i - (2 + 2 * nrow_hlo + 1)
    if (snp_row < 1) next

    snp1 <- raw_model[[1]][snp_row]
    snp2 <- raw_model[[2]][snp_row]
    if (is.na(snp1) || is.na(snp2)) next
    pair_name <- paste0(snp1, "_", snp2)

    if (!(snp1 %in% colnames(base_data)) || !(snp2 %in% colnames(base_data))) next

    geno1 <- base_data[[snp1]]
    geno2 <- base_data[[snp2]]
    observed1 <- sort(unique(stats::na.omit(as.integer(geno1))))
    observed2 <- sort(unique(stats::na.omit(as.integer(geno2))))

    full_matrix <- matrix("O", nrow = 3, ncol = 3)
    for (r_idx in seq_len(nrow_hlo)) {
      for (c_idx in seq_len(ncol_hlo)) {
        val <- hlo_matrix[r_idx, c_idx]
        if (val %in% c("H", "L", "O")) {
          if (r_idx <= length(observed1) && c_idx <= length(observed2)) {
            g1 <- observed1[r_idx]
            g2 <- observed2[c_idx]
            full_matrix[g1 + 1, g2 + 1] <- val
          }
        }
      }
    }

    hlo_numeric <- apply(full_matrix, c(1, 2), function(x) dplyr::recode(x, "H" = 1L, "L" = -1L, "O" = 0L))

    if (!filter_zero || any(hlo_numeric != 0)) {
      mdr_model[[pair_name]] <- hlo_numeric
    }
  }

  elapsed <- Sys.time() - start_time
  message(sprintf("✓ Completed %d pairs in %.2f seconds.", length(mdr_model), as.numeric(elapsed, units = "secs")))
  return(mdr_model)
}
