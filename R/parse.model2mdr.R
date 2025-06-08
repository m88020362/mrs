#' @title Parse MB-MDR .model File to Structured MDR Model List
#' @description Parses a raw `.model.txt` file into a structured list containing affected/unaffected counts and HLO matrices.
#'
#' @param raw_model A `data.frame` read from a `.model.txt` file.
#' @param row_per_pair Number of rows per SNP pair. Default is 14.
#' @return A named list where each element is a SNP pair and includes:
#' \describe{
#'   \item{affected}{3-row matrix for affected counts.}
#'   \item{unaffected}{3-row matrix for unaffected counts.}
#'   \item{HLO}{3×3 matrix of -1, 0, 1 encoding.}
#' }
#' @examples
#' \dontrun{
#' raw_model <- fread("single_model_2D.txt", fill = TRUE, header = FALSE)
#' mdr_model <- parse_model_to_mdr(raw_model)
#' }
#' @importFrom dplyr mutate across recode
#' @importFrom magrittr %>%
#' @export
parse_model_to_mdr <- function(raw_model, row_per_pair = 14) {
  n_pairs <- floor(nrow(raw_model) / row_per_pair)
  mdr_model <- list()

  for (i in seq_len(n_pairs)) {
    idx <- (i - 1) * row_per_pair + 1
    pair_name <- paste(raw_model[idx, 1:2], collapse = "_")

    hlo_matrix <- raw_model[(idx + 10):(idx + 12), ] %>%
      as.data.frame() %>%
      mutate(across(everything(), as.character)) %>%
      mutate(across(everything(), ~ recode(.x, "L" = -1, "H" = 1, "O" = 0))) %>%
      as.matrix()

    mdr_model[[pair_name]] <- list(
      affected = raw_model[(idx + 2):(idx + 4), ],
      unaffected = raw_model[(idx + 6):(idx + 8), ],
      HLO = hlo_matrix
    )
  }

  return(mdr_model)
}

