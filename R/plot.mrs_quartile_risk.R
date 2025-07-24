#' @title Quartile Odds Ratio Plot for MRS
#' @description
#' Visualizes odds ratios across MRS quartiles in the target set.
#' Supports reverse ranking for L-dominant encodings.
#'
#' @param validate_result Output from \code{mrs_splitvalidate()}.
#' @param reverse_rank Logical. If TRUE, rank by \code{-MRS} instead of \code{MRS}. Default: TRUE.
#'
#' @return Invisibly returns the ggplot object (invisible(g)), but renders it by default.
#' @export
plot.mrs_quartile_risk <- function(validate_result, reverse_rank = TRUE) {
  df <- validate_result$results_table %>%
    dplyr::filter(split == 2, !is.na(best_pgs)) %>%
    dplyr::mutate(
      quartile = if (reverse_rank) ntile(-1 * best_pgs, 4) else ntile(best_pgs, 4),
      quartile = factor(quartile, labels = paste0("Q", 1:4))
    )

  model <- glm(pheno ~ quartile, data = df, family = binomial())
  or_table <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(grepl("quartile", term)) %>%
    dplyr::mutate(
      quartile = sub("quartile", "", term),
      quartile = factor(quartile, levels = c("Q2", "Q3", "Q4"))
    )

  or_table <- dplyr::bind_rows(
    tibble::tibble(
      term = "quartileQ1",
      estimate = 1,
      conf.low = 1,
      conf.high = 1,
      quartile = factor("Q1", levels = paste0("Q", 1:4))
    ),
    or_table
  )

  g <- ggplot2::ggplot(or_table, ggplot2::aes(x = quartile, y = estimate)) +
    ggplot2::geom_point(size = 4, shape = 21, fill = "black") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = conf.low, ymax = conf.high), width = 0.15, size = 1.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::scale_y_continuous(name = "Odds Ratio (MCI vs Normal)", trans = "log10") +
    ggplot2::xlab("MRS Quartile") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))
    ) +
    ggplot2::ggtitle("MRS Quartiles vs Odds of MCI")

  print(g)
  invisible(g)
}
