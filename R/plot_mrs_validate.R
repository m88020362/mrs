#' @title Plot validation results from mrs_splitvalidate
#' @description Plots lambda vs correlation.
#'
#' @param result A result list from `mrs_splitvalidate()`.
#' @param metric One of `"correlation"` (default), `"auc"`, `"accuracy"`, `"f1"`.
#' @param log_lambda Logical, whether to log-transform x-axis.
#'
#' @return A base R plot of validation metric vs. lambda.
#'
#' @importFrom graphics abline legend
#' @export
plot_mrs_validate <- function(result, metric = "correlation", log_lambda = TRUE) {
  stopifnot(metric %in% c("correlation", "auc", "accuracy", "f1"))
  lambda <- result$lambda %||% seq_along(result$validation_correlation)

  if (metric == "correlation") {
    y <- result$validation_correlation
    ylab <- "Validation Pearson Correlation"
  } else if (metric == "auc") {
    stop("AUC is only calculated for best_lambda on test set. No per-lambda AUC available.")
  } else {
    stop("Only one test set value exists for accuracy/F1. Not defined across lambdas.")
  }

  plot(lambda, y, type = "b", pch = 19,
       xlab = "Lambda", ylab = ylab,
       log = if (log_lambda) "x" else "")
  abline(v = result$best_lambda, col = "red", lty = 2)
  legend("bottomright", legend = paste0("Best lambda = ", round(result$best_lambda, 4)),
         col = "red", lty = 2, bty = "n")
}

