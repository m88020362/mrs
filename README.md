# lassosumMRS

This R package extends the [lassosum](https://github.com/tshmak/lassosum) framework to support **Multilocus Risk Score (MRS)** modeling with SNP–SNP interactions, based on MB-MDR outputs. It is designed for researchers interested in epistasis detection and risk prediction in complex disease genetics.

## Features

-   Import MB-MDR outputs (`.output.txt`, `.model.txt`)
-   Parse HLO encodings into interpretable 3×3 matrices
-   Perform LASSO-style shrinkage on F-values
-   Compute multilocus risk scores (MRS)
-   Visualize ROC, risk score distributions, and validation performance
-   Flexible pipeline function for one-step modeling

## Installation

Install from GitHub install.packages("devtools") devtools::install_github("m88020362/lassosumMRS")

## Demo Data

This package includes example files in `inst/extdata/`:
-   `test_data.txt`: Simulated individual-level data
-   `model_2D.txt`: MB-MDR raw model output
-   `output_2D.txt`: MB-MDR F-statistic output

## Example

```{r}
library(lassosumMRS)

# Load demo MB-MDR outputs
f_test_data <- system.file("extdata", "test_data.txt", package = "lassosumMRS")
f_model <- system.file("extdata", "model_2D.txt", package = "lassosumMRS")
f_output <- system.file("extdata", "output_2D.txt", package = "lassosumMRS")

test_data <- data.table::fread(f_test_data)
raw_model <- data.table::fread(f_model, fill = TRUE, header = FALSE)
output_2D <- data.table::fread(f_output, header = FALSE, skip = 3, fill = FALSE,
  col.names = c("snp1", "snp2", "F_value", "p_value")) %>%
  dplyr::mutate(pair = paste(snp1, snp2, sep = "_"))

# One-step modeling
res <- mrs.pipeline(
  base_output = output_2D,
  test_data = test_data,
  raw_model = raw_model,
  seed = 42
)

# Visualizations
plot_mrs_roc(res)
plot_mrs_score_dist(res)
plot_mrs_validate(res, metric = "correlation")
```

## License
This package is released under the MIT License. See LICENSE for details.

## Attribution
Portions of this package are based on original code by [tshmak/lassosum](https://github.com/tshmak/lassosum) (c) 2016, released under the MIT license.

