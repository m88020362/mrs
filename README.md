# mrs

An R package for constructing **Multilocus Risk Scores (MRS)** based on SNP–SNP interaction effects. This package extends the original [lassosum](https://github.com/tshmak/lassosum) concept to incorporate pairwise epistasis using MB-MDR output and LASSO-style shrinkage.

## Features

-   Support for both **binary** and **continuous** phenotypes
-   Direct parsing of MB-MDR `.model` and `.output` files
-   Efficient shrinkage of F-values via soft-thresholding (`mrs_indeplasso`)
-   HLO matrix parsing to encode pairwise interaction patterns (`parse_model_to_mdr`)
-   Construction of individual-level MRS using genotype data (`compute_one_mrs`, `compute_one_mrs_fast`)
    *(use `compute_one_mrs` for safe checks; switch to `compute_one_mrs_fast` for speed once inputs are verified)*
-   Built-in **cross-validation** to select optimal lambda (`mrs_splitvalidate`)
-   Full pipeline wrapper for streamlined analysis (`mrs.pipeline`) → integrates shrinkage, interaction encoding, risk construction, and validation in one call
-   Visualization of **quartile-based risk separation** (`plot.mrs_quartile_risk`)

## Demo Data

-   demo_base.txt — 450 samples (rows); 114 SNP genotypes + 1 binary phenotype (columns)

-   demo_test.txt — 50 samples (rows); 114 SNP genotypes + 1 binary phenotype (columns)

-   demo_output.txt — MB-MDR output for 100 SNP pairs (from 114 SNPs); contains F-values and p-values

-   demo_model.txt — MB-MDR output for the same 100 SNP pairs; contains HLO matrices (pairwise interaction encodings)

Note: the phenotype column is named Class in the demo. If yours differs, adjust the code that drops Class.

## Quick Start: compute MRS with the demo data

This walkthrough shows the function order when you already have MB-MDR outputs and genotype matrices.

## Step-by-step workflow

``` r
library(mrs)
library(data.table)
library(dplyr)

# ---- Load demo data ----
f <- function(n) system.file("extdata", n, package = "mrs")
model_demo <- fread(f("demo_model.txt"))
out_demo   <- fread(f("demo_output.txt"))
base_data  <- fread(f("demo_base.txt"),  check.names = FALSE)
test_data  <- fread(f("demo_test.txt"),  check.names = FALSE)

# ---- Parse interaction patterns (HLO) from .model ----
snp_only_base <- base_data[, !(colnames(base_data) %in% "Class"), drop = FALSE]
mdr <- parse_model_to_mdr(model_demo, base_data = snp_only_base, filter_zero = FALSE)
# Output: 'mdr' is a named list of HLO matrices, one per SNP pair

# ---- Align F-values from .output to the parsed model ----
out_aligned <- out_demo %>%
  filter(pair %in% names(mdr)) %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(match(pair, names(mdr))) %>%
  select(pair, F_value)
# Output: data.frame with columns: pair, F_value (same order as 'mdr')

# ---- Shrink F-values (produces lambda path & coefficients) ----
shr <- mrs_indeplasso(out_aligned)
# Output: list(beta, pair_names, lambda)

# ---- Construct risk matrix [N sample × N lambda] ----
risk <- compute_one_mrs_fast(
  snp_data   = test_data[, !(colnames(test_data) %in% "Class"), drop = FALSE],
  mdr_model  = mdr,
  beta       = shr$beta,
  pair_names = shr$pair_names
)
# Output: numeric matrix 'risk' (rows = samples, cols = lambdas)

# ---- Split-validate and select best lambda by AUC ----
res <- mrs_splitvalidate(
  risk_mat  = risk,
  test_data = test_data,
  lambda    = shr$lambda,
  seed      = 42,
  select_by = "auc"
)
# Output: list(best_lambda, test_auc, results_table)
res$best_lambda
res$test_auc
head(res$results_table, 10)
```

## One-stop pipeline

``` r
# One-stop wrapper: shrinkage + interaction encoding + risk construction + CV in one call
res2 <- mrs.pipeline(
  base_output = out_aligned,  # aligned (pair, F_value)
  test_data   = test_data,    # contains 'Class' + SNPs
  raw_model   = model_demo,   # MB-MDR .model
  lambda      = shr$lambda,   # reuse the same lambda grid for exact comparability
  base_data   = snp_only_base,
  select_by   = "auc",
  seed        = 42
)

# Check if both approaches give identical MRS results (they should be the same)
identical(res$results_table, res2$results_table)  # TRUE
```

## Installation

To install from GitHub:

``` r
install.packages("devtools")  # if not already installed
devtools::install_github("m88020362/mrs")
```

## License

This package is released under the MIT License. See [`LICENSE`](LICENSE) for details.

## Attribution

Portions of this package are based on:

-   [tshmak/lassosum](https://github.com/tshmak/lassosum)\
-   [lelaboratoire/rethink-prs](https://github.com/lelaboratoire/rethink-prs)
