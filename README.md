# mrs

An R package for constructing **Multilocus Risk Scores (MRS)** based on SNP–SNP interaction effects. This package extends the original [lassosum](https://github.com/tshmak/lassosum) concept to incorporate pairwise epistasis using MB-MDR output and LASSO-style shrinkage.

## Features

- Support for both **binary** and **continuous** phenotypes  
- Direct parsing of MB-MDR `.model` and `.output` files  
- Efficient shrinkage of F-values via soft-thresholding (`mrs_indeplasso`)  
- HLO matrix parsing to encode pairwise interaction patterns (`parse_model_to_mdr`)  
- Construction of individual-level MRS using genotype data (`compute_one_mrs`, `compute_one_mrs_fast`)  
- Built-in **cross-validation** to select optimal lambda (`mrs_splitvalidate`)  
- Full pipeline wrapper for streamlined analysis (`mrs.pipeline`)
→ integrates shrinkage, interaction encoding, risk construction, and validation in one call  
- Visualization of **quartile-based risk separation** (`plot.mrs_quartile_risk`)

## Installation

To install from GitHub:

```r
install.packages("devtools")  # if not already installed
devtools::install_github("m88020362/mrs")
```

## License

This package is released under the MIT License. See [`LICENSE`](LICENSE) for details.

## Attribution

Portions of this package are based on:

- [tshmak/lassosum](https://github.com/tshmak/lassosum)  
- [lelaboratoire/rethink-prs](https://github.com/lelaboratoire/rethink-prs)

