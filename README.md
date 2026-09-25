# synpdx

`synpdx` is an R package for nonlinear mixed-effects modeling and interaction
analysis in *in vivo* patient-derived xenograft (PDX) combination studies.
It fits tumor-growth and drug-effect models, computes subject-level
interaction distances (D_i), and performs chi-square and cluster-bootstrap
inference for synergism, additivity, or antagonism.

---

## Installation

You can install the development version of synpdx from GitHub using remotes:

```r
remotes::install_github("mothquitoman-kor/synpdx")
```

---

## Input data

This example repository includes a curated sample data file at
`inst/extdata/PDX_CRC_BYL719+binimetinib_curated.csv`. Your own analyses can
use any data frame with one row per observation and these columns:

`subject`, `time`, `treatment`, and `logvol_norm`.

`logvol_norm` is the normalized log-volume, typically `log(V(t)/V(0))`,
calculated within each subject and treatment arm. The four treatment labels
are passed explicitly to `synpdx_fit`.

## Example

```r
library(synpdx)

# Use the included sample data. Replace this path with your own CSV for a new study.
csv_path <- system.file(
  "extdata", "PDX_CRC_BYL719+binimetinib_curated.csv", package = "synpdx"
)
if (!nzchar(csv_path)) csv_path <- file.path("inst", "extdata", "PDX_CRC_BYL719+binimetinib_curated.csv")
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# Quick example run. Auto structural/random-effect selection can be enabled for
# a publication analysis, but it requires substantially more fitting time.
res <- synpdx_fit(
  data = df,
  control = "control",
  drug_a = "BYL719",
  drug_b = "binimetinib",
  combo = "combination",
  B = 0,
  seed = 2025,
  cutoff = "auto",
  model = "exp_constA_constB",
  select_random_effects = FALSE,
  out_dir = "synpdx_out",
  distance = "marginal"
)

# Inspect results
res$fit$selected_model_name   # chosen model
res$chisq$overall             # chi-square test result
res$bootstrap$summary         # bootstrap summary
```

---

## Statistical defaults

The default distance is the manuscript version: pointwise log-scale Jacobian
variance plus the fitted residual variance. Structural model selection is
structure-first, with `re_alpha = 0.01`. Bootstrap refits keep the selected
structural model and use the full diagonal random-effects covariance without
repeating random-effects selection. This bootstrap behavior is intentional.

Set `distance = "residual"` only to reproduce the older residual-MSE
denominator.

## Output files

When `out_dir` is specified, the following files are written:

| Output | Description |
|--------|-------------|
| `individual_prediction_plot.png` | Fitted vs observed tumor growth curves |
| `combo_pred_vs_obs.csv` | Predicted vs observed values for the combination arm |
| `Di_by_subject.csv` | Subject-level interaction distances (Dᵢ) |
| `chi_square_result.csv`, `chi_square_subject_results.csv`, `chi_square_test_plot.png` | Chi-square test outputs |
| `bootstrap_summary_result.csv`, `bootstrap_subject_Si.csv`, `bootstrap_ci_overlay.png` | Bootstrap test outputs |

---

## Citation

The work is not yet published
