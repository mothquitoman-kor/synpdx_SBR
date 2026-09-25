# Run this script from the repository root after installing synpdx.
# The sample data are included at inst/extdata/.

if (!requireNamespace("synpdx", quietly = TRUE)) {
  stop("Install the package first, for example with remotes::install_local('.').")
}

csv_path <- system.file(
  "extdata", "PDX_CRC_BYL719+binimetinib_curated.csv", package = "synpdx"
)
if (!nzchar(csv_path)) {
  csv_path <- file.path("inst", "extdata", "PDX_CRC_BYL719+binimetinib_curated.csv")
}
if (!file.exists(csv_path)) stop("Sample data not found: ", csv_path)

dat <- read.csv(csv_path, stringsAsFactors = FALSE)

# B = 0 gives a quick example run. Use B >= 2 for cluster bootstrap inference.
fit <- synpdx::synpdx_fit(
  data = dat,
  control = "control",
  drug_a = "BYL719",
  drug_b = "binimetinib",
  combo = "combination",
  B = 0,
  model = "exp_constA_constB",
  select_random_effects = FALSE,
  distance = "marginal",
  out_dir = "synpdx_sample_results"
)

print(fit$fit)
print(fit$chisq$overall)
print(fit$Di$Di)
