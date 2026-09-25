#' Fit control/monotherapy NLME models and assess combination interaction.
#' @export
synpdx_fit <- function(dat = NULL, control = NULL, drug_a = NULL, drug_b = NULL, combo = NULL,
                       B = 100, out_dir = NULL, seed = 990202,
                       cutoff = c("auto", "none"), model = "auto",
                       select_random_effects = TRUE,
                       distance = c("marginal", "residual"),
                       selection_order = c("structure_first", "random_first"),
                       re_alpha = 0.01, alpha = 0.05, saemix_control = list(), verbose = TRUE,
                       data = NULL) {
  if (is.null(dat)) dat <- data
  else if (!is.null(data)) stop("Supply either `dat` or the alias `data`, not both.")
  if (is.null(dat)) stop("Supply a data.frame through `dat` or `data`.")
  if (is.logical(cutoff) && length(cutoff) == 1L && !is.na(cutoff))
    cutoff <- if (cutoff) "auto" else "none"
  cutoff <- match.arg(cutoff); distance <- match.arg(distance)
  selection_order <- match.arg(selection_order)
  .synpdx_integer(B, "B"); .synpdx_integer(seed, "seed")
  .synpdx_flag(verbose, "verbose"); .synpdx_flag(select_random_effects, "select_random_effects")
  .synpdx_probability(alpha, "alpha"); .synpdx_probability(re_alpha, "re_alpha")
  if (B == 1L) stop("B must be 0 (skip bootstrap) or at least 2.")
  dat <- .synpdx_data(dat, control, drug_a, drug_b, combo)
  .synpdx_output(out_dir)
  fit <- synpdx_model_fit(dat, control, drug_a, drug_b, combo, seed = seed, model = model,
                          select_random_effects = select_random_effects, selection_order = selection_order,
                          re_alpha = re_alpha, saemix_control = saemix_control)
  Di <- synpdx_compute_Di(fit, dat, control, drug_a, drug_b, combo, cutoff = cutoff,
                           out_dir = out_dir, distance = distance)
  cs <- synpdx_chisq_test(Di$Di_vec, Di$combo_df, alpha = alpha, out_dir = out_dir)
  bt <- if (B == 0L) NULL else synpdx_bootstrap(fit, dat, control, drug_a, drug_b, combo,
                            B = B, seed = seed, cutoff = cutoff, out_dir = out_dir,
                            distance = distance, alpha = alpha, saemix_control = saemix_control, verbose = verbose)
  p <- synpdx_plot_individual(fit, dat, control, drug_a, drug_b, combo, cutoff = cutoff, out_dir = out_dir)
  if (!is.null(out_dir)) {
    utils::write.csv(fit$model_comparison, file.path(out_dir, "model_comparison.csv"), row.names = FALSE)
    if (nrow(fit$selection_trace)) utils::write.csv(fit$selection_trace,
                                   file.path(out_dir, "random_effect_selection.csv"), row.names = FALSE)
  }
  list(fit = fit, plot = p, Di = Di, chisq = cs, bootstrap = bt,
       settings = list(cutoff = cutoff, distance = distance, alpha = alpha,
                       selection_order = selection_order, re_alpha = re_alpha, B = B, seed = seed))
}
