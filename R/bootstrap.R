#' @export
synpdx_bootstrap <- function(fit, data, control, drug_a, drug_b, combo,
                             B = 100, seed = 202, cutoff = c("auto", "none"), out_dir = NULL,
                             distance = c("marginal", "residual"), alpha = 0.05,
                             saemix_control = NULL, verbose = TRUE) {
  cutoff <- match.arg(cutoff); distance <- match.arg(distance)
  .synpdx_integer(B, "B", 2); .synpdx_integer(seed, "seed")
  .synpdx_probability(alpha, "alpha"); .synpdx_flag(verbose, "verbose")
  if (!inherits(fit, "synpdx_fit")) stop("fit must be a synpdx_fit object.")
  model_name <- fit$selected_model_name
  if (length(model_name) != 1L || !model_name %in% names(.synpdx_models())) stop("Invalid selected_model_name.")
  dat <- .synpdx_data(data, control, drug_a, drug_b, combo)
  subjects <- unique(dat$subject); n <- length(subjects)
  if (n < 2L) stop("Bootstrap needs at least two subjects.")
  if (is.null(saemix_control)) saemix_control <- fit$settings$saemix_control
  if (is.null(saemix_control)) saemix_control <- list()
  .synpdx_control(seed, saemix_control)
  obs <- .synpdx_signed(synpdx_compute_Di(fit, dat, control, drug_a, drug_b, combo,
                                        cutoff = cutoff, distance = distance))
  original <- stats::median(obs$S_norm)
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  draws <- replicate(B, sample.int(n, n, replace = TRUE), simplify = FALSE)
  fit_seeds <- sample.int(.Machine$integer.max, B)
  medians <- rep(NA_real_, B); all_subjects <- vector("list", B)
  failures <- rep(NA_character_, B); warnings_by_rep <- vector("list", B)
  for (b in seq_len(B)) {
    if (verbose) message(sprintf("[bootstrap %s] %d/%d", model_name, b, B))
    boot <- do.call(rbind, lapply(seq_len(n), function(k) {
      d <- dat[dat$subject == subjects[draws[[b]][k]], , drop = FALSE]
      d$subject <- sprintf("boot_%06d", k)
      d
    }))
    ans <- tryCatch(withCallingHandlers({
      fb <- synpdx_model_fit(boot, control, drug_a, drug_b, combo, seed = fit_seeds[b],
                            model = model_name, select_random_effects = FALSE,
                            saemix_control = saemix_control)
      res <- .synpdx_signed(synpdx_compute_Di(fb, boot, control, drug_a, drug_b, combo,
                                             cutoff = cutoff, distance = distance))
      res$replicate <- b; res$model <- model_name
      res$source_subject <- subjects[draws[[b]][match(res$subject, sprintf("boot_%06d", seq_len(n)))]]
      res
    }, warning = function(w) {
      warnings_by_rep[[b]] <<- c(warnings_by_rep[[b]], conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
    if (inherits(ans, "error")) failures[b] <- conditionMessage(ans) else {
      all_subjects[[b]] <- ans; medians[b] <- stats::median(ans$S_norm)
    }
  }
  ok <- is.finite(medians)
  valid <- all(ok)
  lo <- hi <- rep(NA_real_, 2)
  if (valid) {
    q <- stats::quantile(medians, c(alpha/2, 1 - alpha/2), names = FALSE)
    lo <- c(q[1], 2 * original - q[2]); hi <- c(q[2], 2 * original - q[1])
  } else warning(sum(!ok), " bootstrap fits failed. CIs withheld; inspect replicate_results and increase fitting effort or inspect data.", call. = FALSE)
  labels <- if (valid) ifelse(lo > 0, "synergistic", ifelse(hi < 0, "antagonistic", "additive")) else rep("undetermined", 2)
  summary <- data.frame(method = c("percentile", "pivotal"), model = model_name,
                         S_norm_median = original, CI_lower = lo, CI_upper = hi,
                         conf_level = 1 - alpha, overall_interaction = labels,
                         B_requested = B, B_success = sum(ok), inference_valid = valid,
                         random_effects = "full_diagonal_no_selection", distance = distance)
  if (alpha == 0.05) { summary$CI_95_lower <- lo; summary$CI_95_upper <- hi }
  rep_results <- data.frame(replicate = seq_len(B), seed = fit_seeds, median = medians,
                            error = failures, warnings = vapply(warnings_by_rep, paste, character(1), collapse = " | "))
  subject_results <- if (any(ok)) do.call(rbind, all_subjects[ok]) else data.frame()
  plot_data <- data.frame(x = medians[ok])
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x))
  if (nrow(plot_data) >= 2L && length(unique(plot_data$x)) > 1L) plot <- plot + ggplot2::geom_density()
  else if (nrow(plot_data)) plot <- plot + ggplot2::geom_rug()
  plot <- plot + ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = original, linetype = "dashed") +
    ggplot2::labs(x = "Bootstrap median S_norm", y = "Density",
                   subtitle = if (valid) paste("Successful replicates:", B) else "Incomplete bootstrap: confidence intervals withheld") + ggplot2::theme_bw()
  .synpdx_output(out_dir)
  if (!is.null(out_dir)) {
    utils::write.csv(summary, file.path(out_dir, "bootstrap_summary_result.csv"), row.names = FALSE)
    utils::write.csv(subject_results, file.path(out_dir, "bootstrap_subject_Si.csv"), row.names = FALSE)
    utils::write.csv(rep_results, file.path(out_dir, "bootstrap_replicates.csv"), row.names = FALSE)
    ggplot2::ggsave(file.path(out_dir, "bootstrap_ci_overlay.png"), plot, width = 9, height = 7, dpi = 300)
  }
  list(summary = summary, subject_Si = subject_results, plot = plot, replicate_results = rep_results,
       resampling_plan = lapply(draws, function(x) subjects[x]), observed_subjects = obs)
}
