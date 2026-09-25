#' @export
synpdx_chisq_test <- function(Di_vec, combo_df, alpha = 0.05, out_dir = NULL) {
  .synpdx_probability(alpha, "alpha")
  need <- c("subject", "time", "DV", "predicted_logvol", "cutoff")
  if (!is.data.frame(combo_df) || !nrow(combo_df) || !all(need %in% names(combo_df))) stop("Invalid combo_df.")
  if (!is.numeric(Di_vec) || !length(Di_vec) || is.null(names(Di_vec)) || anyDuplicated(names(Di_vec)) ||
      anyNA(names(Di_vec)) || any(!nzchar(names(Di_vec))) || any(!is.finite(Di_vec)) || any(Di_vec < 0))
    stop("Di_vec must be finite, non-negative and uniquely named by subject.")
  combo_df$subject <- as.character(combo_df$subject)
  if (anyNA(combo_df$subject) || any(!nzchar(combo_df$subject)) ||
      !setequal(names(Di_vec), unique(combo_df$subject))) stop("Subject IDs in Di_vec and combo_df must match exactly.")
  if (anyDuplicated(combo_df[c("subject", "time")])) stop("Duplicate combination observations.")
  for (nm in c("time", "DV", "predicted_logvol")) {
    if (!is.numeric(combo_df[[nm]]) || any(!is.finite(combo_df[[nm]]))) stop("Non-finite or non-numeric ", nm, ".")
  }
  s <- .synpdx_signed(list(Di_vec = Di_vec, combo_df = combo_df))
  result <- data.frame(subject = s$subject, n_obs = s$n_i, Di = s$Di,
                       Si = sign(s$mean_resid) * sqrt(s$Di), mean_resid = s$mean_resid,
                       PIT = stats::pchisq(s$Di, s$n_i),
                       chisq_threshold = stats::qchisq(1 - alpha, s$n_i))
  result$significant <- result$Di > result$chisq_threshold
  result$interaction <- ifelse(result$significant & result$mean_resid > 0, "synergistic",
                                ifelse(result$significant & result$mean_resid < 0, "antagonistic", "additive"))
  total <- sum(s$Di); degrees <- sum(s$n_i)
  pvalue <- stats::pchisq(total, degrees, lower.tail = FALSE)
  syn <- sum(result$interaction == "synergistic")
  ant <- sum(result$interaction == "antagonistic")
  overall <- if (pvalue >= alpha) "additive" else if (syn == ant) "undetermined" else if (syn > ant) "synergistic" else "antagonistic"
  pit <- pmin(pmax(result$PIT, 1e-12), 1 - 1e-12)
  ad_stat <- ad_p <- NA_real_
  if (length(pit) >= 3L && length(unique(pit)) >= 3L) {
    if (requireNamespace("goftest", quietly = TRUE)) {
      ad <- goftest::ad.test(pit, null = "punif")
      ad_stat <- as.numeric(ad$statistic); ad_p <- ad$p.value
    } else warning("Install goftest to compute the optional PIT diagnostic.", call. = FALSE)
  }
  summary <- data.frame(chi2_statistic = total, df = degrees, p_value = pvalue,
                        overall_interaction = overall, AD_statistic_goodness_of_fit = ad_stat,
                        AD_p_value_goodness_of_fit = ad_p, alpha = alpha)
  x <- seq(0, max(total, stats::qchisq(0.999, degrees)), length.out = 500)
  curve <- data.frame(x = x, y = stats::dchisq(x, degrees)); curve <- curve[is.finite(curve$y), ]
  plot <- ggplot2::ggplot(curve, ggplot2::aes(x = x, y = y)) + ggplot2::geom_line() +
    ggplot2::geom_area(data = curve[curve$x >= total, ], alpha = 0.3) +
    ggplot2::geom_vline(xintercept = total, linetype = "dashed") +
    ggplot2::labs(title = "Asymptotic chi-square test for group-level interaction",
                   subtitle = paste0("df = ", degrees, "; p = ", format.pval(pvalue)),
                   x = "Chi-square statistic", y = "Density") + ggplot2::theme_bw()
  .synpdx_output(out_dir)
  if (!is.null(out_dir)) {
    utils::write.csv(summary, file.path(out_dir, "chi_square_result.csv"), row.names = FALSE)
    utils::write.csv(result, file.path(out_dir, "chi_square_subject_results.csv"), row.names = FALSE)
    ggplot2::ggsave(file.path(out_dir, "chi_square_test_plot.png"), plot, width = 8, height = 6, dpi = 300)
  }
  list(overall = summary, per_subject = result, plot = plot)
}
