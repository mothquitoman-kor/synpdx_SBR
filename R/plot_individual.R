#' Plot observations and MAP predictions with a common analysis cutoff.
#' @export
synpdx_plot_individual <- function(fit, data, control, drug_a, drug_b, combo,
                                   cutoff = c("auto", "none"), out_dir = NULL,
                                   width = 15, height = 12, dpi = 300) {
  cutoff <- match.arg(cutoff)
  dat <- .synpdx_data(data, control, drug_a, drug_b, combo)
  pars <- .synpdx_parameters(fit); cuts <- .synpdx_cutoffs(dat)
  predictions <- lapply(split(dat, dat$subject), function(d) {
    cut <- cuts$tmax[match(d$subject[1], cuts$subject)]
    times <- sort(unique(c(seq(0, max(d$time), by = 0.5), max(d$time), cut)))
    grid <- expand.grid(time = times, curve = c("C", "A", "B", "null"), stringsAsFactors = FALSE)
    grid$subject <- d$subject[1]
    grid$I_A <- as.integer(grid$curve %in% c("A", "null"))
    grid$I_B <- as.integer(grid$curve %in% c("B", "null"))
    grid$predicted_logvol <- .synpdx_predict(fit, grid, pars)
    grid$phase <- if (cutoff == "auto") ifelse(grid$time <= cut, "within", "beyond") else "within"
    if (cutoff == "auto" && max(times) > cut) {
      boundary <- grid[grid$time == cut, ]; boundary$phase <- "beyond"
      grid <- rbind(grid, boundary)
    }
    grid
  })
  pred <- do.call(rbind, predictions)
  dat$curve <- c("C", "A", "B", "AB")[match(dat$treatment, c(control, drug_a, drug_b, combo))]
  colors <- c(C = "#1b9e77", A = "#d95f02", B = "#7570b3", AB = "#333333", null = "black")
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = dat, ggplot2::aes(x = time, y = logvol_norm, color = curve, shape = curve), size = 1.8) +
    ggplot2::geom_line(data = pred, ggplot2::aes(x = time, y = predicted_logvol, color = curve,
                          linetype = phase, group = interaction(subject, curve, phase)), linewidth = 0.7) +
    ggplot2::facet_wrap(~subject, scales = "free_y") +
    ggplot2::scale_color_manual(name = NULL, values = colors,
                                 breaks = c("C", "A", "B", "AB", "null"),
                                 labels = c(control, drug_a, drug_b, combo, "Additive (predicted)")) +
    ggplot2::scale_shape_manual(values = c(C = 19, A = 19, B = 19, AB = 1), guide = "none") +
    ggplot2::scale_linetype_manual(name = NULL, values = c(within = "solid", beyond = "dashed"),
                                    breaks = if (cutoff == "auto") c("within", "beyond") else "within",
                                    labels = if (cutoff == "auto") c("Within cutoff", "Beyond cutoff") else "Prediction") +
    ggplot2::labs(x = "Time (days)", y = "Normalized log-volume") +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
  if (cutoff == "auto") p <- p + ggplot2::geom_vline(data = cuts, ggplot2::aes(xintercept = tmax),
                                                    color = "red", linetype = "dotted", linewidth = 0.6)
  .synpdx_output(out_dir)
  if (!is.null(out_dir)) ggplot2::ggsave(file.path(out_dir, "individual_prediction_plot.png"), p,
                                        width = width, height = height, dpi = dpi)
  attr(p, "synpdx_predictions") <- pred
  p
}
