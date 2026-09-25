#' @export
synpdx_compute_Di <- function(fit, data, control, drug_a, drug_b, combo,
                              cutoff = c("auto", "none"), out_dir = NULL,
                              distance = c("marginal", "residual")) {
  cutoff <- match.arg(cutoff); distance <- match.arg(distance)
  dat <- .synpdx_data(data, control, drug_a, drug_b, combo)
  obj <- .synpdx_saem(fit); pars <- .synpdx_parameters(fit)
  train <- dat[dat$treatment != combo, , drop = FALSE]
  if (distance == "residual") {
    sigma2 <- mean((.synpdx_predict(fit, train, pars) - train$logvol_norm)^2)
  } else {
    if (!identical(as.character(obj@model@error.model), "constant")) stop("A constant residual error model is required.")
    sigma2 <- obj@results@respar[1]^2
  }
  if (length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0)
    stop("Residual variance must be finite and positive; no outcome-variance fallback is used.")
  d <- dat[dat$treatment == combo, , drop = FALSE]
  cuts <- .synpdx_cutoffs(dat)
  d$cutoff <- if (cutoff == "auto") cuts$tmax[match(d$subject, cuts$subject)] else NA_real_
  if (cutoff == "auto") d <- d[d$time <= d$cutoff, , drop = FALSE]
  missing <- setdiff(unique(dat$subject), unique(d$subject))
  if (length(missing)) stop("No combination observations remain within cutoff for: ", paste(missing, collapse = ", "))
  d$DV <- d$logvol_norm
  d$predicted_logvol <- .synpdx_predict(fit, d, pars)
  d$variance <- sigma2
  if (distance == "marginal") {
    omega <- obj@results@omega
    p <- ncol(pars)
    if (!identical(dim(omega), c(p, p)) || any(!is.finite(omega)) ||
        !isTRUE(all.equal(omega, t(omega), check.attributes = FALSE)) ||
        min(eigen(omega, symmetric = TRUE, only.values = TRUE)$values) < -1e-10)
      stop("Invalid random-effect covariance matrix.")
    for (sid in unique(d$subject)) {
      rows <- which(d$subject == sid); theta <- pars[sid, ]
      J <- matrix(0, length(rows), p)
      for (k in seq_len(p)) {
        h <- .Machine$double.eps^(1/3) * max(1, abs(theta[k]))
        plus <- minus <- pars[sid, , drop = FALSE]
        plus[1, k] <- theta[k] + h; minus[1, k] <- theta[k] - h
        J[, k] <- (.synpdx_predict(fit, d[rows, ], plus) -
                     .synpdx_predict(fit, d[rows, ], minus)) / (2 * h)
      }
      d$variance[rows] <- sigma2 + pmax(0, rowSums((J %*% omega) * J))
    }
  }
  if (any(!is.finite(d$variance)) || any(d$variance <= 0)) stop("Invalid pointwise variance.")
  d$distance_contribution <- (d$predicted_logvol - d$DV)^2 / d$variance
  if (any(!is.finite(d$distance_contribution))) stop("Non-finite distance contribution.")
  Di_vec <- vapply(split(d$distance_contribution, d$subject), sum, numeric(1))
  combo_df <- d[, c("subject", "time", "DV", "predicted_logvol", "I_A", "I_B", "cutoff",
                    "variance", "distance_contribution")]
  .synpdx_output(out_dir)
  if (!is.null(out_dir)) {
    utils::write.csv(combo_df, file.path(out_dir, "combo_pred_vs_obs.csv"), row.names = FALSE)
    utils::write.csv(data.frame(subject = names(Di_vec), Di = unname(Di_vec)),
                     file.path(out_dir, "Di_by_subject.csv"), row.names = FALSE)
  }
  list(Di_vec = Di_vec, combo_df = combo_df, sigma2_used = sigma2, distance = distance)
}
