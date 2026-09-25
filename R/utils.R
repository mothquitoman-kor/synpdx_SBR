utils::globalVariables(c("x", "y", "time", "logvol_norm", "curve", "predicted_logvol",
                         "phase", "subject", "tmax"))
.synpdx_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop(name, " must be a number strictly between 0 and 1.")
}

.synpdx_integer <- function(x, name, min = 0) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < min || x != floor(x) || x > .Machine$integer.max)
    stop(name, " must be an integer >= ", min, " within R's integer range.")
}

.synpdx_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(name, " must be TRUE or FALSE.")
}

.synpdx_data <- function(dat, control, drug_a, drug_b, combo) {
  arms <- c(control, drug_a, drug_b, combo)
  if (length(arms) != 4L || !is.character(arms) || anyNA(arms) ||
      any(!nzchar(arms)) || anyDuplicated(arms)) stop("Supply four distinct, non-empty arm labels.")
  if (!is.data.frame(dat) || !nrow(dat)) stop("dat must be a non-empty data.frame.")
  need <- c("subject", "time", "treatment", "logvol_norm")
  if (!all(need %in% names(dat))) stop("Missing columns: ", paste(setdiff(need, names(dat)), collapse = ", "))
  dat <- as.data.frame(dat[, need])
  dat$subject <- as.character(dat$subject)
  dat$treatment <- as.character(dat$treatment)
  if (anyNA(dat$subject) || any(!nzchar(trimws(dat$subject))) || anyNA(dat$treatment))
    stop("subject and treatment cannot be missing or blank.")
  dat <- dat[dat$treatment %in% arms, , drop = FALSE]
  if (!nrow(dat)) stop("No observations for the requested arms.")
  if (!is.numeric(dat$time) || !is.numeric(dat$logvol_norm) ||
      any(!is.finite(dat$time)) || any(!is.finite(dat$logvol_norm)) || any(dat$time < 0))
    stop("time and logvol_norm must be finite numeric values; time must be non-negative.")
  if (anyDuplicated(dat[c("subject", "treatment", "time")]))
    stop("Duplicate subject/treatment/time rows: this implementation requires one mouse per arm (1x1x1).")
  complete <- vapply(split(dat$treatment, dat$subject), function(x) all(arms %in% x), logical(1))
  if (any(!complete)) stop("All four arms are required for each subject. Incomplete: ",
                           paste(names(complete)[!complete], collapse = ", "))
  dat <- dat[order(dat$subject, dat$time, match(dat$treatment, arms)), , drop = FALSE]
  rownames(dat) <- NULL
  dat$I_A <- as.integer(dat$treatment %in% c(drug_a, combo))
  dat$I_B <- as.integer(dat$treatment %in% c(drug_b, combo))
  dat
}

.synpdx_cutoffs <- function(dat) {
  do.call(rbind, lapply(split(dat, dat$subject), function(d) {
    data.frame(subject = d$subject[1], tmax = min(tapply(d$time, d$treatment, max)),
               stringsAsFactors = FALSE)
  }))
}

.synpdx_saem <- function(fit) {
  if (inherits(fit, "SaemixObject")) return(fit)
  if (inherits(fit, "synpdx_fit") && inherits(fit$saemix, "SaemixObject")) return(fit$saemix)
  stop("Expected a synpdx_fit or SaemixObject.")
}

.synpdx_parameters <- function(fit) {
  obj <- .synpdx_saem(fit)
  pnames <- obj@model@name.modpar
  if (any(obj@model@transform.par != 0)) stop("Expected internally log-parameterized synpdx models.")
  # psi(type='mode') has exactly the named parameter columns, not an ID column.
  pars <- as.data.frame(saemix::psi(obj, type = "mode"))
  ids <- unique(as.character(obj@data@data[[obj@data@name.group]]))
  if (!all(pnames %in% names(pars)) || nrow(pars) != length(ids))
    stop("Cannot align fitted MAP parameters with fitted subject IDs.")
  pars <- as.matrix(pars[, pnames, drop = FALSE])
  if (any(!is.finite(pars))) stop("Non-finite individual parameter estimates.")
  rownames(pars) <- ids
  pars
}

.synpdx_predict <- function(fit, dat, pars = .synpdx_parameters(fit)) {
  obj <- .synpdx_saem(fit)
  id <- match(as.character(dat$subject), rownames(pars))
  if (anyNA(id)) stop("Prediction contains subjects absent from the fitted model.")
  y <- as.numeric(obj@model@model(pars, id, dat[, c("time", "I_A", "I_B"), drop = FALSE]))
  if (length(y) != nrow(dat) || any(!is.finite(y))) stop("Non-finite or incorrectly sized model predictions.")
  y
}

.synpdx_output <- function(out_dir) {
  if (!is.null(out_dir)) {
    if (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir))
      stop("out_dir must be NULL or a non-empty path.")
    if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) stop("Cannot create out_dir.")
  }
}

.synpdx_signed <- function(Di) {
  d <- Di$combo_df
  ids <- names(Di$Di_vec)
  n <- vapply(ids, function(s) sum(d$subject == s), integer(1))
  avg <- vapply(ids, function(s) mean((d$predicted_logvol - d$DV)[d$subject == s]), numeric(1))
  data.frame(subject = ids, n_i = n, Di = unname(Di$Di_vec), mean_resid = avg,
             S_norm = sign(avg) * unname(Di$Di_vec) / n, row.names = NULL)
}
