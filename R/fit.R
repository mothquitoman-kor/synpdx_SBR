.synpdx_models <- function() {
  nms <- as.vector(outer(c("exp", "gomp", "logistic"),
                        c("constA_constB", "constA_decayB", "decayA_constB", "decayA_decayB"), paste, sep = "_"))
  ans <- lapply(nms, get, envir = environment(exp_constA_constB))
  names(ans) <- nms
  ans
}

get_bic_saemix <- function(obj) {
  val <- as.numeric(stats::BIC(obj, method = "is"))
  if (length(val) != 1L || !is.finite(val)) stop("Finite importance-sampling BIC is unavailable.")
  val
}

.synpdx_ll <- function(obj) {
  val <- obj@results@ll.is
  if (length(val) != 1L || !is.finite(val)) stop("Finite importance-sampling log likelihood is unavailable.")
  as.numeric(val)
}

.synpdx_control <- function(seed, control = list()) {
  .synpdx_integer(seed, "seed")
  if (!is.list(control)) stop("saemix_control must be a uniquely named list.")
  if (length(control) && (is.null(names(control)) || anyNA(names(control)) ||
                          any(!nzchar(names(control))) || anyDuplicated(names(control))))
    stop("saemix_control must be a uniquely named list.")
  valid <- names(formals(saemix::saemixControl))
  if (!all(names(control) %in% valid)) stop("Unknown saemix_control option: ", paste(setdiff(names(control), valid), collapse = ", "))
  opts <- utils::modifyList(list(fim = FALSE, print = FALSE, displayProgress = FALSE,
                                 save = FALSE, save.graphs = FALSE), control)
  opts$seed <- seed; opts$fix.seed <- TRUE; opts$map <- TRUE; opts$ll.is <- TRUE
  do.call(saemix::saemixControl, opts)
}

.synpdx_fit_one <- function(fun, dat, mask, ctrl, psi0 = NULL) {
  pn <- attr(fun, "parnames")
  if (is.null(psi0)) psi0 <- matrix(attr(fun, "theta0"), nrow = 1, dimnames = list(NULL, pn))
  mod <- saemix::saemixModel(model = fun, psi0 = psi0, name.modpar = pn,
                             transform.par = rep(0, length(pn)), covariance.model = mask,
                             error.model = "constant", verbose = FALSE)
  if (!inherits(mod, "SaemixModel")) stop("saemix failed to construct the model.")
  ans <- saemix::saemix(mod, dat, control = ctrl)
  if (!inherits(ans, "SaemixObject")) stop("saemix did not return a fitted object.")
  .synpdx_ll(ans)
  .synpdx_parameters(ans)
  ans
}

backward_selection_saemix <- function(par_names, saemix_data, saemix_model_fun, psi0,
                                     alpha = 0.01, seed = 990202, saemix_control = list(),
                                     initial_fit = NULL) {
  .synpdx_probability(alpha, "alpha")
  p <- length(par_names)
  if (ncol(psi0) != p) stop("psi0 and par_names differ in length.")
  keep <- par_names
  mask_from <- function(k) {
    m <- diag(as.numeric(par_names %in% k), nrow = p)
    dimnames(m) <- list(par_names, par_names)
    m
  }
  ctrl <- .synpdx_control(seed, saemix_control)
  full <- if (is.null(initial_fit)) .synpdx_fit_one(saemix_model_fun, saemix_data, mask_from(keep), ctrl, psi0) else initial_fit
  history <- list(); step <- 0L
  while (length(keep) > 1L) {
    step <- step + 1L
    trials <- vector("list", length(keep)); pv <- rep(NA_real_, length(keep))
    rows <- vector("list", length(keep))
    for (j in seq_along(keep)) {
      fit_j <- tryCatch(.synpdx_fit_one(saemix_model_fun, saemix_data,
                          mask_from(setdiff(keep, keep[j])), ctrl, psi0), error = identity)
      lr <- NA_real_; status <- "ok"
      if (inherits(fit_j, "error")) status <- conditionMessage(fit_j) else {
        lr <- 2 * (.synpdx_ll(full) - .synpdx_ll(fit_j))
        if (lr < -1e-6) status <- "negative_LR: check SAEM/IS stability" else {
          pv[j] <- if (lr <= 0) 1 else 0.5 * stats::pchisq(lr, 1, lower.tail = FALSE)
          trials[[j]] <- fit_j
        }
      }
      rows[[j]] <- data.frame(step = step, removed = keep[j], LR = lr, p_value = pv[j],
                               accepted = FALSE, status = status)
    }
    tab <- do.call(rbind, rows)
    eligible <- which(is.finite(pv) & pv > alpha)
    if (!length(eligible)) { history[[step]] <- tab; break }
    j <- eligible[which.max(pv[eligible])]
    tab$accepted[j] <- TRUE; history[[step]] <- tab
    keep <- setdiff(keep, keep[j])
    full <- trials[[j]]  
  }
  trace <- if (length(history)) do.call(rbind, history) else data.frame()
  limited <- length(keep) == 1L
  if (limited) warning("One random effect retained because saemix requires IIV; it was not tested against a fixed-effects-only model.", call. = FALSE)
  if (nrow(trace) && any(trace$status != "ok")) warning("Some random-effect comparisons were unavailable; inspect selection_trace.", call. = FALSE)
  list(selected = keep, covariance.model = mask_from(keep), fit = full,
       trace = trace, minimum_random_effects_reached = limited)
}

synpdx_model_fit <- function(dat, control, drug_a, drug_b, combo, seed = 990202,
                             model = "auto", select_random_effects = FALSE,
                             selection_order = c("structure_first", "random_first"),
                             re_alpha = 0.01, saemix_control = list()) {
  .synpdx_flag(select_random_effects, "select_random_effects")
  .synpdx_probability(re_alpha, "re_alpha")
  selection_order <- match.arg(selection_order)
  models <- .synpdx_models()
  model <- match.arg(model, c("auto", names(models)))
  dat <- .synpdx_data(dat, control, drug_a, drug_b, combo)
  train <- dat[dat$treatment != combo, , drop = FALSE]
  if (length(unique(train$subject)) < 2L) stop("At least two patient-derived models are required.")
  sd <- saemix::saemixData(name.data = train, name.group = "subject",
                          name.predictors = c("time", "I_A", "I_B"),
                          name.response = "logvol_norm", verbose = FALSE)
  ctrl <- .synpdx_control(seed, saemix_control)
  fit_one <- function(nm, select) {
    f <- models[[nm]]; pn <- attr(f, "parnames")
    mask <- diag(1, length(pn)); dimnames(mask) <- list(pn, pn)
    obj <- .synpdx_fit_one(f, sd, mask, ctrl)
    sel <- NULL
    if (select) {
      sel <- backward_selection_saemix(pn, sd, f, matrix(attr(f, "theta0"), nrow = 1,
                         dimnames = list(NULL, pn)), re_alpha, seed, saemix_control, obj)
      obj <- sel$fit
    }
    list(object = obj, selection = sel, bic = get_bic_saemix(obj))
  }
  candidates <- if (model == "auto") names(models) else model
  refine_candidates <- select_random_effects && (model != "auto" || selection_order == "random_first")
  fits <- vector("list", length(candidates)); names(fits) <- candidates
  table <- data.frame(model = candidates, BIC = NA_real_, status = "ok", stringsAsFactors = FALSE)
  for (i in seq_along(candidates)) {
    ans <- tryCatch(fit_one(candidates[i], refine_candidates), error = identity)
    if (inherits(ans, "error")) table$status[i] <- conditionMessage(ans) else {
      fits[[i]] <- ans; table$BIC[i] <- ans$bic
    }
  }
  valid <- which(is.finite(table$BIC))
  if (!length(valid)) stop("All candidate models failed: ", paste(paste(table$model, table$status, sep = ": "), collapse = "; "))
  if (length(valid) != length(candidates)) warning("Some candidate models failed; inspect model_comparison.", call. = FALSE)
  chosen <- valid[which.min(table$BIC[valid])]
  nm <- candidates[chosen]; best <- fits[[chosen]]
  if (model == "auto" && select_random_effects && selection_order == "structure_first") {
    f <- models[[nm]]; pn <- attr(f, "parnames")
    sel <- backward_selection_saemix(pn, sd, f, matrix(attr(f, "theta0"), nrow = 1,
                 dimnames = list(NULL, pn)), re_alpha, seed, saemix_control, best$object)
    best <- list(object = sel$fit, selection = sel, bic = get_bic_saemix(sel$fit))
  }
  out <- list(saemix = best$object, selected_model = best$object, selected_model_name = nm,
              data_processed = dat, tmax_tbl = .synpdx_cutoffs(dat), model_comparison = table,
              final_BIC = best$bic, selection_trace = if (is.null(best$selection)) data.frame() else best$selection$trace,
              minimum_random_effects_reached = !is.null(best$selection) && best$selection$minimum_random_effects_reached,
              candidate_selections = lapply(fits, function(x) {
                if (is.null(x) || is.null(x$selection)) return(NULL)
                x$selection[setdiff(names(x$selection), "fit")]
              }),
              settings = list(seed = seed, select_random_effects = select_random_effects,
                              selection_order = selection_order, re_alpha = re_alpha,
                              saemix_control = saemix_control))
  class(out) <- c("synpdx_fit", "list")
  out
}

#' @export
print.synpdx_fit <- function(x, ...) {
  cat("synpdx_fit object\n model:", x$selected_model_name, "\n")
  if (!is.null(x$final_BIC)) cat(" BIC (importance sampling):", format(x$final_BIC), "\n")
  if (isTRUE(x$minimum_random_effects_reached)) cat(" Note: at least one random effect retained (saemix restriction).\n")
  invisible(x)
}
