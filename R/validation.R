# Validation, tuning, comparison and recommendation -------------------------

.cca_folds <- function(n, k, seed = NULL, min_assess = 1L) {
  n <- as.integer(n); k <- as.integer(k); min_assess <- as.integer(min_assess)
  if (n < 2L) .cca_stop("At least two observations are required to form folds.")
  if (min_assess < 1L) min_assess <- 1L
  k <- min(k, floor(n / min_assess))
  if (k < 2L) .cca_stop("The requested validation requires at least two folds with enough held-out observations.")
  .cca_with_seed(seed, sample(rep(seq_len(k), length.out = n)))
}

.cca_clean_params <- function(x) {
  x <- as.list(x)
  lapply(x, function(z) if (is.factor(z)) as.character(z) else z)
}

.cca_test_cor <- function(fit, Xtest, Ytest) {
  if (nrow(Xtest) < 3L || nrow(Ytest) < 3L) .cca_stop("At least three held-out paired observations are required to estimate an out-of-sample correlation.")
  pr <- predict(fit, newdata_x = Xtest, newdata_y = Ytest)
  rho <- .cca_cor_cols(as.matrix(pr$x), as.matrix(pr$y))
  data.frame(component = seq_along(rho), correlation = rho,
             absolute_correlation = abs(rho))
}

#' Repeated k-fold out-of-sample validation
#' @export
cca_validate <- function(X, Y, method = "classic", params = list(), folds = 5,
                         repeats = 5, seed = 1) {
  if (identical(method, "deep") && !is.null(params$model) && is.null(params$model_factory) && is.null(params$fit_fun)) {
    .cca_stop("Deep CCA resampling must use model_factory or a stateless fit_fun so that each fold starts from a fresh model.")
  }
  pair <- .cca_pair(X, Y, finite = TRUE); X <- pair$X; Y <- pair$Y
  folds <- min(as.integer(folds), floor(nrow(X) / 3L))
  if (folds < 2L) .cca_stop("Repeated k-fold validation requires enough observations for at least two assessment folds of size >= 3.")
  results <- list(); pos <- 1L
  for (r in seq_len(repeats)) {
    foldid <- .cca_folds(nrow(X), folds, if (is.null(seed)) NULL else seed + r - 1L, min_assess = 3L)
    for (f in seq_len(max(foldid))) {
      te <- which(foldid == f); tr <- setdiff(seq_len(nrow(X)), te)
      fit <- tryCatch(.cca_fit_method(method, X[tr,,drop=FALSE], Y[tr,,drop=FALSE], params), error = function(e) NULL)
      if (is.null(fit)) next
      rr <- tryCatch(.cca_test_cor(fit, X[te,,drop=FALSE], Y[te,,drop=FALSE]), error = function(e) NULL)
      if (is.null(rr)) next
      rr[["repeat"]] <- r; rr$fold <- f; rr$n_train <- length(tr); rr$n_test <- length(te)
      results[[pos]] <- rr; pos <- pos + 1L
    }
  }
  tab <- if (length(results)) do.call(rbind, results) else data.frame()
  out <- list(method = method, params = params, results = tab, folds = folds,
              repeats = repeats, seed = seed, X = X, Y = Y)
  class(out) <- "cca_validation"
  out
}

#' Rolling-origin validation for time-ordered data
#' @export
cca_validate_time <- function(X, Y, method = "ridge", params = list(),
                              initial, assess = 5, skip = 0, cumulative = TRUE) {
  pair <- .cca_pair(X, Y, finite = TRUE); X <- pair$X; Y <- pair$Y; n <- nrow(X)
  initial <- as.integer(initial); assess <- as.integer(assess); skip <- as.integer(skip)
  if (assess < 3L) .cca_stop("assess must contain at least three held-out observations so that correlation is estimable.")
  if (initial < 3 || initial + assess > n) .cca_stop("initial/assess leave no valid rolling split.")
  starts <- seq.int(initial + 1L, n - assess + 1L, by = assess + skip)
  res <- list(); splits <- list(); pos <- 1L
  for (s in starts) {
    tr <- if (cumulative) seq_len(s - 1L) else seq.int(max(1L, s - initial), s - 1L)
    te <- seq.int(s, min(n, s + assess - 1L))
    ff <- tryCatch(.cca_fit_method(method, X[tr,,drop=FALSE], Y[tr,,drop=FALSE], params), error = function(e) NULL)
    if (is.null(ff)) next
    rr <- tryCatch(.cca_test_cor(ff, X[te,,drop=FALSE], Y[te,,drop=FALSE]), error = function(e) NULL)
    if (is.null(rr)) next
    rr$split <- pos; rr$train_start <- min(tr); rr$train_end <- max(tr)
    rr$test_start <- min(te); rr$test_end <- max(te)
    res[[pos]] <- rr; splits[[pos]] <- list(train = tr, test = te); pos <- pos + 1L
  }
  out <- list(method = method, params = params,
              results = if (length(res)) do.call(rbind, res) else data.frame(), splits = splits,
              initial = initial, assess = assess, skip = skip, cumulative = cumulative, X = X, Y = Y)
  class(out) <- c("cca_validation_time", "cca_validation")
  out
}

#' Group-preserving cross-validation
#' @export
cca_validate_grouped <- function(X, Y, groups, method = "ridge", params = list(),
                                 folds = 5, repeats = 3, seed = 1) {
  pair <- .cca_pair(X, Y, finite = TRUE); X <- pair$X; Y <- pair$Y
  if (length(groups) != nrow(X)) .cca_stop("groups must have one value per row.")
  ug <- unique(groups); folds <- min(as.integer(folds), length(ug), floor(nrow(X) / 3L))
  if (folds < 2L) .cca_stop("Grouped validation requires at least two folds and enough rows for assessment correlations.")
  res <- list(); splits <- list(); pos <- 1L
  for (r in seq_len(repeats)) {
    gf <- .cca_folds(length(ug), folds, if (is.null(seed)) NULL else seed + r - 1L, min_assess = 1L)
    names(gf) <- as.character(ug)
    rowfold <- gf[match(as.character(groups), names(gf))]
    for (f in seq_len(folds)) {
      te <- which(rowfold == f); tr <- which(rowfold != f)
      ff <- tryCatch(.cca_fit_method(method, X[tr,,drop=FALSE], Y[tr,,drop=FALSE], params), error = function(e) NULL)
      if (is.null(ff)) next
      rr <- tryCatch(.cca_test_cor(ff, X[te,,drop=FALSE], Y[te,,drop=FALSE]), error = function(e) NULL)
      if (is.null(rr)) next
      rr[["repeat"]] <- r; rr$fold <- f; rr$n_groups_test <- length(unique(groups[te]))
      res[[pos]] <- rr
      splits[[pos]] <- list(train = tr, test = te, train_groups = unique(groups[tr]), test_groups = unique(groups[te]))
      pos <- pos + 1L
    }
  }
  out <- list(method = method, params = params, groups = groups,
              results = if (length(res)) do.call(rbind, res) else data.frame(), splits = splits,
              folds = folds, repeats = repeats, seed = seed, X = X, Y = Y)
  class(out) <- c("cca_validation_grouped", "cca_validation")
  out
}

#' Hyperparameter grid search using held-out canonical association
#' @export
cca_tune <- function(X, Y, method, grid, fixed_params = list(), folds = 5,
                     repeats = 3, seed = 1, objective = c("first", "mean")) {
  objective <- match.arg(objective)
  if (!is.data.frame(grid) || !nrow(grid)) .cca_stop("grid must be a non-empty data.frame.")
  rows <- vector("list", nrow(grid)); fits <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    pars <- c(.cca_clean_params(grid[i, , drop = FALSE]), fixed_params)
    # All candidate rows use identical resampling partitions for a fair comparison.
    cv <- tryCatch(cca_validate(X, Y, method, pars, folds, repeats, seed), error = function(e) NULL)
    fits[[i]] <- cv
    if (is.null(cv) || !nrow(cv$results)) {
      score <- NA_real_; sdev <- NA_real_
    } else {
      z <- cv$results
      if (objective == "first") z <- z[z$component == 1, , drop = FALSE]
      vals <- z$correlation; score <- mean(vals, na.rm = TRUE); sdev <- stats::sd(vals, na.rm = TRUE)
    }
    rows[[i]] <- cbind(grid[i, , drop = FALSE], mean_validation = score, sd_validation = sdev)
  }
  tab <- do.call(rbind, rows)
  best <- if (all(!is.finite(tab$mean_validation))) NA_integer_ else which.max(tab$mean_validation)
  out <- list(method = method, grid = grid, results = tab, best_index = best,
              best_params = if (is.na(best)) NULL else .cca_clean_params(grid[best,,drop=FALSE]),
              validations = fits, objective = objective)
  class(out) <- "cca_tuning"
  out
}

#' Nested cross-validation for tuned CCA models
#' @export
cca_nested_cv <- function(X, Y, method, grid, fixed_params = list(), outer_folds = 5,
                          inner_folds = 4, inner_repeats = 2, seed = 1,
                          objective = c("first", "mean")) {
  objective <- match.arg(objective)
  pair <- .cca_pair(X, Y, finite = TRUE); X <- pair$X; Y <- pair$Y
  outer_folds <- min(as.integer(outer_folds), floor(nrow(X) / 3L))
  if (outer_folds < 2L) .cca_stop("Nested CV requires enough observations for at least two outer assessment folds of size >= 3.")
  outer <- .cca_folds(nrow(X), outer_folds, seed, min_assess = 3L)
  res <- list(); chosen <- list(); pos <- 1L
  for (f in seq_len(max(outer))) {
    te <- which(outer == f); tr <- which(outer != f)
    tune <- cca_tune(X[tr,,drop=FALSE], Y[tr,,drop=FALSE], method, grid, fixed_params,
                     folds = inner_folds, repeats = inner_repeats, seed = seed + 1000L * f,
                     objective = objective)
    if (is.null(tune$best_params)) next
    pars <- c(tune$best_params, fixed_params)
    ff <- tryCatch(.cca_fit_method(method, X[tr,,drop=FALSE], Y[tr,,drop=FALSE], pars), error = function(e) NULL)
    if (is.null(ff)) next
    rr <- tryCatch(.cca_test_cor(ff, X[te,,drop=FALSE], Y[te,,drop=FALSE]), error = function(e) NULL)
    if (is.null(rr)) next
    rr$outer_fold <- f; res[[pos]] <- rr; chosen[[pos]] <- tune$best_params; pos <- pos + 1L
  }
  tab <- if (length(res)) do.call(rbind, res) else data.frame()
  if (nrow(tab)) tab$test_correlation <- tab$correlation
  out <- list(method = method, results = tab, test_correlation = tab$correlation %||% numeric(),
              chosen_params = chosen, grid = grid, outer_folds = outer_folds,
              inner_folds = inner_folds, inner_repeats = inner_repeats, seed = seed)
  class(out) <- "cca_nested_cv"
  out
}

#' Compare several CCA estimators on common criteria
#' @export
cca_compare <- function(X, Y,
                        methods = c("classic", "ridge", "shrinkage", "sparse", "copula", "kernel"),
                        specs = list(), validate = TRUE, folds = 5, repeats = 3, seed = 1) {
  rows <- list(); models <- list(); validations <- list(); pos <- 1L
  for (mth in methods) {
    pars <- specs[[mth]] %||% list()
    ff <- tryCatch(.cca_fit_method(mth, X, Y, pars), error = function(e) structure(list(error = conditionMessage(e)), class = "cca_fit_error"))
    models[[mth]] <- ff
    if (inherits(ff, "cca_fit_error")) {
      rows[[pos]] <- data.frame(method = mth, training_rho1 = NA_real_, validation_rho1 = NA_real_,
                                validation_abs_rho1 = NA_real_, validation_sd = NA_real_,
                                selected_x = NA_integer_, selected_y = NA_integer_,
                                status = ff$error); pos <- pos + 1L; next
    }
    cv <- if (validate) tryCatch(cca_validate(X, Y, mth, pars, folds, repeats, seed), error = function(e) NULL) else NULL
    validations[[mth]] <- cv
    vv <- if (!is.null(cv) && nrow(cv$results)) cv$results$correlation[cv$results$component == 1] else NA_real_
    sx <- if (!is.null(ff$weights_x) && nrow(ff$weights_x) == ncol(.cca_as_matrix(X))) sum(abs(ff$weights_x[,1]) > 1e-8) else NA_integer_
    sy <- if (!is.null(ff$weights_y) && nrow(ff$weights_y) == ncol(.cca_as_matrix(Y))) sum(abs(ff$weights_y[,1]) > 1e-8) else NA_integer_
    rows[[pos]] <- data.frame(method = mth, training_rho1 = ff$correlations[1],
                              validation_rho1 = mean(vv, na.rm = TRUE),
                              validation_abs_rho1 = mean(abs(vv), na.rm = TRUE),
                              validation_sd = stats::sd(vv, na.rm = TRUE),
                              selected_x = sx, selected_y = sy, status = "ok")
    pos <- pos + 1L
  }
  tab <- do.call(rbind, rows)
  out <- list(table = tab, models = models, validations = validations, specs = specs)
  class(out) <- "cca_comparison"
  out
}

#' Recommend CCA method families from data characteristics
#' @export
cca_recommend <- function(X, Y, nonlinear = FALSE, mixed = FALSE, functional = FALSE,
                          multiblock = FALSE, outliers = NULL) {
  au <- cca_audit(X, Y)
  n <- au$n; p <- ncol(.cca_as_matrix(X)); q <- ncol(.cca_as_matrix(Y))
  conds <- au$summary$condition_number
  high_col <- any(is.finite(conds) & conds > 100)
  high_dim <- max(p, q) >= n || (p + q) / n > 0.5
  if (is.null(outliers)) {
    Z <- cbind(.cca_as_matrix(X), .cca_as_matrix(Y))
    med <- apply(Z, 2, stats::median, na.rm = TRUE); md <- apply(Z, 2, stats::mad, na.rm = TRUE)
    md[md == 0 | !is.finite(md)] <- 1
    rz <- abs(sweep(sweep(Z, 2, med, "-"), 2, md, "/"))
    outliers <- mean(apply(rz, 1, max, na.rm = TRUE) > 5) > 0.02
  }
  rec <- data.frame(method = character(), priority = character(), rationale = character())
  add <- function(method, priority, rationale) rec <<- rbind(rec, data.frame(method, priority, rationale))
  if (multiblock) add("RGCCA/SGCCA or multiCCA", "high", "Three or more blocks require an explicit multiview formulation.")
  if (functional) add("Functional CCA / multiCCA functional", "high", "Curves or repeated trajectories should preserve functional structure.")
  if (mixed) add("mixedCCA", "high", "Mixed continuous/binary/truncated data benefit from latent semiparametric association modeling.")
  if (nonlinear) add("Kernel CCA", "high", "A scientifically plausible nonlinear cross-block relation was specified.")
  if (high_dim) add("Ridge/Shrinkage or Sparse CCA", "high", "Dimensionality is large relative to sample size; unregularized sample CCA can overfit.")
  if (high_col) add("Ridge/Shrinkage CCA", "high", "Severe within-block conditioning indicates unstable unregularized directions.")
  if (outliers) add("Robust CCA / Robust Sparse CCA", "high", "Extreme observations may materially alter covariance-based directions.")
  if (!nrow(rec)) add("Classical CCA", "primary", "The dimensionality and conditioning audit does not by itself require regularization or resistance.")
  add("Out-of-sample validation", "required for model comparison", "Choose among competing estimators by generalization and stability, not training correlation alone.")
  list(audit = au, recommendations = rec)
}

#' @export
print.cca_validation <- function(x, ...) {
  cat("CCAflow validation\n  Method:", x$method, "\n")
  if (nrow(x$results)) print(aggregate(correlation ~ component, x$results, function(z) c(mean=mean(z,na.rm=TRUE), sd=stats::sd(z,na.rm=TRUE))))
  else cat("  No successful validation splits.\n")
  invisible(x)
}

#' @export
print.cca_tuning <- function(x, ...) {
  cat("CCAflow tuning:", x$method, "\n"); print(x$results, row.names = FALSE)
  if (!is.null(x$best_params)) cat("Best row:", x$best_index, "\n")
  invisible(x)
}

#' @export
print.cca_nested_cv <- function(x, ...) {
  cat("CCAflow nested cross-validation:", x$method, "\n")
  if (nrow(x$results)) print(aggregate(correlation ~ component, x$results, mean, na.rm = TRUE), row.names = FALSE)
  invisible(x)
}

#' @export
print.cca_comparison <- function(x, ...) { print(x$table, row.names = FALSE); invisible(x) }
