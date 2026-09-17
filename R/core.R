# Core CCA ------------------------------------------------------------------

#' Prepare paired data blocks for CCA
#'
#' Creates an explicit preprocessing recipe, with transparent missing-data,
#' centering and scaling decisions.
#' @export
cca_prepare <- function(X, Y, center = TRUE, scale = TRUE,
                        na_action = c("fail", "complete")) {
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  class(prep) <- "cca_prepared"
  prep
}

#' Unified two-block CCA dispatcher
#' @export
cca_fit <- function(X, Y,
                    method = c("classic", "ridge", "shrinkage", "sparse", "elastic",
                               "group_sparse", "robust", "mcd", "rank", "copula",
                               "kernel", "robust_sparse", "robust_kernel", "probabilistic",
                               "bayes", "mixed", "deep"),
                    ...) {
  method <- match.arg(method)
  .cca_fit_method(method, X, Y, list(...))
}

#' Classical canonical correlation analysis
#' @export
cca_classic <- function(X, Y, ncomp = NULL, center = TRUE, scale = TRUE,
                        na_action = c("fail", "complete"), tol = 1e-10) {
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  n <- nrow(prep$Xs)
  Sxx <- crossprod(prep$Xs) / (n - 1)
  Syy <- crossprod(prep$Ys) / (n - 1)
  Sxy <- crossprod(prep$Xs, prep$Ys) / (n - 1)
  sol <- .cca_cov_solution(Sxx, Syy, Sxy, ncomp, tol)
  .cca_make_fit(
    "classic", prep, sol$xcoef, sol$ycoef,
    engine = "internal-whitened-svd", call = match.call(),
    extras = list(cov_x = Sxx, cov_y = Syy, cov_xy = Sxy),
    tuning = list(ncomp = ncol(sol$xcoef), center = center, scale = scale,
                  na_action = match.arg(na_action), tol = tol),
    inference_scope = "classical"
  ) -> fit
  fit$refit_params <- fit$tuning
  fit
}

#' Ridge-regularized CCA
#' @export
cca_ridge <- function(X, Y, lambda_x = 0.1, lambda_y = 0.1, ncomp = NULL,
                      center = TRUE, scale = TRUE,
                      na_action = c("fail", "complete"), tol = 1e-10) {
  if (lambda_x < 0 || lambda_y < 0) .cca_stop("Ridge penalties must be non-negative.")
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  n <- nrow(prep$Xs)
  Sxx0 <- crossprod(prep$Xs) / (n - 1)
  Syy0 <- crossprod(prep$Ys) / (n - 1)
  Sxy <- crossprod(prep$Xs, prep$Ys) / (n - 1)
  Sxx <- Sxx0 + lambda_x * diag(ncol(prep$Xs))
  Syy <- Syy0 + lambda_y * diag(ncol(prep$Ys))
  sol <- .cca_cov_solution(Sxx, Syy, Sxy, ncomp, tol)
  fit <- .cca_make_fit(
    "ridge", prep, sol$xcoef, sol$ycoef, engine = "internal-ridge-svd",
    tuning = list(lambda_x = lambda_x, lambda_y = lambda_y,
                  ncomp = ncol(sol$xcoef), center = center, scale = scale,
                  na_action = match.arg(na_action), tol = tol),
    call = match.call(),
    extras = list(cov_x = Sxx0, cov_y = Syy0, cov_xy = Sxy,
                  regularized_cov_x = Sxx, regularized_cov_y = Syy),
    inference_scope = "resampling"
  )
  fit$refit_params <- fit$tuning
  fit
}

#' Shrinkage canonical correlation analysis
#' @export
cca_shrinkage <- function(X, Y, alpha_x = NULL, alpha_y = NULL,
                          ncomp = NULL, engine = c("internal", "corpcor", "auto"),
                          center = TRUE, scale = TRUE,
                          na_action = c("fail", "complete"), tol = 1e-10) {
  engine <- match.arg(engine)
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  n <- nrow(prep$Xs)
  Sxy <- crossprod(prep$Xs, prep$Ys) / (n - 1)
  if (engine == "auto") engine <- if (requireNamespace("corpcor", quietly = TRUE)) "corpcor" else "internal"
  if (engine == "corpcor") {
    .cca_require("corpcor", "cca_shrinkage(engine='corpcor')")
    Sxx <- corpcor::cov.shrink(prep$Xs, verbose = FALSE)
    Syy <- corpcor::cov.shrink(prep$Ys, verbose = FALSE)
    ax <- attr(Sxx, "lambda") %||% NA_real_
    ay <- attr(Syy, "lambda") %||% NA_real_
  } else {
    Sxx0 <- crossprod(prep$Xs) / (n - 1)
    Syy0 <- crossprod(prep$Ys) / (n - 1)
    if (is.null(alpha_x)) alpha_x <- min(1, ncol(prep$Xs) / max(n - 1, 1)) * 0.1
    if (is.null(alpha_y)) alpha_y <- min(1, ncol(prep$Ys) / max(n - 1, 1)) * 0.1
    if (alpha_x < 0 || alpha_x > 1 || alpha_y < 0 || alpha_y > 1) .cca_stop("alpha_x and alpha_y must lie in [0, 1].")
    Sxx <- (1 - alpha_x) * Sxx0 + alpha_x * diag(diag(Sxx0))
    Syy <- (1 - alpha_y) * Syy0 + alpha_y * diag(diag(Syy0))
    ax <- alpha_x; ay <- alpha_y
  }
  sol <- .cca_cov_solution(Sxx, Syy, Sxy, ncomp, tol)
  fit <- .cca_make_fit(
    "shrinkage", prep, sol$xcoef, sol$ycoef,
    engine = paste0(engine, "-shrinkage"),
    tuning = list(alpha_x = ax, alpha_y = ay, ncomp = ncol(sol$xcoef),
                  engine = engine, center = center, scale = scale,
                  na_action = match.arg(na_action), tol = tol),
    call = match.call(), extras = list(shrink_cov_x = Sxx, shrink_cov_y = Syy, cov_xy = Sxy),
    inference_scope = "resampling"
  )
  fit$refit_params <- fit$tuning
  fit
}

#' Add identifiers, groups and supplementary metadata to a fitted model
#' @export
cca_set_metadata <- function(object, ids = NULL, groups = NULL, supplementary = NULL) {
  if (!inherits(object, "cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  n <- object$n
  if (!is.null(ids) && length(ids) != n) .cca_stop("ids must have one value per fitted observation.")
  if (!is.null(groups) && length(groups) != n) .cca_stop("groups must have one value per fitted observation.")
  if (!is.null(supplementary) && NROW(supplementary) != n) .cca_stop("supplementary must have one row per fitted observation.")
  object$metadata <- list(ids = ids, groups = groups, supplementary = supplementary)
  object
}

#' @export
print.cca_fit <- function(x, ...) {
  cat("CCAflow fit\n")
  cat("  Method:", x$method, "\n")
  cat("  Engine:", x$engine, "\n")
  cat("  Paired observations:", x$n, "\n")
  cat("  Variables: X =", x$p, ", Y =", x$q, "\n")
  cat("  Canonical correlations:", paste(format(round(x$correlations, 4), nsmall = 4), collapse = ", "), "\n")
  if (!identical(x$inference_scope, "classical")) {
    cat("  Inference scope: resampling/method-specific inference is recommended.\n")
  }
  invisible(x)
}

#' @export
summary.cca_fit <- function(object, ...) {
  out <- list(
    method = object$method,
    engine = object$engine,
    dimensions = c(n = object$n, p = object$p, q = object$q),
    correlations = object$correlations,
    variance_extracted_x = object$variance_extracted_x,
    variance_extracted_y = object$variance_extracted_y,
    redundancy_x = object$redundancy_x,
    redundancy_y = object$redundancy_y,
    inference_scope = object$inference_scope,
    tuning = object$tuning
  )
  class(out) <- "summary.cca_fit"
  out
}

#' @export
print.summary.cca_fit <- function(x, ...) {
  cat("CCAflow summary\n")
  cat("  Method:", x$method, "| Engine:", x$engine, "\n")
  print(round(cbind(
    correlation = x$correlations,
    variance_X = x$variance_extracted_x,
    variance_Y = x$variance_extracted_y,
    redundancy_X = x$redundancy_x,
    redundancy_Y = x$redundancy_y
  ), 4))
  invisible(x)
}

#' Predict canonical scores for new paired or single-block data
#' @export
predict.cca_fit <- function(object, newdata_x = NULL, newdata_y = NULL, ...) {
  .cca_check_paired_newdata(newdata_x, newdata_y)
  ans <- list()
  if (!is.null(newdata_x)) {
    X <- .cca_as_matrix(newdata_x, "newdata_x")
    if (ncol(X) != object$p)
      .cca_stop(sprintf("newdata_x has %d variables; the fitted model expects %d.", ncol(X), object$p))
    Xs <- .cca_center_scale(X, center_values = object$center_x, scale_values = object$scale_x)$data
    ans$x <- Xs %*% object$weights_x
  }
  if (!is.null(newdata_y)) {
    Y <- .cca_as_matrix(newdata_y, "newdata_y")
    if (ncol(Y) != object$q)
      .cca_stop(sprintf("newdata_y has %d variables; the fitted model expects %d.", ncol(Y), object$q))
    Ys <- .cca_center_scale(Y, center_values = object$center_y, scale_values = object$scale_y)$data
    ans$y <- Ys %*% object$weights_y
  }
  if (!length(ans)) .cca_stop("Provide newdata_x and/or newdata_y.")
  ans
}
