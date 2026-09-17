# Sparse, structured and resistant CCA --------------------------------------

.cca_sparse_internal <- function(X, Y, ncomp = 1L, penalty_x = 0.1, penalty_y = 0.1,
                                 ridge_x = 0, ridge_y = 0, groups_x = NULL, groups_y = NULL,
                                 max_iter = 200L, tol = 1e-7) {
  Xd <- X; Yd <- Y
  p <- ncol(X); q <- ncol(Y)
  ncomp <- min(as.integer(ncomp), p, q)
  A <- matrix(0, p, ncomp); B <- matrix(0, q, ncomp)
  soft_group <- function(z, lambda, groups) {
    if (is.null(groups)) return(.cca_soft(z, lambda))
    if (length(groups) != length(z)) .cca_stop("Group vector length does not match the block dimension.")
    out <- numeric(length(z))
    for (g in unique(groups)) {
      ii <- which(groups == g); zg <- z[ii]; ng <- sqrt(sum(zg^2))
      mult <- if (ng > 0) max(0, 1 - lambda * sqrt(length(ii)) / ng) else 0
      out[ii] <- mult * zg
    }
    out
  }
  for (k in seq_len(ncomp)) {
    sv <- svd(crossprod(Xd, Yd), nu = 1, nv = 1)
    a <- sv$u[, 1]; b <- sv$v[, 1]
    for (it in seq_len(max_iter)) {
      a0 <- a; b0 <- b
      b <- crossprod(Yd, Xd %*% a) / nrow(Xd)
      b <- soft_group(as.vector(b), penalty_y, groups_y) / (1 + ridge_y)
      b <- .cca_norm(b)
      a <- crossprod(Xd, Yd %*% b) / nrow(Xd)
      a <- soft_group(as.vector(a), penalty_x, groups_x) / (1 + ridge_x)
      a <- .cca_norm(a)
      if (sum(a^2) <= 0 || sum(b^2) <= 0) break
      if (max(abs(c(a - a0, b - b0))) < tol) break
    }
    A[, k] <- a; B[, k] <- b
    u <- as.vector(Xd %*% a); v <- as.vector(Yd %*% b)
    if (sum(u^2) > 0) Xd <- Xd - (u %*% crossprod(u, Xd)) / sum(u^2)
    if (sum(v^2) > 0) Yd <- Yd - (v %*% crossprod(v, Yd)) / sum(v^2)
  }
  if (!any(A != 0) || !any(B != 0)) {
    if (!is.null(groups_x) || !is.null(groups_y)) {
      .cca_stop(sprintf(paste0("Sparsity removed all weights (penalty_x = %g, penalty_y = %g). ",
                               "Reduce the penalties or the group thresholds and retry."),
                        penalty_x, penalty_y))
    }
    .cca_stop(sprintf(paste0("Sparsity removed all weights (penalty_x = %g, penalty_y = %g). ",
                             "Reduce the penalties, or regularise instead of thresholding with ",
                             "cca_ridge()/cca_elastic()."), penalty_x, penalty_y))
  }
  list(A = A, B = B)
}

#' Sparse canonical correlation analysis
#' @export
cca_sparse <- function(X, Y, ncomp = 1, penalty_x = 0.1, penalty_y = 0.1,
                       engine = c("internal", "PMA", "auto"),
                       center = TRUE, scale = TRUE,
                       na_action = c("fail", "complete"), max_iter = 200, tol = 1e-7) {
  engine <- match.arg(engine)
  if (engine == "auto") {
    engine <- if (requireNamespace("PMA", quietly = TRUE)) "PMA" else "internal"
    .cca_warn(sprintf(paste0("cca_sparse(engine = 'auto') selected the '%s' engine. The ",
                             "penalty scales differ between engines (internal thresholds are ",
                             "absolute; PMA penalties lie in (0, 1]). Set engine explicitly when ",
                             "comparing penalties across machines."), engine))
  }
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  if (engine == "PMA") {
    .cca_require("PMA", "cca_sparse(engine='PMA')")
    px <- if (missing(penalty_x)) 0.7 else penalty_x
    py <- if (missing(penalty_y)) 0.7 else penalty_y
    if (px <= 0 || px > 1 || py <= 0 || py > 1) .cca_stop("For PMA, penalty_x and penalty_y must lie in (0, 1]. Larger values impose less sparsity.")
    # Store the effective PMA penalties so bootstrap/CV refits reproduce the
    # model actually fitted, including the PMA-specific defaults.
    penalty_x <- px; penalty_y <- py
    pm <- PMA::CCA(prep$Xs, prep$Ys, typex = "standard", typez = "standard",
                   penaltyx = px, penaltyz = py, K = ncomp,
                   niter = max(15L, min(as.integer(max_iter), 100L)),
                   trace = FALSE, standardize = FALSE)
    A <- as.matrix(pm$u); B <- as.matrix(pm$v)
    eng <- pm
  } else {
    if (penalty_x < 0 || penalty_y < 0) .cca_stop("Internal sparsity thresholds must be non-negative.")
    sp <- .cca_sparse_internal(prep$Xs, prep$Ys, ncomp, penalty_x, penalty_y,
                               max_iter = max_iter, tol = tol)
    A <- sp$A; B <- sp$B; eng <- NULL
  }
  fit <- .cca_make_fit(
    "sparse", prep, A, B, engine = if (engine == "PMA") "PMA::CCA" else "internal-alternating-threshold",
    tuning = list(ncomp = ncol(A), penalty_x = penalty_x, penalty_y = penalty_y,
                  engine = engine, center = center, scale = scale,
                  na_action = match.arg(na_action), max_iter = max_iter, tol = tol),
    call = match.call(), extras = list(engine_fit = eng), inference_scope = "resampling"
  )
  fit$refit_params <- fit$tuning
  fit
}

#' Elastic-net-inspired sparse CCA
#' @export
cca_elastic <- function(X, Y, ncomp = 1, l1_x = 0.1, l1_y = 0.1,
                        l2_x = 0.1, l2_y = 0.1, center = TRUE, scale = TRUE,
                        na_action = c("fail", "complete"), max_iter = 200, tol = 1e-7) {
  if (min(l1_x, l1_y, l2_x, l2_y) < 0) .cca_stop("Elastic-net penalties must be non-negative.")
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  sp <- .cca_sparse_internal(prep$Xs, prep$Ys, ncomp, l1_x, l1_y,
                             ridge_x = l2_x, ridge_y = l2_y,
                             max_iter = max_iter, tol = tol)
  fit <- .cca_make_fit(
    "elastic", prep, sp$A, sp$B, engine = "internal-elastic-alternating",
    tuning = list(ncomp = ncomp, l1_x = l1_x, l1_y = l1_y, l2_x = l2_x, l2_y = l2_y,
                  center = center, scale = scale, na_action = match.arg(na_action),
                  max_iter = max_iter, tol = tol), call = match.call(), inference_scope = "resampling"
  )
  fit$refit_params <- fit$tuning
  fit
}

#' Group-structured sparse CCA
#' @export
cca_group_sparse <- function(X, Y, groups_x, groups_y, ncomp = 1,
                             penalty_x = 0.1, penalty_y = 0.1,
                             center = TRUE, scale = TRUE,
                             na_action = c("fail", "complete"), max_iter = 200, tol = 1e-7) {
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  if (length(groups_x) != ncol(prep$Xs) || length(groups_y) != ncol(prep$Ys)) {
    .cca_stop("groups_x and groups_y must identify one group per variable in their respective blocks.")
  }
  sp <- .cca_sparse_internal(prep$Xs, prep$Ys, ncomp, penalty_x, penalty_y,
                             groups_x = groups_x, groups_y = groups_y,
                             max_iter = max_iter, tol = tol)
  fit <- .cca_make_fit(
    "group_sparse", prep, sp$A, sp$B, engine = "internal-group-threshold",
    tuning = list(groups_x = groups_x, groups_y = groups_y, ncomp = ncomp,
                  penalty_x = penalty_x, penalty_y = penalty_y,
                  center = center, scale = scale, na_action = match.arg(na_action),
                  max_iter = max_iter, tol = tol), call = match.call(), inference_scope = "resampling"
  )
  fit$refit_params <- fit$tuning
  fit
}

.cca_robust_preprocess <- function(X, Y, engine = c("winsor", "robustbase", "MASS", "auto"), winsor = 4.5) {
  pair <- .cca_pair(X, Y, finite = TRUE)
  X <- pair$X; Y <- pair$Y; Z <- cbind(X, Y); p <- ncol(X)
  engine <- match.arg(engine)
  if (engine == "auto") {
    engine <- if (requireNamespace("robustbase", quietly = TRUE)) "robustbase" else if (requireNamespace("MASS", quietly = TRUE)) "MASS" else "winsor"
  }
  if (engine == "robustbase") {
    .cca_require("robustbase", "robust CCA")
    rr <- robustbase::covMcd(Z)
    ctr <- rr$center; S <- rr$cov; distances <- rr$mah
  } else if (engine == "MASS") {
    .cca_require("MASS", "robust CCA")
    rr <- MASS::cov.rob(Z, method = "mcd")
    ctr <- rr$center; S <- rr$cov
    distances <- stats::mahalanobis(Z, ctr, S)
  } else {
    ctr <- apply(Z, 2, stats::median)
    scl0 <- apply(Z, 2, stats::mad, constant = 1.4826)
    bad_scl <- !is.finite(scl0) | scl0 <= sqrt(.Machine$double.eps)
    scl0[bad_scl] <- .cca_sd(Z)[bad_scl]
    scl0[!is.finite(scl0) | scl0 <= sqrt(.Machine$double.eps)] <- 1
    Zw <- sweep(Z, 2, ctr, "-"); Zw <- sweep(Zw, 2, scl0, "/")
    Zw <- pmax(pmin(Zw, winsor), -winsor)
    Sstd <- stats::cov(Zw)
    S <- diag(scl0) %*% Sstd %*% diag(scl0)
    distances <- rowSums(Zw^2)
    rr <- list(center = ctr, cov = S, winsorized = Zw)
  }
  S <- (S + t(S)) / 2
  scl <- sqrt(pmax(diag(S), .Machine$double.eps))
  Xs <- sweep(sweep(X, 2, ctr[seq_len(p)], "-"), 2, scl[seq_len(p)], "/")
  jy <- p + seq_len(ncol(Y))
  Ys <- sweep(sweep(Y, 2, ctr[jy], "-"), 2, scl[jy], "/")
  D <- diag(1 / scl)
  Sstd <- D %*% S %*% D
  list(
    X = X, Y = Y, Xs = Xs, Ys = Ys,
    center_x = ctr[seq_len(p)], center_y = ctr[jy],
    scale_x = scl[seq_len(p)], scale_y = scl[jy],
    cov = Sstd, raw_cov = S, robust_center = ctr, robust_scale = scl,
    robust_distances = distances, engine = engine, engine_fit = rr
  )
}

#' Resistant canonical correlation analysis
#' @export
cca_robust <- function(X, Y, ncomp = NULL,
                       engine = c("winsor", "robustbase", "MASS", "auto"),
                       winsor = 4.5, tol = 1e-10, seed = NULL) {
  ## The MCD engines draw random subsets: seed makes the fit reproducible and
  ## recordable. The internal winsor engine is deterministic either way.
  rp <- .cca_with_seed(seed, .cca_robust_preprocess(X, Y, match.arg(engine), winsor))
  p <- ncol(rp$X); q <- ncol(rp$Y)
  Sxx <- rp$cov[seq_len(p), seq_len(p), drop = FALSE]
  Syy <- rp$cov[p + seq_len(q), p + seq_len(q), drop = FALSE]
  Sxy <- rp$cov[seq_len(p), p + seq_len(q), drop = FALSE]
  sol <- .cca_cov_solution(Sxx, Syy, Sxy, ncomp, tol)
  prep <- rp[c("X", "Y", "Xs", "Ys", "center_x", "center_y", "scale_x", "scale_y")]
  fit <- .cca_make_fit(
    "robust", prep, sol$xcoef, sol$ycoef, engine = paste0("robust-", rp$engine),
    tuning = list(ncomp = ncol(sol$xcoef), engine = rp$engine, winsor = winsor, tol = tol, seed = seed),
    call = match.call(),
    extras = list(robust_cov = rp$raw_cov, robust_center = rp$robust_center,
                  robust_scale = rp$robust_scale, robust_distances = rp$robust_distances,
                  engine_fit = rp$engine_fit), inference_scope = "resampling"
  )
  fit$refit_params <- fit$tuning
  fit
}

#' Minimum covariance determinant CCA convenience wrapper
#' @export
cca_mcd <- function(X, Y, ncomp = NULL, engine = c("auto", "robustbase", "MASS"),
                  tol = 1e-10, seed = NULL) {
  engine <- match.arg(engine)
  if (engine == "auto") engine <- if (requireNamespace("robustbase", quietly = TRUE)) "robustbase" else "MASS"
  fit <- cca_robust(X, Y, ncomp = ncomp, engine = engine, tol = tol, seed = seed)
  fit$method <- "mcd"
  class(fit) <- c("cca_mcd", "cca_robust", "cca_fit")
  fit$refit_params <- list(ncomp = ncol(fit$weights_x), engine = engine, tol = tol, seed = seed)
  fit
}

#' Resistant sparse CCA
#' @export
cca_robust_sparse <- function(X, Y, ncomp = 1, penalty_x = 0.1, penalty_y = 0.1,
                              robust_engine = c("winsor", "robustbase", "MASS", "auto"),
                              winsor = 4.5, max_iter = 200, tol = 1e-7, seed = NULL) {
  rp <- .cca_with_seed(seed, .cca_robust_preprocess(X, Y, match.arg(robust_engine), winsor))
  sp <- .cca_sparse_internal(rp$Xs, rp$Ys, ncomp, penalty_x, penalty_y,
                             max_iter = max_iter, tol = tol)
  prep <- rp[c("X", "Y", "Xs", "Ys", "center_x", "center_y", "scale_x", "scale_y")]
  fit <- .cca_make_fit(
    "robust_sparse", prep, sp$A, sp$B, engine = paste0("robust-", rp$engine, "+internal-sparse"),
    tuning = list(ncomp = ncomp, penalty_x = penalty_x, penalty_y = penalty_y,
                  robust_engine = rp$engine, winsor = winsor, max_iter = max_iter, tol = tol, seed = seed),
    call = match.call(), extras = list(robust_cov = rp$raw_cov, robust_center = rp$robust_center,
                                       robust_scale = rp$robust_scale, robust_distances = rp$robust_distances),
    inference_scope = "resampling"
  )
  fit$refit_params <- fit$tuning
  fit
}

#' Rank-transformed CCA
#' @export
cca_rank <- function(X, Y, ncomp = NULL, center = TRUE, scale = TRUE,
                     na_action = c("fail", "complete"), tol = 1e-10) {
  pair <- .cca_pair(X, Y, finite = FALSE)
  if (anyNA(pair$X) || anyNA(pair$Y)) {
    if (match.arg(na_action) == "fail") .cca_stop("Missing values found.")
    keep <- stats::complete.cases(pair$X, pair$Y)
    pair$X <- pair$X[keep, , drop = FALSE]; pair$Y <- pair$Y[keep, , drop = FALSE]
  }
  mx <- .cca_empirical_map_fit(pair$X, gaussian = FALSE)
  my <- .cca_empirical_map_fit(pair$Y, gaussian = FALSE)
  fit <- cca_classic(mx$train, my$train, ncomp, center, scale, "fail", tol)
  fit$method <- "rank"; fit$engine <- "empirical-rank+classical"
  fit$original_X <- pair$X; fit$original_Y <- pair$Y; fit$X <- pair$X; fit$Y <- pair$Y; fit$transform_x <- mx; fit$transform_y <- my
  fit$inference_scope <- "resampling"
  fit$refit_params <- list(ncomp = ncol(fit$weights_x), center = center, scale = scale,
                           na_action = "fail", tol = tol)
  class(fit) <- c("cca_rank", "cca_fit")
  fit
}

#' Gaussian-copula normal-score CCA
#' @export
cca_copula <- function(X, Y, ncomp = NULL, center = TRUE, scale = TRUE,
                       na_action = c("fail", "complete"), tol = 1e-10) {
  pair <- .cca_pair(X, Y, finite = FALSE)
  if (anyNA(pair$X) || anyNA(pair$Y)) {
    if (match.arg(na_action) == "fail") .cca_stop("Missing values found.")
    keep <- stats::complete.cases(pair$X, pair$Y)
    pair$X <- pair$X[keep, , drop = FALSE]; pair$Y <- pair$Y[keep, , drop = FALSE]
  }
  mx <- .cca_empirical_map_fit(pair$X, gaussian = TRUE)
  my <- .cca_empirical_map_fit(pair$Y, gaussian = TRUE)
  fit <- cca_classic(mx$train, my$train, ncomp, center, scale, "fail", tol)
  fit$method <- "copula"; fit$engine <- "empirical-gaussian-copula+classical"
  fit$original_X <- pair$X; fit$original_Y <- pair$Y; fit$X <- pair$X; fit$Y <- pair$Y; fit$transform_x <- mx; fit$transform_y <- my
  fit$inference_scope <- "resampling"
  fit$refit_params <- list(ncomp = ncol(fit$weights_x), center = center, scale = scale,
                           na_action = "fail", tol = tol)
  class(fit) <- c("cca_copula", "cca_fit")
  fit
}

#' @export
predict.cca_rank <- function(object, newdata_x = NULL, newdata_y = NULL, ...) {
  .cca_check_paired_newdata(newdata_x, newdata_y)
  ans <- list()
  if (!is.null(newdata_x)) {
    Xtr <- .cca_empirical_map_apply(newdata_x, object$transform_x)
    Xs <- .cca_center_scale(Xtr, center_values = object$center_x, scale_values = object$scale_x)$data
    ans$x <- Xs %*% object$weights_x
  }
  if (!is.null(newdata_y)) {
    Ytr <- .cca_empirical_map_apply(newdata_y, object$transform_y)
    Ys <- .cca_center_scale(Ytr, center_values = object$center_y, scale_values = object$scale_y)$data
    ans$y <- Ys %*% object$weights_y
  }
  if (!length(ans)) .cca_stop("Provide newdata_x and/or newdata_y.")
  ans
}

#' @export
predict.cca_copula <- function(object, newdata_x = NULL, newdata_y = NULL, ...) predict.cca_rank(object, newdata_x, newdata_y, ...)

#' Direct comparison of classical and resistant CCA
#' @export
cca_robust_compare <- function(X, Y, ncomp = NULL, robust_engine = "winsor", ...) {
  classical <- cca_classic(X, Y, ncomp = ncomp)
  resistant <- cca_robust(X, Y, ncomp = ncomp, engine = robust_engine, ...)
  m <- min(length(classical$correlations), length(resistant$correlations))
  delta_loading_x <- sqrt(colSums((abs(classical$loadings_x[, seq_len(m), drop = FALSE]) -
                                   abs(resistant$loadings_x[, seq_len(m), drop = FALSE]))^2))
  delta_loading_y <- sqrt(colSums((abs(classical$loadings_y[, seq_len(m), drop = FALSE]) -
                                   abs(resistant$loadings_y[, seq_len(m), drop = FALSE]))^2))
  tab <- data.frame(
    component = seq_len(m),
    classical = classical$correlations[seq_len(m)],
    resistant = resistant$correlations[seq_len(m)],
    delta_rho = classical$correlations[seq_len(m)] - resistant$correlations[seq_len(m)],
    loading_shift_x = delta_loading_x,
    loading_shift_y = delta_loading_y
  )
  out <- list(classical = classical, resistant = resistant, comparison = tab)
  class(out) <- "cca_robust_comparison"
  out
}

#' Observation-level CCA diagnostics
#' @export
cca_diagnose <- function(object, cutoff = 0.975) {
  if (!inherits(object, "cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  Z <- cbind(object$Xs, object$Ys)
  ctr <- colMeans(Z)
  S <- stats::cov(Z)
  dist_classic <- tryCatch(stats::mahalanobis(Z, ctr, S), error = function(e) rep(NA_real_, nrow(Z)))
  rdist <- object$robust_distances %||% rep(NA_real_, nrow(Z))
  score_energy <- rowSums(object$scores_x^2) + rowSums(object$scores_y^2)
  df <- ncol(Z)
  threshold <- stats::qchisq(cutoff, df = max(df, 1))
  d <- data.frame(
    observation = seq_len(nrow(Z)),
    classical_distance = dist_classic,
    robust_distance = rdist,
    score_energy = score_energy,
    flagged_classical = is.finite(dist_classic) & dist_classic > threshold,
    flagged_robust = is.finite(rdist) & rdist > threshold
  )
  out <- list(data = d, cutoff = cutoff, threshold = threshold, fit = object)
  class(out) <- "cca_diagnostics"
  out
}
