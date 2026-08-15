# Internal utilities --------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.cca_stop <- function(...) stop(..., call. = FALSE)
.cca_warn <- function(...) warning(..., call. = FALSE)

.cca_require <- function(pkg, feature = pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .cca_stop(sprintf("'%s' requires the optional package '%s'. Install it explicitly before using this feature.", feature, pkg))
  }
  invisible(TRUE)
}

.cca_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

.cca_as_matrix <- function(x, name = "X") {
  if (is.data.frame(x)) {
    bad <- vapply(x, function(z) is.factor(z) || is.character(z) || is.logical(z), logical(1))
    if (any(bad)) {
      .cca_stop(sprintf("%s contains non-numeric columns (%s). Encode categorical variables explicitly or use cca_mixed() where appropriate.",
                        name, paste(names(x)[bad], collapse = ", ")))
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) .cca_stop(sprintf("%s must be numeric.", name))
  storage.mode(x) <- "double"
  if (is.null(colnames(x))) colnames(x) <- paste0(tolower(name), seq_len(ncol(x)))
  x
}

.cca_pair <- function(X, Y, finite = TRUE) {
  X <- .cca_as_matrix(X, "X")
  Y <- .cca_as_matrix(Y, "Y")
  if (nrow(X) != nrow(Y)) .cca_stop("X and Y must have the same number of rows (paired observations).")
  if (nrow(X) < 3L) .cca_stop("At least three paired observations are required.")
  if (ncol(X) < 1L || ncol(Y) < 1L) .cca_stop("Each block must contain at least one variable.")
  if (finite && (!all(is.finite(X)) || !all(is.finite(Y)))) .cca_stop("X and Y must contain only finite values after missing-data handling.")
  list(X = X, Y = Y)
}

.cca_sd <- function(x) {
  s <- apply(x, 2, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s <= sqrt(.Machine$double.eps)] <- 1
  s
}

.cca_center_scale <- function(X, center = TRUE, scale = TRUE, center_values = NULL, scale_values = NULL) {
  X <- .cca_as_matrix(X, "X")
  ctr <- if (!is.null(center_values)) center_values else if (isTRUE(center)) colMeans(X) else rep(0, ncol(X))
  scl <- if (!is.null(scale_values)) scale_values else if (isTRUE(scale)) .cca_sd(X) else rep(1, ncol(X))
  if (length(ctr) != ncol(X) || length(scl) != ncol(X)) .cca_stop("Preprocessing vectors are incompatible with the number of variables.")
  Z <- sweep(X, 2, ctr, "-")
  Z <- sweep(Z, 2, scl, "/")
  colnames(Z) <- colnames(X)
  list(data = Z, center = ctr, scale = scl)
}

.cca_inv_sqrt <- function(S, tol = 1e-10) {
  S <- (S + t(S)) / 2
  ee <- eigen(S, symmetric = TRUE)
  keep <- ee$values > max(ee$values, 1) * tol
  if (!any(keep)) .cca_stop("Covariance metric is numerically singular after tolerance filtering.")
  V <- ee$vectors[, keep, drop = FALSE]
  d <- ee$values[keep]
  V %*% (diag(1 / sqrt(d), nrow = length(d))) %*% t(V)
}

.cca_safe_solve <- function(A, B = NULL, tol = 1e-10) {
  ee <- eigen((A + t(A)) / 2, symmetric = TRUE)
  keep <- ee$values > max(abs(ee$values), 1) * tol
  if (!any(keep)) .cca_stop("Matrix is numerically singular.")
  inv <- ee$vectors[, keep, drop = FALSE] %*% diag(1 / ee$values[keep], sum(keep)) %*% t(ee$vectors[, keep, drop = FALSE])
  if (is.null(B)) inv else inv %*% B
}

.cca_cov_solution <- function(Sxx, Syy, Sxy, ncomp = NULL, tol = 1e-10) {
  p <- nrow(Sxx); q <- nrow(Syy)
  if (is.null(ncomp)) ncomp <- min(p, q)
  ncomp <- min(as.integer(ncomp), p, q)
  Wx <- .cca_inv_sqrt(Sxx, tol)
  Wy <- .cca_inv_sqrt(Syy, tol)
  M <- Wx %*% Sxy %*% Wy
  sv <- svd(M, nu = ncomp, nv = ncomp)
  list(
    cor = pmin(pmax(sv$d[seq_len(ncomp)], 0), 1),
    xcoef = Wx %*% sv$u[, seq_len(ncomp), drop = FALSE],
    ycoef = Wy %*% sv$v[, seq_len(ncomp), drop = FALSE]
  )
}

.cca_cor_cols <- function(A, B) {
  m <- min(ncol(A), ncol(B))
  vapply(seq_len(m), function(k) {
    z <- suppressWarnings(stats::cor(A[, k], B[, k], use = "pairwise.complete.obs"))
    if (is.finite(z)) z else NA_real_
  }, numeric(1))
}

.cca_align_sign <- function(A, B, U, V) {
  r <- .cca_cor_cols(U, V)
  for (k in seq_along(r)) {
    if (is.finite(r[k]) && r[k] < 0) {
      B[, k] <- -B[, k]
      V[, k] <- -V[, k]
      r[k] <- -r[k]
    }
  }
  list(A = A, B = B, U = U, V = V, cor = r)
}

.cca_loadings <- function(Xs, Ys, U, V) {
  lx <- stats::cor(Xs, U, use = "pairwise.complete.obs")
  ly <- stats::cor(Ys, V, use = "pairwise.complete.obs")
  clx <- stats::cor(Xs, V, use = "pairwise.complete.obs")
  cly <- stats::cor(Ys, U, use = "pairwise.complete.obs")
  list(loadings_x = lx, loadings_y = ly, cross_x = clx, cross_y = cly)
}

.cca_interpretation <- function(Xs, Ys, U, V, rho) {
  ld <- .cca_loadings(Xs, Ys, U, V)
  vx <- colMeans(ld$loadings_x^2, na.rm = TRUE)
  vy <- colMeans(ld$loadings_y^2, na.rm = TRUE)
  rx <- vx * rho^2
  ry <- vy * rho^2
  comm_x <- rowSums(ld$loadings_x^2, na.rm = TRUE)
  comm_y <- rowSums(ld$loadings_y^2, na.rm = TRUE)
  list(
    loadings_x = ld$loadings_x,
    loadings_y = ld$loadings_y,
    cross_loadings_x = ld$cross_x,
    cross_loadings_y = ld$cross_y,
    variance_extracted_x = vx,
    variance_extracted_y = vy,
    redundancy_x = rx,
    redundancy_y = ry,
    communality_x = comm_x,
    communality_y = comm_y
  )
}

.cca_make_fit <- function(method, prep, A, B, engine = "internal", tuning = list(), call = NULL,
                          extras = list(), inference_scope = "classical") {
  U <- prep$Xs %*% A
  V <- prep$Ys %*% B
  al <- .cca_align_sign(A, B, U, V)
  A <- al$A; B <- al$B; U <- al$U; V <- al$V; rho <- al$cor
  good <- is.finite(rho)
  A <- A[, good, drop = FALSE]; B <- B[, good, drop = FALSE]
  U <- U[, good, drop = FALSE]; V <- V[, good, drop = FALSE]
  rho <- abs(rho[good])
  if (!length(rho)) .cca_stop("No finite canonical component could be constructed.")
  cn <- paste0("CV", seq_along(rho))
  colnames(A) <- colnames(B) <- colnames(U) <- colnames(V) <- cn
  rownames(A) <- colnames(prep$X); rownames(B) <- colnames(prep$Y)
  interp <- .cca_interpretation(prep$Xs, prep$Ys, U, V, rho)
  raw_A <- sweep(A, 1, prep$scale_x, "/")
  raw_B <- sweep(B, 1, prep$scale_y, "/")
  out <- c(list(
    call = call,
    method = method,
    engine = engine,
    n = nrow(prep$X),
    p = ncol(prep$X),
    q = ncol(prep$Y),
    X = prep$X,
    Y = prep$Y,
    Xs = prep$Xs,
    Ys = prep$Ys,
    center_x = prep$center_x,
    center_y = prep$center_y,
    scale_x = prep$scale_x,
    scale_y = prep$scale_y,
    weights_x = A,
    weights_y = B,
    raw_weights_x = raw_A,
    raw_weights_y = raw_B,
    scores_x = U,
    scores_y = V,
    correlations = rho,
    tuning = tuning,
    inference_scope = inference_scope,
    metadata = list()
  ), interp, extras)
  class(out) <- c(paste0("cca_", method), "cca_fit")
  out
}

.cca_preprocess_pair <- function(X, Y, center = TRUE, scale = TRUE, na_action = c("fail", "complete")) {
  na_action <- match.arg(na_action)
  pair <- .cca_pair(X, Y, finite = FALSE)
  X <- pair$X; Y <- pair$Y
  keep <- rep(TRUE, nrow(X))
  if (anyNA(X) || anyNA(Y) || !all(is.finite(X)) || !all(is.finite(Y))) {
    if (na_action == "fail") .cca_stop("Missing/non-finite values found. Use na_action = 'complete' or preprocess explicitly.")
    keep <- stats::complete.cases(X, Y) & apply(is.finite(X), 1, all) & apply(is.finite(Y), 1, all)
    X <- X[keep, , drop = FALSE]; Y <- Y[keep, , drop = FALSE]
    if (nrow(X) < 3L) .cca_stop("Fewer than three complete paired observations remain.")
  }
  px <- .cca_center_scale(X, center, scale)
  py <- .cca_center_scale(Y, center, scale)
  list(X = X, Y = Y, Xs = px$data, Ys = py$data,
       center_x = px$center, center_y = py$center,
       scale_x = px$scale, scale_y = py$scale, row_keep = keep)
}

.cca_fit_function <- function(method) {
  switch(method,
    classic = cca_classic,
    ridge = cca_ridge,
    shrinkage = cca_shrinkage,
    sparse = cca_sparse,
    elastic = cca_elastic,
    group_sparse = cca_group_sparse,
    robust = cca_robust,
    mcd = cca_mcd,
    rank = cca_rank,
    copula = cca_copula,
    kernel = cca_kernel,
    robust_sparse = cca_robust_sparse,
    robust_kernel = cca_robust_kernel,
    probabilistic = cca_probabilistic,
    bayes = cca_bayes,
    mixed = cca_mixed,
    deep = cca_deep,
    zoo = cca_zoo,
    NULL
  )
}

.cca_filter_call_args <- function(fn, params) {
  if (!length(params)) return(list())
  fml <- names(formals(fn))
  if ("..." %in% fml) return(params)
  params[intersect(names(params), fml)]
}

.cca_fit_method <- function(method, X, Y, params = list()) {
  allowed <- c("classic", "ridge", "shrinkage", "sparse", "elastic", "group_sparse",
               "robust", "mcd", "rank", "copula", "kernel", "robust_sparse",
               "robust_kernel", "probabilistic", "bayes", "mixed", "deep", "zoo")
  method <- match.arg(method, allowed)
  fn <- .cca_fit_function(method)
  do.call(fn, c(list(X = X, Y = Y), params))
}

.cca_refit_like <- function(fit, X, Y, overrides = list()) {
  pars <- fit$refit_params %||% fit$tuning %||% list()
  pars <- utils::modifyList(pars, overrides)
  method <- fit$method
  # Structured wrappers are resampled in the representation stored in the fit.
  # Their method-specific preprocessing (e.g., spline basis or lag alignment)
  # has already produced X and Y in that analysis space.
  if (method %in% c("functional", "dynamic", "repeated", "multiway")) {
    method <- fit$base_method %||% fit$functional_method %||% "ridge"
  }
  fn <- .cca_fit_function(method)
  if (is.null(fn)) .cca_stop(sprintf("No refit engine is registered for method '%s'.", method))
  pars <- .cca_filter_call_args(fn, pars)
  .cca_fit_method(method, X, Y, pars)
}

.cca_long_matrix <- function(M, block, quantity) {
  if (is.null(M)) return(data.frame())
  d <- expand.grid(variable = rownames(M), component = colnames(M), stringsAsFactors = FALSE)
  d$value <- as.vector(M)
  d$block <- block
  d$quantity <- quantity
  d
}

.cca_empirical_map_fit <- function(X, gaussian = FALSE) {
  lapply(seq_len(ncol(X)), function(j) sort(X[, j])) |> setNames(colnames(X)) -> maps
  transform_train <- apply(X, 2, function(z) {
    r <- rank(z, ties.method = "average")
    u <- (r - 0.5) / length(z)
    if (gaussian) stats::qnorm(u) else u
  })
  if (is.vector(transform_train)) transform_train <- matrix(transform_train, ncol = ncol(X))
  colnames(transform_train) <- colnames(X)
  list(maps = maps, train = transform_train, gaussian = gaussian)
}

.cca_empirical_map_apply <- function(X, map) {
  X <- .cca_as_matrix(X, "newdata")
  if (ncol(X) != length(map$maps)) .cca_stop("New data have an incompatible number of variables.")
  out <- matrix(NA_real_, nrow(X), ncol(X), dimnames = list(rownames(X), names(map$maps)))
  for (j in seq_len(ncol(X))) {
    s <- map$maps[[j]]; n <- length(s)
    u <- (findInterval(X[, j], s) + 0.5) / (n + 1)
    u <- pmin(pmax(u, 1 / (2 * n + 2)), 1 - 1 / (2 * n + 2))
    out[, j] <- if (isTRUE(map$gaussian)) stats::qnorm(u) else u
  }
  out
}

.cca_soft <- function(z, lambda) sign(z) * pmax(abs(z) - lambda, 0)
.cca_norm <- function(z) {
  n <- sqrt(sum(z^2))
  if (!is.finite(n) || n <= sqrt(.Machine$double.eps)) return(rep(0, length(z)))
  z / n
}

.cca_top <- function(x, n = 15L) {
  if (length(x) <= n) seq_along(x) else order(abs(x), decreasing = TRUE)[seq_len(n)]
}
