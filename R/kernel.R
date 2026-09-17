# Kernel and nonlinear CCA --------------------------------------------------

.cca_kernel_matrix <- function(X, Z = NULL, kernel = c("rbf", "linear", "polynomial"),
                               gamma = NULL, degree = 2, coef0 = 1) {
  kernel <- match.arg(kernel)
  X <- .cca_as_matrix(X, "X")
  if (is.null(Z)) Z <- X else Z <- .cca_as_matrix(Z, "Z")
  if (ncol(X) != ncol(Z)) .cca_stop("Kernel inputs must have the same number of columns.")
  if (kernel == "linear") return(X %*% t(Z))
  if (kernel == "polynomial") return((X %*% t(Z) + coef0)^degree)
  if (is.null(gamma)) gamma <- 1 / max(ncol(X), 1)
  xx <- rowSums(X^2); zz <- rowSums(Z^2)
  D2 <- outer(xx, zz, "+") - 2 * X %*% t(Z)
  D2[D2 < 0] <- 0
  exp(-gamma * D2)
}

.cca_kernel_center_train <- function(K) {
  cm <- colMeans(K); gm <- mean(K)
  Kc <- sweep(K, 1, rowMeans(K), "-")
  Kc <- sweep(Kc, 2, cm, "-") + gm
  list(Kc = (Kc + t(Kc)) / 2, col_mean = cm, grand_mean = gm)
}

.cca_kernel_center_new <- function(Knew, train_col_mean, train_grand_mean) {
  sweep(sweep(Knew, 1, rowMeans(Knew), "-"), 2, train_col_mean, "-") + train_grand_mean
}

.cca_kpca_fit <- function(Kc, tol = 1e-10) {
  ee <- eigen((Kc + t(Kc)) / 2, symmetric = TRUE)
  keep <- ee$values > max(ee$values, 1) * tol
  if (!any(keep)) .cca_stop("Centered kernel has no positive eigenvalues above tolerance.")
  U <- ee$vectors[, keep, drop = FALSE]
  d <- ee$values[keep]
  scores <- U %*% diag(sqrt(d), length(d))
  list(vectors = U, values = d, scores = scores)
}

.cca_kpca_project <- function(Knew_c, kpca) {
  Knew_c %*% kpca$vectors %*% diag(1 / sqrt(kpca$values), length(kpca$values))
}

#' Kernel canonical correlation analysis
#' @export
cca_kernel <- function(X, Y, kernel = c("rbf", "linear", "polynomial"),
                       gamma_x = NULL, gamma_y = NULL, degree = 2, coef0 = 1,
                       lambda_x = 0.1, lambda_y = 0.1, ncomp = 2,
                       center = TRUE, scale = TRUE,
                       na_action = c("fail", "complete"), tol = 1e-10) {
  kernel <- match.arg(kernel)
  prep0 <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  Kx <- .cca_kernel_matrix(prep0$Xs, kernel = kernel, gamma = gamma_x, degree = degree, coef0 = coef0)
  Ky <- .cca_kernel_matrix(prep0$Ys, kernel = kernel, gamma = gamma_y, degree = degree, coef0 = coef0)
  cx <- .cca_kernel_center_train(Kx); cy <- .cca_kernel_center_train(Ky)
  kpx <- .cca_kpca_fit(cx$Kc, tol); kpy <- .cca_kpca_fit(cy$Kc, tol)
  # CCA in finite kernel-PCA coordinates; prediction uses the exact training embedding.
  zfit <- cca_ridge(kpx$scores, kpy$scores, lambda_x, lambda_y,
                    ncomp = min(ncomp, ncol(kpx$scores), ncol(kpy$scores)),
                    center = FALSE, scale = FALSE, na_action = "fail", tol = tol)
  out <- zfit
  out$method <- "kernel"
  out$engine <- "internal-kernel-pca+ridge-cca"
  out$feature_weights_x <- out$weights_x; out$feature_weights_y <- out$weights_y
  out$weight_space <- "kernel_pca"
  out$original_X <- prep0$X; out$original_Y <- prep0$Y
  out$X <- prep0$X; out$Y <- prep0$Y; out$Xs <- prep0$Xs; out$Ys <- prep0$Ys
  out$n <- nrow(prep0$X); out$p <- ncol(prep0$X); out$q <- ncol(prep0$Y)
  out$loadings_x <- stats::cor(prep0$Xs, out$scores_x, use = "pairwise.complete.obs")
  out$loadings_y <- stats::cor(prep0$Ys, out$scores_y, use = "pairwise.complete.obs")
  out$cross_loadings_x <- stats::cor(prep0$Xs, out$scores_y, use = "pairwise.complete.obs")
  out$cross_loadings_y <- stats::cor(prep0$Ys, out$scores_x, use = "pairwise.complete.obs")
  out$variance_extracted_x <- colMeans(out$loadings_x^2, na.rm = TRUE)
  out$variance_extracted_y <- colMeans(out$loadings_y^2, na.rm = TRUE)
  out$redundancy_x <- out$variance_extracted_x * out$correlations^2
  out$redundancy_y <- out$variance_extracted_y * out$correlations^2
  out$communality_x <- rowSums(out$loadings_x^2, na.rm = TRUE)
  out$communality_y <- rowSums(out$loadings_y^2, na.rm = TRUE)
  out$input_center_x <- prep0$center_x; out$input_center_y <- prep0$center_y
  out$input_scale_x <- prep0$scale_x; out$input_scale_y <- prep0$scale_y
  out$kernel <- kernel
  out$kernel_params <- list(gamma_x = gamma_x, gamma_y = gamma_y, degree = degree, coef0 = coef0)
  out$kernel_train_x <- prep0$Xs; out$kernel_train_y <- prep0$Ys
  out$kernel_center_x <- cx[c("col_mean", "grand_mean")]
  out$kernel_center_y <- cy[c("col_mean", "grand_mean")]
  out$kpca_x <- kpx; out$kpca_y <- kpy
  out$tuning <- list(kernel = kernel, gamma_x = gamma_x, gamma_y = gamma_y,
                     degree = degree, coef0 = coef0, lambda_x = lambda_x,
                     lambda_y = lambda_y, ncomp = ncol(out$weights_x),
                     center = center, scale = scale, na_action = match.arg(na_action), tol = tol)
  out$refit_params <- out$tuning
  out$call <- match.call()
  out$inference_scope <- "resampling"
  ## A radial kernel with a large gamma sends the Gram matrix towards the identity
  ## and the fit interpolates: the training correlation approaches 1 even when the
  ## blocks are independent. The tolerance is loose on purpose (1e-4): the reported
  ## symptom was rho = 0.99999, not exactly 1, and a genuine association never sits
  ## that close to the ceiling at these sample sizes.
  if (any(out$correlations > 1 - 1e-4)) {
    k <- which.max(out$correlations)
    .cca_warn(sprintf(paste0("Canonical correlation %d is %.6f, within 1e-4 of the ceiling. ",
                             "With a wide radial kernel this indicates memorisation rather than ",
                             "association: inspect out-of-sample performance with cca_validate() ",
                             "before interpreting the fit, and consider a smaller gamma or a ",
                             "larger lambda."), k, out$correlations[k]))
  }
  class(out) <- c("cca_kernel", "cca_fit")
  out
}

#' @export
predict.cca_kernel <- function(object, newdata_x = NULL, newdata_y = NULL, ...) {
  .cca_check_paired_newdata(newdata_x, newdata_y)
  ans <- list()
  if (!is.null(newdata_x)) {
    X <- .cca_as_matrix(newdata_x, "newdata_x")
    ## Checked before scaling: otherwise the failure happens inside the
    ## preprocessing and the message does not name the expected width.
    if (ncol(X) != object$p)
      .cca_stop(sprintf(paste0("newdata_x has %d variables; the fitted kernel model ",
                             "expects %d."), ncol(X), object$p))
    Xs <- .cca_center_scale(X, center_values = object$input_center_x, scale_values = object$input_scale_x)$data
    K <- .cca_kernel_matrix(Xs, object$kernel_train_x, object$kernel,
                            object$kernel_params$gamma_x, object$kernel_params$degree,
                            object$kernel_params$coef0)
    Kc <- .cca_kernel_center_new(K, object$kernel_center_x$col_mean, object$kernel_center_x$grand_mean)
    Z <- .cca_kpca_project(Kc, object$kpca_x)
    ans$x <- Z %*% (object$feature_weights_x %||% object$weights_x)
  }
  if (!is.null(newdata_y)) {
    Y <- .cca_as_matrix(newdata_y, "newdata_y")
    if (ncol(Y) != object$q)
      .cca_stop(sprintf(paste0("newdata_y has %d variables; the fitted kernel model ",
                             "expects %d."), ncol(Y), object$q))
    Ys <- .cca_center_scale(Y, center_values = object$input_center_y, scale_values = object$input_scale_y)$data
    K <- .cca_kernel_matrix(Ys, object$kernel_train_y, object$kernel,
                            object$kernel_params$gamma_y, object$kernel_params$degree,
                            object$kernel_params$coef0)
    Kc <- .cca_kernel_center_new(K, object$kernel_center_y$col_mean, object$kernel_center_y$grand_mean)
    Z <- .cca_kpca_project(Kc, object$kpca_y)
    ans$y <- Z %*% (object$feature_weights_y %||% object$weights_y)
  }
  if (!length(ans)) .cca_stop("Provide newdata_x and/or newdata_y.")
  ans
}

#' Resistant preprocessing followed by kernel CCA
#' @export
cca_robust_kernel <- function(X, Y, kernel = c("rbf", "linear", "polynomial"),
                              gamma_x = NULL, gamma_y = NULL, degree = 2, coef0 = 1,
                              lambda_x = 0.1, lambda_y = 0.1, ncomp = 2,
                              robust_engine = c("winsor", "robustbase", "MASS", "auto"),
                              winsor = 4.5, tol = 1e-10, seed = NULL) {
  rp <- .cca_with_seed(seed, .cca_robust_preprocess(X, Y, match.arg(robust_engine), winsor))
  fit <- cca_kernel(rp$Xs, rp$Ys, kernel = match.arg(kernel), gamma_x = gamma_x, gamma_y = gamma_y,
                    degree = degree, coef0 = coef0, lambda_x = lambda_x, lambda_y = lambda_y,
                    ncomp = ncomp, center = FALSE, scale = FALSE, tol = tol)
  fit$method <- "robust_kernel"
  fit$engine <- paste0("robust-", rp$engine, "+", fit$engine)
  fit$robust_input_center_x <- rp$center_x; fit$robust_input_center_y <- rp$center_y
  fit$robust_input_scale_x <- rp$scale_x; fit$robust_input_scale_y <- rp$scale_y
  fit$robust_distances <- rp$robust_distances; fit$robust_cov <- rp$raw_cov
  fit$original_X <- rp$X; fit$original_Y <- rp$Y; fit$X <- rp$X; fit$Y <- rp$Y; fit$p <- ncol(rp$X); fit$q <- ncol(rp$Y); fit$n <- nrow(rp$X)
  fit$tuning$robust_engine <- rp$engine; fit$tuning$winsor <- winsor; fit$tuning$seed <- seed
  fit$refit_params <- fit$tuning
  class(fit) <- c("cca_robust_kernel", "cca_kernel", "cca_fit")
  fit
}

#' @export
predict.cca_robust_kernel <- function(object, newdata_x = NULL, newdata_y = NULL, ...) {
  .cca_check_paired_newdata(newdata_x, newdata_y)
  ## Counts are checked here, before any scaling, so the message names the
  ## expectation instead of failing inside the preprocessing.
  if (!is.null(newdata_x) && ncol(.cca_as_matrix(newdata_x, "newdata_x")) != object$p)
    .cca_stop(sprintf(paste0("newdata_x has %d variables; the fitted kernel model ",
                           "expects %d."), ncol(.cca_as_matrix(newdata_x, "newdata_x")), object$p))
  if (!is.null(newdata_y) && ncol(.cca_as_matrix(newdata_y, "newdata_y")) != object$q)
    .cca_stop(sprintf(paste0("newdata_y has %d variables; the fitted kernel model ",
                           "expects %d."), ncol(.cca_as_matrix(newdata_y, "newdata_y")), object$q))
  # Convert raw new observations to the resistant standardized scale that was
  # supplied to cca_kernel(), whose own input centering/scaling is identity.
  nx <- if (is.null(newdata_x)) NULL else .cca_center_scale(.cca_as_matrix(newdata_x, "newdata_x"),
      center_values = object$robust_input_center_x, scale_values = object$robust_input_scale_x)$data
  ny <- if (is.null(newdata_y)) NULL else .cca_center_scale(.cca_as_matrix(newdata_y, "newdata_y"),
      center_values = object$robust_input_center_y, scale_values = object$robust_input_scale_y)$data
  predict.cca_kernel(object, newdata_x = nx, newdata_y = ny, ...)
}

#' Nonlinear CCA dispatcher
#' @export
cca_nonlinear <- function(X, Y, method = c("kernel", "rank", "copula", "deep"), ...) {
  method <- match.arg(method)
  switch(method,
         kernel = cca_kernel(X, Y, ...),
         rank = cca_rank(X, Y, ...),
         copula = cca_copula(X, Y, ...),
         deep = cca_deep(X, Y, ...))
}
