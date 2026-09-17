# Extended latent, mixed, functional, temporal and multiblock CCA ------------

.cca_loglik_gaussian <- function(S, C, n) {
  ee <- eigen((C + t(C)) / 2, symmetric = TRUE)
  vals <- pmax(ee$values, 1e-12)
  invC <- ee$vectors %*% diag(1 / vals, length(vals)) %*% t(ee$vectors)
  -0.5 * n * (nrow(C) * log(2 * pi) + sum(log(vals)) + sum(diag(S %*% invC)))
}

#' Probabilistic CCA fitted by EM
#' @export
cca_probabilistic <- function(X, Y, ncomp = 2, center = TRUE, scale = TRUE,
                              na_action = c("fail", "complete"), max_iter = 500,
                              tol = 1e-6, psi_floor = 1e-6, seed = 1) {
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  Z <- cbind(prep$Xs, prep$Ys)
  n <- nrow(Z); p <- ncol(prep$Xs); q <- ncol(prep$Ys); d <- p + q
  k <- min(as.integer(ncomp), p, q)
  S <- crossprod(Z) / n
  ee <- eigen((S + t(S)) / 2, symmetric = TRUE)
  W <- ee$vectors[, seq_len(k), drop = FALSE] %*%
    diag(sqrt(pmax(ee$values[seq_len(k)] - mean(diag(S)) * 0.25, psi_floor)), k)
  Psi_x <- diag(pmax(diag(S)[seq_len(p)] * 0.5, psi_floor), p)
  Psi_y <- diag(pmax(diag(S)[p + seq_len(q)] * 0.5, psi_floor), q)
  ll <- numeric(max_iter)
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    Psi <- matrix(0, d, d)
    Psi[seq_len(p), seq_len(p)] <- Psi_x
    Psi[p + seq_len(q), p + seq_len(q)] <- Psi_y
    Pinv <- .cca_safe_solve(Psi)
    M <- diag(k) + t(W) %*% Pinv %*% W
    Minv <- .cca_safe_solve(M)
    Ez <- t(Minv %*% t(W) %*% Pinv %*% t(Z))
    Ezz <- n * Minv + crossprod(Ez)
    Wnew <- crossprod(Z, Ez) %*% .cca_safe_solve(Ezz)
    Cxz <- crossprod(Z, Ez) / n
    R <- S - Wnew %*% t(Cxz)
    R <- (R + t(R)) / 2
    Psi_x_new <- R[seq_len(p), seq_len(p), drop = FALSE]
    Psi_y_new <- R[p + seq_len(q), p + seq_len(q), drop = FALSE]
    diag(Psi_x_new) <- pmax(diag(Psi_x_new), psi_floor)
    diag(Psi_y_new) <- pmax(diag(Psi_y_new), psi_floor)
    # Stabilize each residual block if numerical noise makes it indefinite.
    ex <- eigen((Psi_x_new + t(Psi_x_new)) / 2, symmetric = TRUE)
    ey <- eigen((Psi_y_new + t(Psi_y_new)) / 2, symmetric = TRUE)
    Psi_x_new <- ex$vectors %*% diag(pmax(ex$values, psi_floor), p) %*% t(ex$vectors)
    Psi_y_new <- ey$vectors %*% diag(pmax(ey$values, psi_floor), q) %*% t(ey$vectors)
    W <- Wnew; Psi_x <- Psi_x_new; Psi_y <- Psi_y_new
    C <- W %*% t(W)
    C[seq_len(p), seq_len(p)] <- C[seq_len(p), seq_len(p)] + Psi_x
    C[p + seq_len(q), p + seq_len(q)] <- C[p + seq_len(q), p + seq_len(q)] + Psi_y
    ll[it] <- .cca_loglik_gaussian(S, C, n)
    if (it > 1 && abs(ll[it] - ll[it - 1]) <= tol * (1 + abs(ll[it - 1]))) {
      converged <- TRUE
      ll <- ll[seq_len(it)]
      break
    }
  }
  Wx <- W[seq_len(p), , drop = FALSE]
  Wy <- W[p + seq_len(q), , drop = FALSE]
  Sxx <- Wx %*% t(Wx) + Psi_x
  Syy <- Wy %*% t(Wy) + Psi_y
  Sxy <- Wx %*% t(Wy)
  sol <- .cca_cov_solution(Sxx, Syy, Sxy, ncomp = k)
  fit <- .cca_make_fit(
    "probabilistic", prep, sol$xcoef, sol$ycoef, engine = "internal-pCCA-EM",
    tuning = list(ncomp = k, center = center, scale = scale,
                  na_action = match.arg(na_action), max_iter = max_iter,
                  tol = tol, psi_floor = psi_floor, seed = seed),
    call = match.call(),
    extras = list(latent_loading_x = Wx, latent_loading_y = Wy,
                  psi_x = Psi_x, psi_y = Psi_y,
                  log_likelihood = ll, converged = converged,
                  model_cov_x = Sxx, model_cov_y = Syy, model_cov_xy = Sxy),
    inference_scope = "model-based"
  )
  fit$refit_params <- fit$tuning
  fit
}

#' Bayesian probabilistic CCA through Stan
#' @export
cca_bayes <- function(X, Y, ncomp = 2, center = TRUE, scale = TRUE,
                      na_action = c("fail", "complete"), chains = 4, iter = 2000,
                      warmup = floor(iter / 2), seed = 1, refresh = 0,
                      max_posterior_draws = 1000, ...) {
  .cca_require("rstan", "cca_bayes")
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  k <- min(as.integer(ncomp), ncol(prep$Xs), ncol(prep$Ys))
  stan_file <- system.file("stan", "pcca.stan", package = "CCAflow")
  if (!nzchar(stan_file)) .cca_stop("The bundled Stan model could not be located.")
  dat <- list(n = nrow(prep$Xs), p = ncol(prep$Xs), q = ncol(prep$Ys), k = k,
              X = prep$Xs, Y = prep$Ys)
  ## Compiling a Stan model depends on the local toolchain and on a StanHeaders
  ## version matching rstan; when that is broken rstan fails with an opaque I/O
  ## error, so the environment is named explicitly here.
  sf <- tryCatch(
    rstan::stan(file = stan_file, data = dat, chains = chains, iter = iter,
                warmup = warmup, seed = seed, refresh = refresh, ...),
    error = function(e) .cca_stop(sprintf(paste0(
      "The Bayesian engine could not compile or sample its Stan model: %s ",
      "This is an environment problem rather than a model problem: rstan needs ",
      "a working C++ toolchain and a StanHeaders version matching rstan. Check ",
      "the installation with rstan::stan_model(model_code = ",
      "\"parameters { real x; }\")."), conditionMessage(e))))
  ex <- rstan::extract(sf, pars = c("Wx", "Wy", "sigma_x", "sigma_y"), permuted = TRUE)
  nd <- dim(ex$Wx)[1]
  use <- seq_len(min(nd, as.integer(max_posterior_draws)))
  rho_draws <- matrix(NA_real_, length(use), k)
  Sxx_sum <- matrix(0, ncol(prep$Xs), ncol(prep$Xs))
  Syy_sum <- matrix(0, ncol(prep$Ys), ncol(prep$Ys))
  Sxy_sum <- matrix(0, ncol(prep$Xs), ncol(prep$Ys))
  for (ii in seq_along(use)) {
    s <- use[ii]
    Wx <- ex$Wx[s, , , drop = FALSE][1, , ]
    Wy <- ex$Wy[s, , , drop = FALSE][1, , ]
    if (k == 1) { Wx <- matrix(Wx, ncol = 1); Wy <- matrix(Wy, ncol = 1) }
    sx <- as.numeric(ex$sigma_x[s, ]); sy <- as.numeric(ex$sigma_y[s, ])
    Sxx <- Wx %*% t(Wx) + diag(sx^2, length(sx))
    Syy <- Wy %*% t(Wy) + diag(sy^2, length(sy))
    Sxy <- Wx %*% t(Wy)
    sol <- .cca_cov_solution(Sxx, Syy, Sxy, ncomp = k)
    rho_draws[ii, seq_along(sol$cor)] <- sol$cor
    Sxx_sum <- Sxx_sum + Sxx; Syy_sum <- Syy_sum + Syy; Sxy_sum <- Sxy_sum + Sxy
  }
  Sxx <- Sxx_sum / length(use); Syy <- Syy_sum / length(use); Sxy <- Sxy_sum / length(use)
  sol <- .cca_cov_solution(Sxx, Syy, Sxy, ncomp = k)
  ci <- t(apply(rho_draws, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  colnames(ci) <- c("lower", "median", "upper")
  fit <- .cca_make_fit(
    "bayes", prep, sol$xcoef, sol$ycoef, engine = "rstan-pCCA",
    tuning = list(ncomp = k, center = center, scale = scale, chains = chains,
                  iter = iter, warmup = warmup, seed = seed), call = match.call(),
    extras = list(stanfit = sf, posterior_correlations = rho_draws,
                  posterior_correlation_ci = ci,
                  posterior_mean_cov_x = Sxx, posterior_mean_cov_y = Syy,
                  posterior_mean_cov_xy = Sxy), inference_scope = "bayesian"
  )
  class(fit) <- c("cca_bayes", "cca_fit")
  fit
}

#' Deep CCA adapter for a user-configured Python/R deep model
#' @export
cca_deep <- function(X, Y, model = NULL, model_factory = NULL,
                     fit_fun = NULL, transform_fun = NULL,
                     ncomp = 2, center = TRUE, scale = TRUE,
                     na_action = c("fail", "complete"), ...) {
  prep <- .cca_preprocess_pair(X, Y, center, scale, match.arg(na_action))
  if (is.null(model) && !is.null(model_factory)) model <- model_factory()
  if (is.null(model) && is.null(fit_fun)) {
    .cca_stop("cca_deep() requires model, model_factory, or fit_fun. Use model_factory for resampling so that each fold receives a fresh deep model.")
  }
  if (!is.null(fit_fun)) {
    fml <- names(formals(fit_fun))
    if ("model" %in% fml) {
      res <- fit_fun(model = model, X = prep$Xs, Y = prep$Ys, ncomp = ncomp, ...)
    } else {
      res <- fit_fun(prep$Xs, prep$Ys, ncomp = ncomp, ...)
    }
    if (!all(c("scores_x", "scores_y") %in% names(res))) .cca_stop("fit_fun must return scores_x and scores_y.")
    fitted_model <- res$model %||% model
    U <- as.matrix(res$scores_x); V <- as.matrix(res$scores_y)
  } else {
    .cca_require("reticulate", "cca_deep(model=...)")
    model$fit(list(prep$Xs, prep$Ys))
    tr <- model$transform(list(prep$Xs, prep$Ys))
    U <- as.matrix(tr[[1]]); V <- as.matrix(tr[[2]])
    fitted_model <- model
  }
  m <- min(ncomp, ncol(U), ncol(V)); U <- U[, seq_len(m), drop = FALSE]; V <- V[, seq_len(m), drop = FALSE]
  rho <- abs(.cca_cor_cols(U, V))
  out <- list(call = match.call(), method = "deep", engine = "user-configured-deep-model",
              n = nrow(prep$X), p = ncol(prep$X), q = ncol(prep$Y), X = prep$X, Y = prep$Y,
              Xs = prep$Xs, Ys = prep$Ys, center_x = prep$center_x, center_y = prep$center_y,
              scale_x = prep$scale_x, scale_y = prep$scale_y, scores_x = U, scores_y = V,
              correlations = rho, weights_x = NULL, weights_y = NULL,
              raw_weights_x = NULL, raw_weights_y = NULL,
              loadings_x = stats::cor(prep$Xs, U), loadings_y = stats::cor(prep$Ys, V),
              cross_loadings_x = stats::cor(prep$Xs, V), cross_loadings_y = stats::cor(prep$Ys, U),
              variance_extracted_x = colMeans(stats::cor(prep$Xs, U)^2),
              variance_extracted_y = colMeans(stats::cor(prep$Ys, V)^2),
              redundancy_x = colMeans(stats::cor(prep$Xs, U)^2) * rho^2,
              redundancy_y = colMeans(stats::cor(prep$Ys, V)^2) * rho^2,
              communality_x = rowSums(stats::cor(prep$Xs, U)^2), communality_y = rowSums(stats::cor(prep$Ys, V)^2),
              model = fitted_model, model_factory = model_factory,
              fit_fun = fit_fun, transform_fun = transform_fun,
              inference_scope = "resampling",
              tuning = list(ncomp = m, center = center, scale = scale),
              refit_params = list(model_factory = model_factory, fit_fun = fit_fun,
                                  transform_fun = transform_fun, ncomp = m,
                                  center = center, scale = scale))
  class(out) <- c("cca_deep", "cca_fit")
  out
}

#' @export
predict.cca_deep <- function(object, newdata_x = NULL, newdata_y = NULL, ...) {
  .cca_check_paired_newdata(newdata_x, newdata_y)
  if (is.null(newdata_x) || is.null(newdata_y)) .cca_stop("Deep CCA prediction currently requires both paired blocks.")
  X <- .cca_center_scale(.cca_as_matrix(newdata_x, "newdata_x"), center_values = object$center_x, scale_values = object$scale_x)$data
  Y <- .cca_center_scale(.cca_as_matrix(newdata_y, "newdata_y"), center_values = object$center_y, scale_values = object$scale_y)$data
  if (!is.null(object$transform_fun)) return(object$transform_fun(object$model, X, Y, ...))
  .cca_require("reticulate", "predict.cca_deep")
  tr <- object$model$transform(list(X, Y))
  list(x = as.matrix(tr[[1]]), y = as.matrix(tr[[2]]))
}


#' CCA-Zoo estimator adapter through reticulate
#'
#' Fits a scikit-learn-style estimator from the optional Python package
#' `cca-zoo`. Linear, nonparametric/kernel, and probabilistic namespaces are
#' supported. Deep models should be supplied through `cca_deep()` because they
#' require explicit encoder and training-loop configuration.
#' @export
cca_zoo <- function(X, Y, model = "CCA",
                    module = c("linear", "nonparametric", "probabilistic"),
                    ncomp = 2, model_args = list(), python = NULL,
                    center = TRUE, scale = TRUE,
                    na_action = c("fail", "complete")) {
  .cca_require("reticulate", "cca_zoo")
  module <- match.arg(module); na_action <- match.arg(na_action)
  if (!is.null(python)) reticulate::use_python(python, required = TRUE)
  has_py <- tryCatch(isTRUE(reticulate::py_available(initialize = FALSE)), error = function(e) FALSE)
  if (!has_py) {
    .cca_stop("No Python interpreter is available in the active reticulate configuration. Configure a Python environment with 'cca-zoo' outside CCAflow, then retry.")
  }
  if (!reticulate::py_module_available("cca_zoo")) {
    .cca_stop("Python package 'cca-zoo' is not available in the active reticulate environment. Install/configure it outside CCAflow, then retry.")
  }
  prep <- .cca_preprocess_pair(X, Y, center, scale, na_action)
  pymod <- reticulate::import(paste0("cca_zoo.", module), convert = TRUE)
  if (!reticulate::py_has_attr(pymod, model)) {
    .cca_stop(sprintf("CCA-Zoo model '%s' was not found in cca_zoo.%s.", model, module))
  }
  cls <- reticulate::py_get_attr(pymod, model)
  args <- c(list(latent_dimensions = as.integer(ncomp)), as.list(model_args))
  estimator <- tryCatch(do.call(cls, args), error = function(e) {
    .cca_stop(paste0("Could not construct the CCA-Zoo estimator: ", conditionMessage(e)))
  })
  views <- reticulate::tuple(prep$Xs, prep$Ys)
  tryCatch(estimator$fit(views), error = function(e) {
    .cca_stop(paste0("CCA-Zoo fit failed: ", conditionMessage(e)))
  })
  tr <- tryCatch(estimator$transform(views), error = function(e) {
    .cca_stop(paste0("CCA-Zoo transform failed: ", conditionMessage(e)))
  })
  U <- as.matrix(tr[[1]]); V <- as.matrix(tr[[2]])
  m <- min(as.integer(ncomp), ncol(U), ncol(V))
  U <- U[, seq_len(m), drop = FALSE]; V <- V[, seq_len(m), drop = FALSE]
  rho_signed <- .cca_cor_cols(U, V)
  sgn <- ifelse(is.finite(rho_signed) & rho_signed < 0, -1, 1)
  V <- sweep(V, 2, sgn, "*")
  rho <- abs(rho_signed)
  cn <- paste0("CV", seq_len(m)); colnames(U) <- colnames(V) <- cn

  Wx <- Wy <- NULL
  if (reticulate::py_has_attr(estimator, "weights_")) {
    ww <- tryCatch(reticulate::py_to_r(reticulate::py_get_attr(estimator, "weights_")), error = function(e) NULL)
    if (is.list(ww) && length(ww) >= 2L) {
      Wx0 <- tryCatch(as.matrix(ww[[1]]), error = function(e) NULL)
      Wy0 <- tryCatch(as.matrix(ww[[2]]), error = function(e) NULL)
      if (!is.null(Wx0) && nrow(Wx0) == ncol(prep$X)) Wx <- Wx0[, seq_len(min(m, ncol(Wx0))), drop = FALSE]
      if (!is.null(Wy0) && nrow(Wy0) == ncol(prep$Y)) {
        Wy <- Wy0[, seq_len(min(m, ncol(Wy0))), drop = FALSE]
        if (ncol(Wy) == m) Wy <- sweep(Wy, 2, sgn, "*")
      }
      if (!is.null(Wx) && ncol(Wx) == m) { rownames(Wx) <- colnames(prep$X); colnames(Wx) <- cn }
      if (!is.null(Wy) && ncol(Wy) == m) { rownames(Wy) <- colnames(prep$Y); colnames(Wy) <- cn }
    }
  }

  Lx <- stats::cor(prep$Xs, U); Ly <- stats::cor(prep$Ys, V)
  Cx <- stats::cor(prep$Xs, V); Cy <- stats::cor(prep$Ys, U)
  vex <- colMeans(Lx^2); vey <- colMeans(Ly^2)
  out <- list(call = match.call(), method = "zoo",
              engine = paste0("CCA-Zoo::", module, "::", model),
              n = nrow(prep$X), p = ncol(prep$X), q = ncol(prep$Y),
              X = prep$X, Y = prep$Y, Xs = prep$Xs, Ys = prep$Ys,
              center_x = prep$center_x, center_y = prep$center_y,
              scale_x = prep$scale_x, scale_y = prep$scale_y,
              weights_x = Wx, weights_y = Wy,
              raw_weights_x = if (!is.null(Wx)) sweep(Wx, 1, prep$scale_x, "/") else NULL,
              raw_weights_y = if (!is.null(Wy)) sweep(Wy, 1, prep$scale_y, "/") else NULL,
              scores_x = U, scores_y = V, correlations = rho,
              loadings_x = Lx, loadings_y = Ly,
              cross_loadings_x = Cx, cross_loadings_y = Cy,
              variance_extracted_x = vex, variance_extracted_y = vey,
              redundancy_x = vex * rho^2, redundancy_y = vey * rho^2,
              communality_x = rowSums(Lx^2), communality_y = rowSums(Ly^2),
              model = estimator, zoo_module = module, zoo_model = model,
              model_args = model_args, inference_scope = "resampling",
              tuning = list(model = model, module = module, ncomp = m,
                            model_args = model_args, center = center, scale = scale),
              refit_params = list(model = model, module = module, ncomp = m,
                                  model_args = model_args, center = center, scale = scale))
  class(out) <- c("cca_zoo", "cca_fit")
  out
}

#' @export
predict.cca_zoo <- function(object, newdata_x = NULL, newdata_y = NULL, ...) {
  if (is.null(newdata_x) || is.null(newdata_y)) .cca_stop("CCA-Zoo prediction requires both paired blocks.")
  .cca_check_paired_newdata(newdata_x, newdata_y)
  X <- .cca_center_scale(.cca_as_matrix(newdata_x, "newdata_x"), center_values = object$center_x, scale_values = object$scale_x)$data
  Y <- .cca_center_scale(.cca_as_matrix(newdata_y, "newdata_y"), center_values = object$center_y, scale_values = object$scale_y)$data
  tr <- object$model$transform(reticulate::tuple(X, Y))
  list(x = as.matrix(tr[[1]]), y = as.matrix(tr[[2]]))
}

#' Semiparametric sparse CCA for mixed data types
#' @export
cca_mixed <- function(X, Y, type_x = "continuous", type_y = "continuous",
                      BICtype = 1, nlamseq = 20, lam.eps = 0.01, ...) {
  .cca_require("mixedCCA", "cca_mixed")
  X <- .cca_as_matrix(X, "X"); Y <- .cca_as_matrix(Y, "Y")
  pair <- .cca_pair(X, Y, finite = TRUE)
  # mixedCCA (>= 1.5) accepts one type per block, not per-variable vectors.
  .cca_collapse_type <- function(ty, name, nc) {
    if (is.factor(ty)) ty <- as.character(ty)
    if (length(ty) == 1L) return(as.character(ty))
    if (length(ty) != nc) .cca_stop(sprintf("%s must be a single type or one type per variable in the block.", name))
    u <- unique(as.character(ty))
    if (length(u) != 1L) .cca_stop(sprintf("The installed mixedCCA engine accepts a single type per block; %s mixes %s.", name, paste(u, collapse = ", ")))
    u
  }
  type1 <- .cca_collapse_type(type_x, "type_x", ncol(pair$X))
  type2 <- .cca_collapse_type(type_y, "type_y", ncol(pair$Y))
  res <- mixedCCA::mixedCCA(pair$X, pair$Y, type1 = type1, type2 = type2,
                            BICtype = BICtype, nlamseq = nlamseq, lam.eps = lam.eps, ...)
  A <- matrix(as.numeric(res$w1), ncol = 1, dimnames = list(colnames(pair$X), "CV1"))
  B <- matrix(as.numeric(res$w2), ncol = 1, dimnames = list(colnames(pair$Y), "CV1"))
  prep <- .cca_preprocess_pair(pair$X, pair$Y, TRUE, TRUE, "fail")
  fit <- .cca_make_fit("mixed", prep, A, B, engine = "mixedCCA::mixedCCA",
                       tuning = list(type_x = type_x, type_y = type_y, BICtype = BICtype,
                                     nlamseq = nlamseq, lam.eps = lam.eps),
                       call = match.call(), extras = list(engine_fit = res, latent_correlation = as.numeric(res$cancor)),
                       inference_scope = "method-specific")
  fit$correlations_observed_score <- fit$correlations
  fit$correlations <- as.numeric(res$cancor)
  class(fit) <- c("cca_mixed", "cca_fit")
  fit
}

.cca_basis_coefficients <- function(block, time, df = 6, degree = 3) {
  if (is.matrix(block) || is.data.frame(block)) block <- list(curve = .cca_as_matrix(block, "functional block"))
  if (!is.list(block) || !length(block)) .cca_stop("Functional blocks must be matrices or lists of matrices.")
  mats <- lapply(block, .cca_as_matrix, name = "functional curve")
  n <- unique(vapply(mats, nrow, integer(1)))
  if (length(n) != 1) .cca_stop("All functional variables in a block must have the same number of subjects.")
  Tn <- unique(vapply(mats, ncol, integer(1)))
  if (length(Tn) != 1) .cca_stop("All functional variables in a block must share the same grid length.")
  if (length(time) != Tn) .cca_stop("time has incompatible length.")
  B <- splines::bs(time, df = df, degree = degree, intercept = TRUE)
  G <- .cca_safe_solve(crossprod(B)) %*% t(B)
  coefs <- lapply(seq_along(mats), function(j) {
    C <- t(G %*% t(mats[[j]]))
    nm <- names(mats)[j] %||% paste0("f", j)
    colnames(C) <- paste0(nm, "_b", seq_len(ncol(C)))
    C
  })
  list(coefs = do.call(cbind, coefs), basis = B, time = time, df = df, degree = degree, names = names(mats))
}

#' Functional CCA using spline basis coefficients
#' @export
cca_functional <- function(X, Y, time_x = seq_len(ncol(as.matrix(if (is.list(X)) X[[1]] else X))),
                           time_y = seq_len(ncol(as.matrix(if (is.list(Y)) Y[[1]] else Y))),
                           basis_df_x = 6, basis_df_y = 6, degree = 3,
                           method = c("ridge", "classic"), ncomp = 2,
                           lambda_x = 0.1, lambda_y = 0.1) {
  method <- match.arg(method)
  fx <- .cca_basis_coefficients(X, time_x, basis_df_x, degree)
  fy <- .cca_basis_coefficients(Y, time_y, basis_df_y, degree)
  if (method == "ridge") fit <- cca_ridge(fx$coefs, fy$coefs, lambda_x, lambda_y, ncomp)
  else fit <- cca_classic(fx$coefs, fy$coefs, ncomp)
  fit$method <- "functional"; fit$engine <- paste0("spline-basis+", method)
  fit$functional_x <- fx; fit$functional_y <- fy; fit$functional_method <- method; fit$base_method <- method
  fit$refit_params <- list(time_x = time_x, time_y = time_y, basis_df_x = basis_df_x,
                           basis_df_y = basis_df_y, degree = degree, method = method,
                           ncomp = ncomp, lambda_x = lambda_x, lambda_y = lambda_y)
  fit$inference_scope <- "resampling"
  class(fit) <- c("cca_functional", "cca_fit")
  fit
}

#' Dynamic lag-aligned CCA
#' @export
cca_dynamic <- function(X, Y, lag = 0L, method = "ridge", ncomp = 2, ...) {
  pair <- .cca_pair(X, Y, finite = TRUE); X <- pair$X; Y <- pair$Y
  lag <- as.integer(lag)
  n <- nrow(X)
  if (abs(lag) >= n - 2) .cca_stop("Absolute lag is too large for the available series length.")
  if (lag > 0) { ix <- seq_len(n - lag); iy <- ix + lag }
  else if (lag < 0) { iy <- seq_len(n + lag); ix <- iy - lag }
  else { ix <- iy <- seq_len(n) }
  fit <- .cca_fit_method(method, X[ix, , drop = FALSE], Y[iy, , drop = FALSE], c(list(ncomp = ncomp), list(...)))
  fit$dynamic <- list(lag = lag, x_index = ix, y_index = iy)
  fit$base_method <- fit$method; fit$method <- "dynamic"
  class(fit) <- c("cca_dynamic", class(fit))
  fit
}

#' Repeated-measures long-to-wide CCA helper
#' @export
cca_repeated <- function(data, id, time, x_vars, y_vars, method = "ridge", ncomp = 2, ...) {
  faltam <- setdiff(c(id, time, x_vars, y_vars), names(data))
  if (length(faltam))
    .cca_stop(sprintf("Column(s) not found in data: %s. Available columns: %s.",
                      paste(faltam, collapse = ", "), paste(names(data), collapse = ", ")))
  if (!is.data.frame(data)) .cca_stop("data must be a data.frame.")
  ids <- data[[id]]; times <- data[[time]]
  uid <- unique(ids); ut <- sort(unique(times))
  make_block <- function(vars) {
    out <- matrix(NA_real_, length(uid), length(vars) * length(ut))
    cn <- character(ncol(out)); k <- 1L
    for (v in vars) for (tt in ut) {
      cn[k] <- paste0(v, "_t", tt)
      for (i in seq_along(uid)) {
        z <- data[ids == uid[i] & times == tt, v]
        if (length(z)) out[i, k] <- mean(as.numeric(z), na.rm = TRUE)
      }
      k <- k + 1L
    }
    colnames(out) <- cn; rownames(out) <- uid; out
  }
  X <- make_block(x_vars); Y <- make_block(y_vars)
  keep <- stats::complete.cases(X, Y)
  fit <- .cca_fit_method(method, X[keep, , drop = FALSE], Y[keep, , drop = FALSE], c(list(ncomp = ncomp), list(...)))
  fit$repeated <- list(id = id, time = time, ids = uid[keep], time_levels = ut,
                       x_vars = x_vars, y_vars = y_vars)
  fit$base_method <- fit$method; fit$method <- "repeated"
  class(fit) <- c("cca_repeated", class(fit))
  fit
}

#' Explicit tensor matricization followed by CCA
#' @export
cca_multiway <- function(X, Y, method = "ridge", ncomp = 2, ...) {
  flatten <- function(A, nm) {
    if (length(dim(A)) < 2L) .cca_stop(sprintf("%s must have observations in the first array dimension.", nm))
    n <- dim(A)[1]; matrix(A, nrow = n, dimnames = list(NULL, paste0(tolower(nm), "_", seq_len(prod(dim(A)[-1])))))
  }
  Xm <- flatten(X, "X"); Ym <- flatten(Y, "Y")
  fit <- .cca_fit_method(method, Xm, Ym, c(list(ncomp = ncomp), list(...)))
  fit$multiway <- list(dim_x = dim(X), dim_y = dim(Y)); fit$base_method <- fit$method; fit$method <- "multiway"
  class(fit) <- c("cca_multiway", class(fit))
  fit
}

#' Generalized CCA through RGCCA
#' @export
cca_gcca <- function(blocks, ncomp = 2, connection = NULL, scheme = "factorial", scale = TRUE, ...) {
  cca_rgcca(blocks, ncomp = ncomp, connection = connection, tau = 1,
            scheme = scheme, scale = scale, method = "gcca", ...)
}

#' Regularized generalized CCA adapter
#' @export
cca_rgcca <- function(blocks, ncomp = 2, connection = NULL, tau = "optimal",
                      sparsity = 1, scheme = "factorial", scale = TRUE,
                      method = "rgcca", ...) {
  .cca_require("RGCCA", "cca_rgcca")
  fit <- RGCCA::rgcca(blocks = blocks, connection = connection, tau = tau,
                      ncomp = ncomp, scheme = scheme, scale = scale,
                      method = method, sparsity = sparsity, ...)
  tau_efetivo <- tryCatch(fit$tau, error = function(e) NULL)
  if (is.null(tau_efetivo)) tau_efetivo <- tau
  ncomp_pedido <- ncomp
  ncomp_efetivo <- tryCatch(ncol(fit$Y[[1]]), error = function(e) NA_integer_)
  if (!is.na(ncomp_efetivo) && !is.null(ncomp_pedido) && ncomp_efetivo < ncomp_pedido)
    .cca_warn(sprintf(paste0("ncomp was reduced from %d to %d, the rank supported by the ",
                             "blocks. The fitted object records the effective value."),
                      as.integer(ncomp_pedido), as.integer(ncomp_efetivo)))
  out <- list(method = method, engine = "RGCCA::rgcca",
              fit = fit, blocks = blocks, components = fit$Y, weights = fit$a,
              adjusted_weights = fit$astar, AVE = fit$AVE, call = match.call(),
              tau = tau_efetivo, ncomp = ncomp_efetivo)
  class(out) <- c("cca_multiblock", "cca_rgcca")
  out
}

#' Sparse generalized CCA adapter
#' @export
cca_sgcca <- function(blocks, ncomp = 2, sparsity = 0.8, connection = NULL,
                      scheme = "factorial", scale = TRUE, ...) {
  cca_rgcca(blocks, ncomp = ncomp, connection = connection, tau = 1,
            sparsity = sparsity, scheme = scheme, scale = scale, method = "sgcca", ...)
}

#' Tensor generalized CCA through the TGCCA engine in RGCCA
#' @export
cca_tgcca <- function(blocks, ncomp = 2, connection = NULL, tau = 1, rank = 1,
                      mode_orth = 1, separable = TRUE, scheme = "factorial",
                      scale = TRUE, ...) {
  .cca_require("RGCCA", "cca_tgcca")
  has_tgcca <- exists("tgcca", envir = asNamespace("RGCCA"), inherits = FALSE) ||
    "tgcca" %in% tryCatch(RGCCA::available_methods(), error = function(e) character())
  if (!isTRUE(has_tgcca)) {
    .cca_stop("The installed RGCCA version does not provide the TGCCA engine (TGCCA was removed from RGCCA 3.x). Install a TGCCA-capable RGCCA version, or use the tensor-flattening workflow through cca_multiway() instead.")
  }
  cca_rgcca(blocks, ncomp = ncomp, connection = connection, tau = tau,
            scheme = scheme, scale = scale, method = "tgcca",
            rank = rank, mode_orth = mode_orth, separable = separable, ...)
}

#' Multiple kernel or functional CCA through multiCCA
#' @export
cca_kernel_gcca <- function(blocks, method = c("kernel", "functional"), ncomp = 2,
                            eps = 0.001, gamma = NULL, nbasis = 5,
                            basis_type = "fourier", argvals = NULL, ...) {
  .cca_require("multiCCA", "cca_kernel_gcca")
  method <- match.arg(method)
  ## multiCCA expects one matrix per subject inside each block when
  ## method = "kernel". A plain matrix would be read as a single subject and
  ## inflate n to rows x columns, silently; method = "functional" keeps its
  ## documented matrix shape.
  ## is read as a single subject and inflates n to rows x columns, silently.
  por_sujeito <- function(b) is.list(b) && length(b) > 0L &&
    all(vapply(b, function(m) is.matrix(m) || is.data.frame(m), logical(1)))
  if (identical(method, "kernel") && !all(vapply(blocks, por_sujeito, logical(1)))) {
    .cca_stop(paste0("For cca_kernel_gcca(), each block must be a list with one matrix per ",
                     "subject (rows = time, columns = variables). A matrix would be read as a ",
                     "single subject and the sample size would be inflated to rows x columns."))
  }
  sujeitos <- if (identical(method, "kernel")) unique(vapply(blocks, length, integer(1))) else unique(vapply(blocks, NROW, integer(1)))
  if (length(sujeitos) > 1L)
    .cca_stop(sprintf(paste0("All blocks must contain the same number of subjects ",
                             "(%s were given)."), paste(sujeitos, collapse = ", ")))
  fit <- multiCCA::mcca_fit(method = method, X = blocks, ncomp = ncomp, eps = eps,
                            gamma = gamma, nbasis = nbasis, basis_type = basis_type,
                            argvals = argvals, ...)
  out <- list(method = paste0(method, "_mcca"), engine = "multiCCA::mcca_fit",
              fit = fit, blocks = blocks, call = match.call())
  class(out) <- c("cca_multiblock", "cca_multiCCA")
  out
}

#' Tune RGCCA hyperparameters by permutation
#' @export
cca_multiblock_tune <- function(blocks, engine = c("RGCCA", "multiCCA"),
                                par_type = c("tau", "sparsity", "ncomp"),
                                par_value = NULL, par_length = 10, n_perms = 50,
                                ncomp = 1, connection = NULL, method = "rgcca",
                                multi_method = c("kernel", "functional"),
                                eps_grid = c(0.0001, 0.001, 0.01), gamma_grid = NULL,
                                nbasis_grid = NULL, ncomp_eval = 2,
                                criterion = c("first", "sumk"), basis_type = "fourier",
                                argvals = NULL, ...) {
  engine <- match.arg(engine)
  if (engine == "RGCCA") {
    .cca_require("RGCCA", "cca_multiblock_tune(engine='RGCCA')")
    par_type <- match.arg(par_type)
    return(RGCCA::rgcca_permutation(blocks = blocks, par_type = par_type, par_value = par_value,
                                    par_length = par_length, n_perms = n_perms, ncomp = ncomp,
                                    connection = connection, method = method, ...))
  }
  .cca_require("multiCCA", "cca_multiblock_tune(engine='multiCCA')")
  multi_method <- match.arg(multi_method); criterion <- match.arg(criterion)
  multiCCA::mcca_grid_search(method = multi_method, X = blocks, eps_grid = eps_grid,
                             gamma_grid = gamma_grid, nbasis_grid = nbasis_grid,
                             ncomp_eval = ncomp_eval, criterion = criterion,
                             basis_type = basis_type, argvals = argvals, ...)
}

#' Supervised RGCCA parameter tuning by cross-validation
#' @export
cca_multiblock_cv <- function(blocks, response, method = "rgcca",
                              par_type = c("tau", "sparsity", "ncomp"),
                              par_value = NULL, par_length = 10,
                              validation = c("kfold", "loo"),
                              prediction_model = "lm", metric = NULL,
                              k = 5, n_run = 1, n_cores = 1,
                              connection = NULL, ...) {
  .cca_require("RGCCA", "cca_multiblock_cv")
  par_type <- match.arg(par_type); validation <- match.arg(validation)
  ## The internal lm.fit warns once per partition about 'connection'; that warning
  ## belongs to the engine, not to the caller's request.
  withCallingHandlers(
  RGCCA::rgcca_cv(blocks = blocks, connection = connection, method = method,
                  response = response, par_type = par_type, par_value = par_value,
                  par_length = par_length, validation = validation,
                  prediction_model = prediction_model, metric = metric, k = k,
                  n_run = n_run, n_cores = n_cores, ...),
  warning = function(w) {
    if (grepl("extra argument", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}

#' RGCCA bootstrap or SGCCA stability adapter
#' @export
cca_multiblock_stability <- function(object, type = c("bootstrap", "stability"),
                                     n_boot = 200, ...) {
  .cca_require("RGCCA", "cca_multiblock_stability")
  type <- match.arg(type)
  fit <- if (inherits(object, "cca_multiblock")) object$fit else object
  if (type == "bootstrap") return(RGCCA::rgcca_bootstrap(fit, n_boot = n_boot, ...))
  estavel <- tryCatch(RGCCA::rgcca_stability(fit, n_boot = n_boot, ...),
                      error = function(e) e)
  if (inherits(estavel, "error"))
    .cca_stop(paste0("RGCCA stability selection requires a sparse fit. Refit with ",
                     "cca_sgcca() (or cca_multiblock_tune(par_type = 'sparsity')) and pass that ",
                     "object; the engine reported: ", conditionMessage(estavel)))
  estavel
}

#' @export
print.cca_multiblock <- function(x, ...) {
  cat("CCAflow multiblock fit\n  Method:", x$method, "\n  Engine:", x$engine, "\n")
  if (!is.null(x$fit)) print(x$fit)
  invisible(x)
}
