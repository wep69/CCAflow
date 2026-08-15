# Audit, inference, resampling and stability --------------------------------

#' Audit paired blocks before fitting CCA
#' @export
cca_audit <- function(X, Y, near_cor = 0.95) {
  pair <- .cca_pair(X, Y, finite = FALSE); X <- pair$X; Y <- pair$Y
  audit_block <- function(Z, label) {
    finite_rows <- apply(is.finite(Z), 1, all)
    s <- apply(Z, 2, stats::sd, na.rm = TRUE)
    zero <- names(s)[!is.finite(s) | s <= sqrt(.Machine$double.eps)]
    complete <- Z[finite_rows & stats::complete.cases(Z), , drop = FALSE]
    cond <- NA_real_; maxcor <- NA_real_; pairs <- character()
    if (nrow(complete) >= 3 && ncol(complete) >= 2 && !length(zero)) {
      S <- stats::cov(scale(complete))
      ev <- eigen((S + t(S))/2, symmetric = TRUE, only.values = TRUE)$values
      pos <- ev[ev > sqrt(.Machine$double.eps)]
      if (length(pos)) cond <- max(pos) / min(pos)
      C <- abs(stats::cor(complete))
      diag(C) <- 0
      maxcor <- max(C, na.rm = TRUE)
      ij <- which(C >= near_cor & upper.tri(C), arr.ind = TRUE)
      if (nrow(ij)) pairs <- apply(ij, 1, function(z) paste(colnames(Z)[z[1]], colnames(Z)[z[2]], sep = " ~ "))
    }
    data.frame(block = label, n = nrow(Z), variables = ncol(Z),
               missing = sum(is.na(Z)), nonfinite = sum(!is.finite(Z) & !is.na(Z)),
               zero_variance = length(zero), condition_number = cond,
               maximum_abs_correlation = maxcor, stringsAsFactors = FALSE) -> summary
    list(summary = summary, zero_variance_variables = zero, high_correlation_pairs = pairs)
  }
  ax <- audit_block(X, "X"); ay <- audit_block(Y, "Y")
  out <- list(summary = rbind(ax$summary, ay$summary),
              details = list(X = ax, Y = ay),
              paired_complete = sum(stats::complete.cases(X, Y) & apply(is.finite(cbind(X, Y)), 1, all)),
              p_total = ncol(X) + ncol(Y), n = nrow(X),
              ratio_variables_per_observation = (ncol(X) + ncol(Y)) / nrow(X))
  class(out) <- "cca_audit"
  out
}

#' @export
print.cca_audit <- function(x, ...) {
  cat("CCAflow data audit\n")
  print(x$summary, row.names = FALSE)
  cat("Paired complete observations:", x$paired_complete, "\n")
  cat("Total-variable / observation ratio:", round(x$ratio_variables_per_observation, 3), "\n")
  invisible(x)
}

#' Sequential classical CCA inference
#' @export
cca_infer <- function(object, rao = FALSE) {
  if (!inherits(object, "cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  if (!identical(object$inference_scope, "classical") || object$method != "classic") {
    .cca_warn("Wilks/Bartlett asymptotic inference is returned only for classical CCA. Use cca_permutation() and/or cca_bootstrap() for this fitted method.")
    return(structure(list(method = object$method, scope = object$inference_scope,
                          message = "Classical asymptotic inference not applied."), class = "cca_inference"))
  }
  r <- pmin(pmax(object$correlations, 0), 1 - 1e-12)
  m <- length(r); p <- object$p; q <- object$q; n <- object$n
  cfac <- n - 1 - (p + q + 1) / 2
  if (!is.finite(cfac) || cfac <= 0) {
    .cca_warn("The Bartlett chi-square approximation is not numerically meaningful for this sample-size/dimension combination; chi-square statistics and p-values are returned as NA.")
  }
  tab <- do.call(rbind, lapply(seq_len(m), function(k) {
    lambda <- prod(1 - r[k:m]^2)
    chisq <- if (is.finite(cfac) && cfac > 0) -cfac * log(lambda) else NA_real_
    df <- (p - k + 1) * (q - k + 1)
    pval <- if (is.finite(chisq)) stats::pchisq(chisq, df = df, lower.tail = FALSE) else NA_real_
    data.frame(component = k, correlation = r[k], wilks_lambda = lambda,
               bartlett_chisq = chisq, df = df, p_value = pval)
  }))
  rao_res <- NULL
  if (isTRUE(rao)) {
    if (requireNamespace("CCP", quietly = TRUE)) {
      rao_res <- tryCatch(CCP::p.asym(r, n, p, q, tstat = "Wilks"),
                          error = function(e) structure(list(error = conditionMessage(e)), class = "cca_engine_error"))
    } else .cca_warn("Package 'CCP' is not installed; Rao/CCP output was skipped.")
  }
  out <- list(table = tab, rao = rao_res, scope = "classical asymptotic")
  class(out) <- "cca_inference"
  out
}

#' @export
print.cca_inference <- function(x, ...) {
  if (!is.null(x$table)) print(x$table, row.names = FALSE) else cat(x$message, "\n")
  invisible(x)
}

.cca_perm_fit <- function(fit, X, Y) {
  tryCatch(.cca_refit_like(fit, X, Y), error = function(e) NULL)
}

#' Permutation inference for canonical association
#'
#' For advanced estimators this is a resampling assessment of fitted canonical
#' association. For classical CCA, `stepwise = TRUE` sequentially deflates
#' earlier observed canonical score directions before testing later dimensions.
#' Exchangeability remains a scientific assumption of the design.
#' @export
cca_permutation <- function(object, B = 999, seed = NULL, stepwise = TRUE,
                            max_stat = TRUE) {
  if (!inherits(object, "cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  B <- as.integer(B); if (!is.finite(B) || B < 1L) .cca_stop("B must be a positive integer.")
  if (B < 99L) .cca_warn("Small B: permutation p-values will have coarse resolution.")
  m <- length(object$correlations)
  null <- matrix(NA_real_, B, m)
  .cca_with_seed(seed, {
    for (b in seq_len(B)) {
      perm <- sample.int(object$n)
      ff <- .cca_perm_fit(object, object$X, object$Y[perm, , drop = FALSE])
      if (!is.null(ff)) null[b, seq_len(min(m, length(ff$correlations)))] <- ff$correlations[seq_len(min(m, length(ff$correlations)))]
    }
  })
  p_raw <- vapply(seq_len(m), function(k) (1 + sum(null[, k] >= object$correlations[k], na.rm = TRUE)) /
                    (1 + sum(is.finite(null[, k]))), numeric(1))
  if (isTRUE(max_stat)) {
    p_adj <- vapply(seq_len(m), function(k) {
      mx <- apply(null[, k:m, drop = FALSE], 1, function(z) {
        z <- z[is.finite(z)]
        if (length(z)) max(z) else NA_real_
      })
      (1 + sum(mx >= object$correlations[k], na.rm = TRUE)) / (1 + sum(is.finite(mx)))
    }, numeric(1))
  } else p_adj <- p_raw
  step <- NULL
  if (isTRUE(stepwise) && object$method == "classic" && m > 1) {
    step <- matrix(NA_real_, B, m)
    Xd <- object$Xs; Yd <- object$Ys
    for (k in seq_len(m)) {
      if (k > 1) {
        u <- object$scores_x[, k - 1]; v <- object$scores_y[, k - 1]
        Xd <- Xd - (u %*% crossprod(u, Xd)) / sum(u^2)
        Yd <- Yd - (v %*% crossprod(v, Yd)) / sum(v^2)
      }
      .cca_with_seed(if (is.null(seed)) NULL else seed + k, {
        for (b in seq_len(B)) {
          perm <- sample.int(nrow(Yd))
          f <- tryCatch(cca_classic(Xd, Yd[perm, , drop = FALSE], ncomp = 1,
                                    center = FALSE, scale = FALSE), error = function(e) NULL)
          if (!is.null(f)) step[b, k] <- f$correlations[1]
        }
      })
    }
  }
  p_step <- if (is.null(step)) rep(NA_real_, m) else vapply(seq_len(m), function(k)
    (1 + sum(step[, k] >= object$correlations[k], na.rm = TRUE)) / (1 + sum(is.finite(step[, k]))), numeric(1))
  out <- list(observed = object$correlations, null = null, p_value = p_raw,
              p_maxT = p_adj, stepwise_null = step, p_stepwise = p_step,
              B = B, seed = seed, fit = object,
              note = "Permutation validity requires exchangeability. With nuisance residualization or complex dependence, use a design-specific permutation scheme.")
  class(out) <- "cca_permutation"
  out
}

#' Circular-shift permutation for time-ordered paired blocks
#' @export
cca_permutation_time <- function(object, B = 999, min_shift = 1, seed = NULL) {
  if (!inherits(object, "cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  B <- as.integer(B); if (!is.finite(B) || B < 1L) .cca_stop("B must be a positive integer.")
  if (B < 99L) .cca_warn("Small B: permutation p-values will have coarse resolution.")
  n <- object$n; m <- length(object$correlations)
  if (n < 6) .cca_stop("Too few observations for circular-shift inference.")
  shifts <- seq.int(min_shift, n - min_shift)
  if (!length(shifts)) .cca_stop("No admissible circular shifts.")
  null <- matrix(NA_real_, B, m)
  used <- integer(B)
  .cca_with_seed(seed, {
    for (b in seq_len(B)) {
      s <- sample(shifts, 1); used[b] <- s
      idx <- c((s + 1):n, seq_len(s))
      ff <- .cca_perm_fit(object, object$X, object$Y[idx, , drop = FALSE])
      if (!is.null(ff)) null[b, seq_len(min(m, length(ff$correlations)))] <- ff$correlations[seq_len(min(m, length(ff$correlations)))]
    }
  })
  pv <- vapply(seq_len(m), function(k) (1 + sum(null[, k] >= object$correlations[k], na.rm = TRUE)) /
                 (1 + sum(is.finite(null[, k]))), numeric(1))
  out <- list(observed = object$correlations, null = null, p_value = pv, shifts = used,
              B = B, seed = seed, fit = object,
              note = "Circular shifts preserve within-series ordering but are appropriate only when this null matches the scientific design.")
  class(out) <- c("cca_permutation_time", "cca_permutation")
  out
}

.cca_match_components <- function(reference, candidate) {
  R <- rbind(reference$loadings_x, reference$loadings_y)
  C <- rbind(candidate$loadings_x, candidate$loadings_y)
  m <- min(ncol(R), ncol(C))
  R <- R[, seq_len(m), drop = FALSE]; C <- C[, seq_len(m), drop = FALSE]
  sim <- abs(stats::cor(R, C, use = "pairwise.complete.obs"))
  ord <- integer(m); available <- seq_len(m)
  for (k in seq_len(m)) {
    j <- available[which.max(sim[k, available])]
    ord[k] <- j; available <- setdiff(available, j)
  }
  signs <- vapply(seq_len(m), function(k) {
    z <- sum(R[, k] * C[, ord[k]], na.rm = TRUE)
    if (z < 0) -1 else 1
  }, numeric(1))
  list(order = ord, sign = signs)
}

.cca_boot_long <- function(M, replicate, quantity, block) {
  if (is.null(M)) return(data.frame())
  d <- .cca_long_matrix(M, block, quantity); d$replicate <- replicate; d
}

#' Bootstrap CCA uncertainty with component alignment
#' @export
cca_bootstrap <- function(object, B = 999, conf = 0.95, seed = NULL) {
  if (!inherits(object, "cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  B <- as.integer(B); if (!is.finite(B) || B < 2L) .cca_stop("B must be at least 2 for bootstrap intervals.")
  m <- length(object$correlations)
  cors <- matrix(NA_real_, B, m)
  longs <- vector("list", B * 6); pos <- 1L
  redx <- matrix(NA_real_, B, m); redy <- matrix(NA_real_, B, m)
  .cca_with_seed(seed, {
    for (b in seq_len(B)) {
      idx <- sample.int(object$n, object$n, replace = TRUE)
      ff <- tryCatch(.cca_refit_like(object, object$X[idx, , drop = FALSE], object$Y[idx, , drop = FALSE]), error = function(e) NULL)
      if (is.null(ff)) next
      mm <- min(m, length(ff$correlations))
      al <- .cca_match_components(object, ff)
      ord <- al$order[seq_len(mm)]; sg <- al$sign[seq_len(mm)]
      cors[b, seq_len(mm)] <- ff$correlations[ord]
      for (k in seq_len(mm)) {
        ff$loadings_x[, ord[k]] <- ff$loadings_x[, ord[k]] * sg[k]
        ff$loadings_y[, ord[k]] <- ff$loadings_y[, ord[k]] * sg[k]
        ff$cross_loadings_x[, ord[k]] <- ff$cross_loadings_x[, ord[k]] * sg[k]
        ff$cross_loadings_y[, ord[k]] <- ff$cross_loadings_y[, ord[k]] * sg[k]
        if (!is.null(ff$weights_x) && nrow(ff$weights_x) == object$p) ff$weights_x[, ord[k]] <- ff$weights_x[, ord[k]] * sg[k]
        if (!is.null(ff$weights_y) && nrow(ff$weights_y) == object$q) ff$weights_y[, ord[k]] <- ff$weights_y[, ord[k]] * sg[k]
      }
      redx[b, seq_len(mm)] <- ff$redundancy_x[ord]
      redy[b, seq_len(mm)] <- ff$redundancy_y[ord]
      ref_names <- paste0("CV", seq_len(mm))
      Lx <- ff$loadings_x[, ord, drop = FALSE]; colnames(Lx) <- ref_names
      Ly <- ff$loadings_y[, ord, drop = FALSE]; colnames(Ly) <- ref_names
      Cx <- ff$cross_loadings_x[, ord, drop = FALSE]; colnames(Cx) <- ref_names
      Cy <- ff$cross_loadings_y[, ord, drop = FALSE]; colnames(Cy) <- ref_names
      longs[[pos]] <- .cca_boot_long(Lx, b, "loading", "x"); pos <- pos + 1L
      longs[[pos]] <- .cca_boot_long(Ly, b, "loading", "y"); pos <- pos + 1L
      longs[[pos]] <- .cca_boot_long(Cx, b, "cross_loading", "x"); pos <- pos + 1L
      longs[[pos]] <- .cca_boot_long(Cy, b, "cross_loading", "y"); pos <- pos + 1L
      if (!is.null(ff$weights_x) && nrow(ff$weights_x) == object$p) { Wx <- ff$weights_x[, ord, drop = FALSE]; colnames(Wx) <- ref_names; longs[[pos]] <- .cca_boot_long(Wx, b, "weight", "x"); pos <- pos + 1L }
      if (!is.null(ff$weights_y) && nrow(ff$weights_y) == object$q) { Wy <- ff$weights_y[, ord, drop = FALSE]; colnames(Wy) <- ref_names; longs[[pos]] <- .cca_boot_long(Wy, b, "weight", "y"); pos <- pos + 1L }
    }
  })
  alpha <- (1 - conf) / 2
  ci <- t(apply(cors, 2, stats::quantile, probs = c(alpha, 0.5, 1 - alpha), na.rm = TRUE))
  colnames(ci) <- c("lower", "median", "upper")
  correlation_ci <- data.frame(component = seq_len(m), estimate = object$correlations,
                               lower = ci[, 1], median = ci[, 2], upper = ci[, 3])
  used_longs <- if (pos > 1L) Filter(function(z) !is.null(z) && nrow(z), longs[seq_len(pos - 1L)]) else list()
  variable_draws <- if (length(used_longs)) do.call(rbind, used_longs) else data.frame()
  out <- list(fit = object, B = B, conf = conf, correlations = cors,
              correlation_ci = correlation_ci, variable_draws = variable_draws,
              redundancy_x = redx, redundancy_y = redy, seed = seed)
  class(out) <- "cca_bootstrap"
  out
}

#' Cluster bootstrap preserving independent experimental units
#' @export
cca_bootstrap_grouped <- function(object, groups, B = 999, conf = 0.95, seed = NULL) {
  if (!inherits(object, "cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  if (length(groups) != object$n) .cca_stop("groups must have one value per observation.")
  B <- as.integer(B); if (!is.finite(B) || B < 2L) .cca_stop("B must be at least 2 for grouped bootstrap intervals.")
  ug <- unique(groups); m <- length(object$correlations); cors <- matrix(NA_real_, B, m)
  .cca_with_seed(seed, {
    for (b in seq_len(B)) {
      gs <- sample(ug, length(ug), replace = TRUE)
      idx <- unlist(lapply(gs, function(g) which(groups == g)), use.names = FALSE)
      ff <- tryCatch(.cca_refit_like(object, object$X[idx, , drop = FALSE], object$Y[idx, , drop = FALSE]), error = function(e) NULL)
      if (!is.null(ff)) {
        mm <- min(m, length(ff$correlations))
        al <- tryCatch(.cca_match_components(object, ff), error = function(e) NULL)
        ord <- if (is.null(al)) seq_len(mm) else al$order[seq_len(mm)]
        cors[b, seq_len(mm)] <- ff$correlations[ord]
      }
    }
  })
  alpha <- (1-conf)/2
  ci <- t(apply(cors, 2, stats::quantile, probs = c(alpha, 0.5, 1-alpha), na.rm = TRUE))
  out <- list(fit = object, groups = groups, B = B, conf = conf, correlations = cors,
              correlation_ci = data.frame(component = seq_len(m), estimate = object$correlations,
                                           lower = ci[,1], median = ci[,2], upper = ci[,3]), seed = seed)
  class(out) <- c("cca_bootstrap_grouped", "cca_bootstrap")
  out
}

#' Bootstrap variable-selection stability
#' @export
cca_stability <- function(object, B = 500, threshold = 1e-8, seed = NULL) {
  bt <- cca_bootstrap(object, B = B, seed = seed)
  d <- bt$variable_draws
  d <- d[d$quantity == "weight", , drop = FALSE]
  if (!nrow(d)) .cca_stop("Weight-space stability is unavailable for this method. Use loading bootstrap intervals instead.")
  key <- interaction(d$block, d$variable, d$component, drop = TRUE)
  sp <- split(d, key)
  tab <- do.call(rbind, lapply(sp, function(z) data.frame(
    block = z$block[1], variable = z$variable[1], component = z$component[1],
    selection_frequency = mean(abs(z$value) > threshold, na.rm = TRUE),
    median_weight = stats::median(z$value, na.rm = TRUE))))
  rownames(tab) <- NULL
  out <- list(table = tab, bootstrap = bt, threshold = threshold)
  class(out) <- "cca_stability"
  out
}

#' Select the number of canonical components
#' @export
cca_select_components <- function(object, method = c("wilks", "permutation", "cv"),
                                  alpha = 0.05, B = 499, X = NULL, Y = NULL,
                                  folds = 5, repeats = 3, seed = NULL,
                                  cv_conf = 0.95) {
  method <- match.arg(method)
  if (method == "wilks") {
    inf <- cca_infer(object)
    if (is.null(inf$table)) .cca_stop("Wilks selection is available only for classical CCA.")
    keep <- which(inf$table$p_value <= alpha)
    ncomp <- if (length(keep)) max(keep) else 0L
    return(list(method = method, ncomp = ncomp, evidence = inf$table))
  }
  if (method == "permutation") {
    pm <- cca_permutation(object, B = B, seed = seed)
    keep <- which(pm$p_maxT <= alpha)
    ncomp <- if (length(keep)) max(keep) else 0L
    return(list(method = method, ncomp = ncomp, evidence = data.frame(component = seq_along(pm$p_maxT), p = pm$p_maxT)))
  }
  if (is.null(X)) X <- object$X; if (is.null(Y)) Y <- object$Y
  cv <- cca_validate(X, Y, method = object$method, params = object$refit_params %||% list(), folds = folds, repeats = repeats, seed = seed)
  if (!nrow(cv$results)) return(list(method = method, ncomp = 0L, evidence = data.frame()))
  alpha_cv <- (1 - cv_conf) / 2
  sp <- split(cv$results$correlation, cv$results$component)
  sm <- do.call(rbind, lapply(names(sp), function(k) {
    z <- sp[[k]]; z <- z[is.finite(z)]
    if (!length(z)) return(data.frame(component = as.integer(k), mean = NA_real_, sd = NA_real_, lower = NA_real_, upper = NA_real_, positive_fraction = NA_real_))
    qq <- stats::quantile(z, probs = c(alpha_cv, 1 - alpha_cv), na.rm = TRUE, names = FALSE)
    data.frame(component = as.integer(k), mean = mean(z), sd = stats::sd(z),
               lower = qq[1], upper = qq[2], positive_fraction = mean(z > 0))
  }))
  sm <- sm[order(sm$component), , drop = FALSE]
  # A sequential CV generalization rule: retain only the leading prefix whose
  # held-out correlations remain positive across the central resampling interval.
  prefix <- cumprod(as.integer(is.finite(sm$lower) & sm$lower > 0)) == 1
  ncomp <- if (any(prefix)) max(sm$component[prefix]) else 0L
  list(method = method, ncomp = ncomp, evidence = sm, cv_conf = cv_conf,
       note = "CV component selection is a generalization heuristic, not a formal hypothesis test.")
}
