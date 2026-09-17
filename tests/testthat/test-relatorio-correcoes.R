# Regression tests for the defects reported by the audited tutorial study -----
# Each block corresponds to one item of RELATORIO-AO-AUTOR.md.

test_that("item 1: unknown method names are rejected, never silently empty", {
  d <- cca_agri_data(n = 60, seed = 1)
  expect_error(cca_validate(d$X, d$Y, method = "Classic", folds = 3, repeats = 1),
               "Unknown method")
  expect_error(cca_validate(d$X, d$Y, method = "ridgee", folds = 3, repeats = 1),
               "Unknown method")
  expect_error(cca_tune(d$X, d$Y, method = "zz",
                        grid = expand.grid(lambda_x = 0.1), folds = 3, repeats = 1),
               "Unknown method")
  expect_error(cca_nested_cv(d$X, d$Y, method = "zz",
                             grid = expand.grid(lambda_x = 0.1),
                             outer_folds = 3, inner_folds = 3),
               "Unknown method")
  expect_error(cca_validate_time(d$X, d$Y, method = "zz", initial = 40, assess = 5),
               "Unknown method")
  expect_error(cca_validate_grouped(d$X, d$Y, groups = rep(1:20, each = 3),
                                    method = "zz"),
               "Unknown method")
  expect_error(cca_compare(d$X, d$Y, methods = c("classic", "zz"), validate = FALSE),
               "Unknown method")
})

test_that("item 1: a method that cannot fit warns instead of returning an empty result", {
  d <- cca_agri_data(n = 60, seed = 1)
  ## a deep fit without a model factory fails in every fold: the result is empty
  ## and the caller must be told instead of receiving a silent empty object
  expect_warning(cca_validate(d$X, d$Y, method = "deep", params = list(model = NULL),
                              folds = 3, repeats = 1),
                 "No fold produced a usable fit")
})

test_that("item 2: resistant fits accept a seed, are reproducible and record it", {
  o <- cca_simulate(n = 80, scenario = "outliers", seed = 4)
  a <- cca_robust(o$X, o$Y, ncomp = 1, engine = "robustbase", seed = 7)
  b <- cca_robust(o$X, o$Y, ncomp = 1, engine = "robustbase", seed = 7)
  expect_equal(a$correlations, b$correlations)
  expect_identical(a$tuning$seed, 7)
  m1 <- cca_mcd(o$X, o$Y, ncomp = 1, engine = "robustbase", seed = 11)
  m2 <- cca_mcd(o$X, o$Y, ncomp = 1, engine = "robustbase", seed = 11)
  expect_equal(m1$correlations, m2$correlations)
  expect_identical(m1$refit_params$seed, 11)
  ## the default (deterministic winsor engine) keeps working without a seed
  expect_true(all(is.finite(cca_robust(o$X, o$Y, ncomp = 1)$correlations)))
})

test_that("item 3: a saturated kernel correlation warns about memorisation", {
  ## reporter's minimal example: null scenario plus a wide RBF kernel
  nu <- cca_simulate(n = 60, scenario = "null", seed = 12)
  expect_warning(cca_kernel(nu$X, nu$Y, kernel = "rbf", gamma_x = 10, gamma_y = 10, ncomp = 1),
                 "memorisation")
  ## control: a moderate gamma on real association stays below the ceiling and silent
  ag <- cca_agri_data(seed = 1)
  k <- suppressWarnings(cca_kernel(ag$X, ag$Y, kernel = "rbf", gamma_x = 0.1, gamma_y = 0.1, ncomp = 1))
  expect_lt(max(k$correlations), 1 - 1e-4)
})

test_that("items 4 and 5: paired prediction requires aligned blocks and says what it expects", {
  d <- cca_agri_data(seed = 1)
  fit <- cca_classic(d$X, d$Y, ncomp = 2)
  expect_error(predict(fit, newdata_x = d$X[1:6, ], newdata_y = d$Y[1:4, ]),
               "same number of rows")
  expect_error(predict(fit, newdata_x = d$X[1:6, 1:2]),
               "expects 6")
  expect_error(predict(fit, newdata_y = d$Y[1:6, 1:2]),
               "expects 5")
  pr <- predict(fit, newdata_x = d$X[1:6, ], newdata_y = d$Y[1:6, ])
  expect_equal(nrow(pr$x), nrow(pr$y))
})

test_that("items 7 and 8: bootstrap figures and reports are actually produced", {
  d <- cca_agri_data(seed = 1)
  fit <- cca_classic(d$X, d$Y, ncomp = 2)
  bt <- cca_bootstrap(fit, B = 19, seed = 1)
  skip_if_not_installed("ggplot2")
  for (q in c("loading", "cross_loading", "weight")) {
    p <- plot(bt, quantity = q, block = "x", component = 1, top = 5)
    expect_s3_class(p, "ggplot")
    f <- tempfile(fileext = ".pdf")
    expect_silent(ggplot2::ggsave(f, p, width = 4, height = 3))
    expect_gt(file.size(f), 0)
  }
  p_scree <- cca_plot(fit, "scree", bootstrap = bt)
  expect_s3_class(p_scree, "ggplot")
  f <- tempfile(fileext = ".pdf")
  expect_silent(ggplot2::ggsave(f, p_scree, width = 4, height = 3))
  md <- tempfile(fileext = ".md")
  expect_silent(cca_report(fit, file = md, bootstrap = bt, include_figures = TRUE))
  expect_true(file.exists(md))
  ## include_figures = TRUE must leave figures on disk and reference them
  figdir <- file.path(dirname(md), paste0(tools::file_path_sans_ext(basename(md)), "_figures"))
  expect_gte(length(list.files(figdir, pattern = "[.]png$")), 3)
  expect_true(any(grepl("_figures/", readLines(md, warn = FALSE))))
  ## a companion object that the requested type cannot use must say so
  expect_warning(cca_plot(fit, "scores", bootstrap = bt), "only used by type = 'scree'")
})

test_that("item 6: inference objects are reachable through cca_plot()", {
  d <- cca_agri_data(seed = 1)
  fit <- cca_classic(d$X, d$Y, ncomp = 2)
  skip_if_not_installed("ggplot2")
  bt <- cca_bootstrap(fit, B = 19, seed = 1)
  pv <- suppressWarnings(cca_validate(d$X, d$Y, folds = 3, repeats = 1, seed = 1))
  expect_s3_class(cca_plot(bt), "ggplot")
  expect_s3_class(cca_plot(pv), "ggplot")
  ## the same figure is reachable in both routes, without recursion
  expect_s3_class(cca_plot(bt, quantity = "loading", block = "x", component = 1), "ggplot")
  ## communality closes the last figure type that had no route
  expect_s3_class(cca_plot(fit, "communality"), "ggplot")
  ## an inference key on a fitted object points at the object to build
  expect_error(cca_plot(fit, "bootstrap"), "inference figure")
  ## an unknown key still fails
  expect_error(cca_plot(fit, "naoexiste"))
})

test_that("item 9: multiCCA blocks must keep the documented per-subject form", {
  skip_if_not_installed("multiCCA")
  set.seed(4); n <- 20; Tn <- 6
  errado <- list(A = matrix(rnorm(n * Tn), n * Tn, 1), B = matrix(rnorm(n * Tn), n * Tn, 1))
  expect_error(cca_kernel_gcca(errado, method = "kernel", ncomp = 1), "list with one matrix")
  certo <- list(A = lapply(seq_len(n), function(i) matrix(rnorm(Tn * 2), Tn, 2)),
                B = lapply(seq_len(n), function(i) matrix(rnorm(Tn * 2), Tn, 2)))
  fit <- cca_kernel_gcca(certo, method = "kernel", ncomp = 1)
  expect_equal(fit$fit$n, n)
  ## blocks with different numbers of subjects stay refused
  expect_error(cca_kernel_gcca(list(A = certo$A, B = certo$A[1:10]), method = "kernel", ncomp = 1),
               "same number of subjects")
})

test_that("interface fixes: grand tour, keyframes, repeated columns, sparse message", {
  d <- cca_agri_data(seed = 1)
  fit <- cca_classic(d$X, d$Y, ncomp = 2)
  skip_if_not_installed("tourr")
  gt <- cca_grand_tour(fit)
  expect_output(print(gt), "grand tour")
  expect_error(plot(gt), "no fixed frame set")
  tr <- cca_tour(fit, "x", frames = 8)
  expect_warning(cca_tour_keyframes(tr, n = 3, strategy = "objective"),
                 "objective")
  lon <- data.frame(id = rep(1:10, each = 3), time = rep(1:3, 10),
                    x1 = rnorm(30), x2 = rnorm(30), y1 = rnorm(30))
  expect_error(cca_repeated(lon, id = "id", time = "time", x_vars = "nao_existe",
                            y_vars = "y1", ncomp = 1),
               "not found")
  hd <- cca_simulate(n = 60, scenario = "highdim", seed = 2)
  expect_error(cca_sparse(hd$X, hd$Y, ncomp = 1, penalty_x = 5, penalty_y = 5),
               "penalt")
})

test_that("engine selection is reported when it changes the penalty scale", {
  hd <- cca_simulate(n = 60, scenario = "highdim", seed = 2)
  skip_if_not_installed("PMA")
  expect_warning(cca_sparse(hd$X, hd$Y, engine = "auto", ncomp = 1), "engine")
  expect_silent(suppressWarnings(cca_sparse(hd$X, hd$Y, engine = "internal", ncomp = 1)))
})
