test_that("rolling-origin validation uses estimable test blocks and keeps order", {
  d <- cca_agri_data(n = 80, seed = 201)
  v <- cca_validate_time(d$X, d$Y, method = "ridge",
                         params = list(lambda_x = .1, lambda_y = .1),
                         initial = 45, assess = 5)
  expect_s3_class(v, "cca_validation")
  expect_true(nrow(v$results) > 0)
  for (sp in v$splits) {
    expect_true(max(sp$train) < min(sp$test))
    expect_gte(length(sp$test), 3)
  }
})

test_that("grouped validation never leaks groups across train/test", {
  d <- cca_agri_data(n = 90, seed = 202)
  groups <- rep(seq_len(30), each = 3)
  v <- cca_validate_grouped(d$X, d$Y, groups, method = "classic",
                            folds = 5, repeats = 2, seed = 3)
  expect_true(length(v$splits) > 0)
  for (sp in v$splits)
    expect_length(intersect(sp$train_groups, sp$test_groups), 0)
})

test_that("deep validation can construct a fresh model in every fold", {
  counter <- new.env(parent = emptyenv()); counter$n <- 0L
  factory <- function() {
    counter$n <- counter$n + 1L
    list(id = counter$n)
  }
  fit_fun <- function(model, X, Y, ...) {
    model$mx <- colMeans(X); model$my <- colMeans(Y)
    list(model = model,
         scores_x = matrix(rowMeans(sweep(X, 2, model$mx)), ncol = 1),
         scores_y = matrix(rowMeans(sweep(Y, 2, model$my)), ncol = 1))
  }
  transform_fun <- function(model, X = NULL, Y = NULL, ...) {
    out <- list()
    if (!is.null(X)) out$x <- matrix(rowMeans(sweep(X, 2, model$mx)), ncol = 1)
    if (!is.null(Y)) out$y <- matrix(rowMeans(sweep(Y, 2, model$my)), ncol = 1)
    out
  }
  d <- cca_agri_data(n = 60, seed = 203)
  v <- cca_validate(d$X, d$Y, method = "deep",
                    params = list(model_factory = factory, fit_fun = fit_fun,
                                  transform_fun = transform_fun, ncomp = 1),
                    folds = 5, repeats = 1, seed = 4)
  expect_true(nrow(v$results) > 0)
  expect_gte(counter$n, 5L)
})
