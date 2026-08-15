test_that("optional RGCCA adapter follows the current object contract", {
  skip_if_not_installed("RGCCA")
  z <- cca_simulate(n = 60, scenario = "three_blocks", seed = 101)
  fit <- cca_rgcca(z$blocks, ncomp = 1, tau = 1)
  expect_s3_class(fit, "cca_multiblock")
  expect_identical(fit$engine, "RGCCA::rgcca")
  expect_length(fit$blocks, 3)
})

test_that("optional mixedCCA adapter returns a common fit object", {
  skip_if_not_installed("mixedCCA")
  set.seed(102)
  n <- 80
  latent <- rnorm(n)
  X <- cbind(x1 = latent + rnorm(n, sd = .5), x2 = rnorm(n))
  Y <- cbind(y1 = latent + rnorm(n, sd = .5), y2 = rnorm(n))
  fit <- cca_mixed(X, Y, type_x = rep("continuous", 2),
                   type_y = rep("continuous", 2), BICtype = 1)
  expect_s3_class(fit, "cca_fit")
  expect_true(all(is.finite(fit$correlations)))
})

test_that("optional multiCCA adapter accepts documented repeated-measures structure", {
  skip_if_not_installed("multiCCA")
  set.seed(103)
  n <- 10; tt <- 8
  blocks <- list(
    lapply(seq_len(n), function(i) matrix(rnorm(tt * 2), tt, 2)),
    lapply(seq_len(n), function(i) matrix(rnorm(tt * 2), tt, 2)),
    lapply(seq_len(n), function(i) matrix(rnorm(tt * 2), tt, 2))
  )
  fit <- cca_kernel_gcca(blocks, method = "kernel", ncomp = 1)
  expect_s3_class(fit, "cca_multiblock")
})

test_that("optional spinifex renderer consumes a canonical basis path", {
  skip_if_not_installed("spinifex")
  skip_if_not_installed("ggplot2")
  d <- cca_agri_data(n = 50, seed = 104)
  fit <- cca_classic(d$X, d$Y, ncomp = 2)
  tr <- cca_tour(fit, block = "x", frames = 5)
  g <- cca_tour_render(tr, engine = "spinifex", animate = "none")
  expect_false(is.null(g))
})

test_that("optional TGCCA wrapper exposes tensor RGCCA without flattening", {
  skip_if_not_installed("RGCCA")
  if (!("tgcca" %in% RGCCA::available_methods())) skip("TGCCA is unavailable in this RGCCA version.")
  set.seed(105)
  blocks <- list(
    A = array(rnorm(24 * 3 * 2), dim = c(24, 3, 2)),
    B = array(rnorm(24 * 2 * 3), dim = c(24, 2, 3))
  )
  fit <- cca_tgcca(blocks, ncomp = 1, rank = 1, mode_orth = 1,
                   separable = TRUE, scale = TRUE)
  expect_s3_class(fit, "cca_multiblock")
  expect_identical(fit$method, "tgcca")
})


test_that("optional CCA-Zoo bridge returns common scores when Python engine is configured", {
  skip_if_not_installed("reticulate")
  if (!reticulate::py_module_available("cca_zoo")) skip("Python package cca-zoo is unavailable.")
  set.seed(106)
  z <- cca_simulate(n = 50, scenario = "linear", seed = 106)
  fit <- cca_zoo(z$X, z$Y, model = "CCA", module = "linear", ncomp = 1)
  expect_s3_class(fit, "cca_fit")
  expect_s3_class(fit, "cca_zoo")
  expect_true(is.finite(fit$correlations[1]))
})
