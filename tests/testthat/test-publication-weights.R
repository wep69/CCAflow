test_that("cca_weight_diagnostics identifies suppressors and handles rank-deficient VIF", {
  d <- cca_agri_data(n = 100, seed = 10)
  fit <- cca_classic(d$X, d$Y, ncomp = 3)
  dg <- cca_weight_diagnostics(fit, components = 1:3)
  expect_s3_class(dg, "cca_weight_diagnostics")
  expect_true(dg$has_weights)
  expect_true(all(c("loading", "cross_loading", "weight", "vif", "suppressor", "collinear_inflation") %in% names(dg$table)))
  expect_true(all(dg$vif$X >= 1) && all(dg$vif$Y >= 1))
  expect_output(print(dg), "Original-variable weights available")
  expect_identical(cca_table(fit, "participation")$variable, cca_weight_diagnostics(fit)$table$variable)

  ## High-dimensional / singular correlation matrix must compute finite VIF without error
  hd <- cca_simulate(n = 25, p = 35, q = 30, scenario = "highdim", seed = 3)
  fit_r <- cca_ridge(hd$X, hd$Y, lambda_x = 0.2, lambda_y = 0.2, ncomp = 2)
  dg_r <- cca_weight_diagnostics(fit_r, components = 1:2)
  expect_true(all(is.finite(dg_r$table$vif)))
  expect_true(all(dg_r$table$vif >= 1))
})

test_that("extended coefficients, sparsity and bootstrap weight plots work cleanly", {
  skip_if_not_installed("ggplot2")
  d <- cca_agri_data(n = 90, seed = 4)
  fit <- cca_classic(d$X, d$Y, ncomp = 3)

  p_coef_both <- cca_plot(fit, "coefficients", block = "both", components = 1:2,
                          style = "lollipop", polarity_color = TRUE)
  expect_s3_class(p_coef_both, "ggplot")

  sp <- cca_sparse(d$X, d$Y, ncomp = 3, penalty_x = 0.05, penalty_y = 0.05)
  expect_s3_class(cca_plot(sp, "sparsity", sparsity_mode = "binary"), "ggplot")
  expect_s3_class(cca_plot(sp, "sparsity", sparsity_mode = "signed"), "ggplot")
  expect_s3_class(cca_plot(sp, "sparsity", sparsity_mode = "magnitude"), "ggplot")

  bt <- cca_bootstrap(fit, B = 15, seed = 2)
  p_bt_both <- plot(bt, quantity = "weight", block = "both", component = 1, highlight_zero = TRUE)
  expect_s3_class(p_bt_both, "ggplot")
})

test_that("suppressor and all four bipolar participation plot styles render and export", {
  skip_if_not_installed("ggplot2")
  d <- cca_agri_data(n = 120, seed = 1)
  fit <- cca_classic(d$X, d$Y, ncomp = 4)
  lbls <- stats::setNames(paste("Var", c(colnames(d$X), colnames(d$Y))), c(colnames(d$X), colnames(d$Y)))

  tipos <- c("suppressor", "bipolar_lollipop", "weight_vs_loading", "tornado",
             "synoptic_heatmap", "heliograph", "bipolar_balance")
  for (tp in tipos) {
    p <- cca_plot(fit, type = tp, components = 1:3, variable_labels = lbls)
    expect_s3_class(p, "ggplot")
    f <- tempfile(fileext = ".pdf")
    expect_no_error(cca_export(p, f, preset = "double"))
    expect_gt(file.size(f), 0)
  }

  expect_s3_class(cca_plot(fit, "synoptic_heatmap", grayscale = TRUE, block = "x"), "ggplot")
  expect_s3_class(cca_weight_diagnostics(fit) |> plot(), "ggplot")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_panel <- cca_plot(fit, type = "heliograph_panel", components = c(2, 3), variable_labels = lbls)
    expect_s3_class(p_panel, "ggplot")
    f_pan <- tempfile(fileext = ".pdf")
    expect_no_error(cca_export(p_panel, f_pan, preset = "double"))
    expect_gt(file.size(f_pan), 0)
  }

  md <- tempfile(fileext = ".md")
  expect_no_error(cca_report(fit, file = md, include_figures = TRUE, extra_figures = TRUE))
  figdir <- file.path(dirname(md), paste0(tools::file_path_sans_ext(basename(md)), "_figures"))
  expect_gte(length(list.files(figdir, pattern = "[.]png$")), 6)
})

test_that("kernel fits refuse weight-only figures and polar layouts reject interactive plotly", {
  skip_if_not_installed("ggplot2")
  nm <- cca_simulate(n = 60, scenario = "nonmonotone", seed = 5)
  fk <- cca_kernel(nm$X, nm$Y, kernel = "rbf", ncomp = 2)
  expect_error(cca_plot(fk, "suppressor"), "original-variable canonical weights")
  expect_error(cca_plot(fk, "weight_vs_loading"), "original-variable canonical weights")
  expect_error(cca_plot(fk, "synoptic_heatmap"), "original-variable canonical weights")
  expect_s3_class(cca_plot(fk, "bipolar_lollipop", components = 1:2), "ggplot")
  expect_s3_class(cca_plot(fk, "heliograph", component = 1), "ggplot")
  expect_s3_class(cca_plot(fk, "bipolar_balance", components = 1:2), "ggplot")

  if (requireNamespace("plotly", quietly = TRUE)) {
    d <- cca_agri_data(n = 80, seed = 2)
    fc <- cca_classic(d$X, d$Y, ncomp = 3)
    expect_error(cca_plot(fc, "heliograph", interactive = TRUE), "not supported for type")
    if (requireNamespace("patchwork", quietly = TRUE)) {
      expect_error(cca_plot(fc, "heliograph_panel", interactive = TRUE), "not supported for type")
    }
  }
})
