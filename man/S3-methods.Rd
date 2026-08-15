\name{predict.cca_fit}
\alias{predict.cca_fit}\alias{predict.cca_rank}\alias{predict.cca_copula}\alias{predict.cca_kernel}\alias{predict.cca_robust_kernel}\alias{predict.cca_deep}\alias{predict.cca_zoo}
\alias{print.cca_fit}\alias{summary.cca_fit}\alias{print.summary.cca_fit}\alias{print.cca_multiblock}\alias{print.cca_audit}\alias{print.cca_inference}
\alias{print.cca_validation}\alias{print.cca_tuning}\alias{print.cca_nested_cv}\alias{print.cca_comparison}\alias{print.cca_tour}\alias{print.cca_teaching_lab}
\alias{plot.cca_fit}\alias{plot.cca_robust_comparison}\alias{plot.cca_bootstrap}\alias{plot.cca_validation}\alias{plot.cca_permutation}\alias{plot.cca_comparison}\alias{plot.cca_tuning}\alias{plot.cca_stability}\alias{plot.cca_tour}
\title{S3 methods for CCAflow objects}
\description{Prediction, printing, summarization and plotting methods used by CCAflow object classes.}
\usage{\method{predict}{cca_fit}(object, newdata_x = NULL, newdata_y = NULL, ...)}
\arguments{\item{object}{A fitted CCAflow object.}\item{newdata_x,newdata_y}{New observations for one or both blocks.}\item{...}{Additional method arguments.}}
\value{Method-dependent scores, summaries, printed output or plot objects.}
\examples{
d  <- cca_agri_data(seed = 1)
ht <- cca_simulate(scenario = "heavy_tail", seed = 2)
out <- cca_simulate(scenario = "outliers", seed = 3)
lin <- cca_simulate(scenario = "linear", seed = 4)
nm <- cca_simulate(scenario = "nonmonotone", seed = 5)
hd <- cca_simulate(scenario = "highdim", seed = 6)

## predict.cca_fit: classic fits, held-out rows ------------------------------
fit <- cca_classic(d$X, d$Y, ncomp = 2)
predict(fit, newdata_x = d$X[1:10, ])
predict(fit, newdata_x = d$X[1:10, ], newdata_y = d$Y[1:10, ])

fit <- cca_classic(iris[, 1:2], iris[, 3:4], ncomp = 2)
predict(fit, newdata_x = iris[1:8, 1:2], newdata_y = iris[1:8, 3:4])

fit <- cca_classic(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2)
predict(fit, newdata_x = mtcars[1:6, 1:3], newdata_y = mtcars[1:6, 4:6])

## predict.cca_rank -----------------------------------------------------------
fit <- cca_rank(ht$X, ht$Y, ncomp = 2)
predict(fit, newdata_x = ht$X[1:10, ], newdata_y = ht$Y[1:10, ])
fit <- cca_rank(out$X, out$Y, ncomp = 2)
predict(fit, newdata_x = out$X[1:10, ], newdata_y = out$Y[1:10, ])
fit <- cca_rank(d$X, d$Y, ncomp = 2)
predict(fit, newdata_x = d$X[1:10, ], newdata_y = d$Y[1:10, ])

## predict.cca_copula ---------------------------------------------------------
fit <- cca_copula(ht$X, ht$Y, ncomp = 2)
predict(fit, newdata_x = ht$X[1:10, ], newdata_y = ht$Y[1:10, ])
fit <- cca_copula(out$X, out$Y, ncomp = 2)
predict(fit, newdata_x = out$X[1:10, ], newdata_y = out$Y[1:10, ])
fit <- cca_copula(lin$X, lin$Y, ncomp = 2)
predict(fit, newdata_x = lin$X[1:10, ], newdata_y = lin$Y[1:10, ])

## predict.cca_kernel ---------------------------------------------------------
fit <- cca_kernel(nm$X, nm$Y, ncomp = 2)
predict(fit, newdata_x = nm$X[1:10, ], newdata_y = nm$Y[1:10, ])
fit <- cca_kernel(lin$X, lin$Y, ncomp = 2)
predict(fit, newdata_x = lin$X[1:10, ], newdata_y = lin$Y[1:10, ])
fit <- cca_kernel(d$X, d$Y, ncomp = 2)
predict(fit, newdata_x = d$X[1:10, ], newdata_y = d$Y[1:10, ])

## predict.cca_robust_kernel --------------------------------------------------
fit <- cca_robust_kernel(out$X, out$Y, ncomp = 2)
predict(fit, newdata_x = out$X[1:10, ], newdata_y = out$Y[1:10, ])
fit <- cca_robust_kernel(nm$X, nm$Y, ncomp = 2)
predict(fit, newdata_x = nm$X[1:10, ], newdata_y = nm$Y[1:10, ])
fit <- cca_robust_kernel(d$X, d$Y, ncomp = 2)
predict(fit, newdata_x = d$X[1:10, ], newdata_y = d$Y[1:10, ])

## predict.cca_deep with pure-R callbacks -------------------------------------
fit_fun <- function(X, Y, ncomp, ...) {
  f <- cca_classic(X, Y, ncomp = ncomp, center = FALSE, scale = FALSE)
  list(scores_x = f$scores_x, scores_y = f$scores_y, model = f)
}
transform_fun <- function(model, X, Y, ...) {
  list(x = X \%*\% model$weights_x, y = Y \%*\% model$weights_y)
}
fit <- cca_deep(d$X, d$Y, fit_fun = fit_fun,
  transform_fun = transform_fun, ncomp = 2)
predict(fit, newdata_x = d$X[1:8, ], newdata_y = d$Y[1:8, ])
fit <- cca_deep(iris[, 1:2], iris[, 3:4], fit_fun = fit_fun,
  transform_fun = transform_fun, ncomp = 2)
predict(fit, newdata_x = iris[1:8, 1:2], newdata_y = iris[1:8, 3:4])
fit <- cca_deep(mtcars[, 1:3], mtcars[, 4:6], fit_fun = fit_fun,
  transform_fun = transform_fun, ncomp = 2)
predict(fit, newdata_x = mtcars[1:6, 1:3], newdata_y = mtcars[1:6, 4:6])

## predict.cca_zoo: skipped when no Python cca-zoo is configured --------------
pyok <- requireNamespace("reticulate", quietly = TRUE) &&
  reticulate::py_available(initialize = FALSE) &&
  reticulate::py_module_available("cca_zoo")
if (pyok) {
  fit <- cca_zoo(d$X, d$Y, model = "CCA", ncomp = 2)
  predict(fit, newdata_x = d$X[1:10, ], newdata_y = d$Y[1:10, ])
  fit <- cca_zoo(iris[, 1:2], iris[, 3:4], model = "CCA", ncomp = 2)
  predict(fit, newdata_x = iris[1:8, 1:2], newdata_y = iris[1:8, 3:4])
  fit <- cca_zoo(mtcars[, 1:3], mtcars[, 4:6], model = "CCA", ncomp = 2)
  predict(fit, newdata_x = mtcars[1:6, 1:3], newdata_y = mtcars[1:6, 4:6])
}

## print.cca_fit --------------------------------------------------------------
print(cca_classic(d$X, d$Y, ncomp = 2))
print(cca_ridge(hd$X, hd$Y, ncomp = 2))
print(cca_kernel(nm$X, nm$Y, ncomp = 2))

## summary.cca_fit ------------------------------------------------------------
summary(cca_classic(d$X, d$Y, ncomp = 2))
summary(cca_ridge(hd$X, hd$Y, ncomp = 2))
summary(cca_kernel(nm$X, nm$Y, ncomp = 2))

## print.summary.cca_fit ------------------------------------------------------
print(summary(cca_classic(d$X, d$Y, ncomp = 2)))
print(summary(cca_ridge(hd$X, hd$Y, ncomp = 2)))
print(summary(cca_kernel(nm$X, nm$Y, ncomp = 2)))

## print.cca_multiblock -------------------------------------------------------
if (requireNamespace("RGCCA", quietly = TRUE)) {
  tb <- cca_simulate(scenario = "three_blocks", seed = 7)$blocks
  print(cca_rgcca(tb, ncomp = 2))

  ib <- list(A = as.matrix(iris[, 1:2]), B = as.matrix(iris[, 3:4]),
             C = as.matrix(scale(iris[, 1:4])[, 1:2]))
  print(cca_sgcca(ib, ncomp = 2))

  mb <- list(A = as.matrix(mtcars[, 1:3]), B = as.matrix(mtcars[, 4:7]),
             C = as.matrix(mtcars[, 8:11]))
  print(cca_rgcca(mb, ncomp = 2))
}

## print.cca_audit ------------------------------------------------------------
print(cca_audit(d$X, d$Y))
print(cca_audit(iris[, 1:2], iris[, 3:4]))
print(cca_audit(mtcars[, 1:3], mtcars[, 4:6]))

## print.cca_inference --------------------------------------------------------
print(cca_infer(cca_classic(d$X, d$Y, ncomp = 2)))
print(cca_infer(cca_classic(iris[, 1:2], iris[, 3:4], ncomp = 2)))
print(cca_infer(cca_classic(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2)))

## print.cca_validation -------------------------------------------------------
print(cca_validate(d$X, d$Y, method = "classic", folds = 3, repeats = 2, seed = 1))
print(cca_validate(hd$X, hd$Y, method = "ridge", params = list(ncomp = 2),
                   folds = 3, repeats = 2, seed = 1))
print(cca_validate(nm$X, nm$Y, method = "kernel", params = list(ncomp = 2),
                   folds = 3, repeats = 2, seed = 1))

## print.cca_tuning -----------------------------------------------------------
grid <- expand.grid(lambda_x = c(0.01, 0.1), lambda_y = c(0.01, 0.1))
print(cca_tune(hd$X, hd$Y, method = "ridge", grid = grid, folds = 3, repeats = 2, seed = 1))
print(cca_tune(d$X, d$Y, method = "ridge", grid = grid, folds = 3, repeats = 2, seed = 1))
print(cca_tune(mtcars[, 1:3], mtcars[, 4:6], method = "ridge", grid = grid,
               folds = 3, repeats = 2, seed = 1))

## print.cca_nested_cv --------------------------------------------------------
print(cca_nested_cv(d$X, d$Y, method = "ridge", grid = grid,
                    outer_folds = 3, inner_folds = 3, inner_repeats = 2, seed = 1))
print(cca_nested_cv(iris[, 1:2], iris[, 3:4], method = "ridge", grid = grid,
                    outer_folds = 3, inner_folds = 3, inner_repeats = 2, seed = 1))
print(cca_nested_cv(mtcars[, 1:3], mtcars[, 4:6], method = "ridge", grid = grid,
                    outer_folds = 3, inner_folds = 3, inner_repeats = 2, seed = 1))

## print.cca_comparison -------------------------------------------------------
print(cca_compare(d$X, d$Y, methods = c("classic", "ridge", "shrinkage"),
                  folds = 3, repeats = 2, seed = 1))
print(cca_compare(iris[, 1:2], iris[, 3:4], methods = c("classic", "kernel", "rank"),
                  folds = 3, repeats = 2, seed = 1))
print(cca_compare(mtcars[, 1:3], mtcars[, 4:6], methods = c("classic", "robust", "copula"),
                  folds = 3, repeats = 2, seed = 1))

## print.cca_tour -------------------------------------------------------------
print(cca_tour(cca_classic(d$X, d$Y, ncomp = 2), frames = 12))
print(cca_tour(cca_classic(iris[, 1:2], iris[, 3:4], ncomp = 2), frames = 12))
print(cca_tour(cca_classic(lin$X, lin$Y, ncomp = 2), frames = 12))

## print.cca_teaching_lab -----------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  print(cca_teach(level = "intro"))
  print(cca_teach(level = "intermediate"))
  print(cca_teach(level = "advanced"))
}

## plot.cca_fit ---------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(cca_classic(d$X, d$Y, ncomp = 2), type = "scree")
  plot(cca_classic(iris[, 1:2], iris[, 3:4], ncomp = 2), type = "scores")
  plot(cca_classic(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2), type = "circle")
}

## plot.cca_robust_comparison -------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(cca_robust_compare(out$X, out$Y))
  plot(cca_robust_compare(ht$X, ht$Y))
  plot(cca_robust_compare(d$X, d$Y))
}

## plot.cca_bootstrap ---------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  bt <- cca_bootstrap(cca_classic(d$X, d$Y, ncomp = 2), B = 19, seed = 1)
  plot(bt)
  plot(bt, quantity = "loading", block = "x")
  plot(bt, quantity = "cross_loading", block = "y")

  plot(cca_bootstrap(cca_classic(iris[, 1:2], iris[, 3:4], ncomp = 2), B = 19, seed = 1))
  plot(cca_bootstrap(cca_classic(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2), B = 19, seed = 1))
}

## plot.cca_validation --------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(cca_validate(d$X, d$Y, method = "classic", folds = 3, repeats = 2, seed = 1))
  plot(cca_validate(iris[, 1:2], iris[, 3:4], method = "ridge", params = list(ncomp = 2),
                    folds = 3, repeats = 2, seed = 1))
  plot(cca_validate(mtcars[, 1:3], mtcars[, 4:6], method = "shrinkage",
                    folds = 3, repeats = 2, seed = 1))
}

## plot.cca_permutation -------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(cca_permutation(cca_classic(d$X, d$Y, ncomp = 2), B = 19, seed = 1))
  plot(cca_permutation(cca_classic(iris[, 1:2], iris[, 3:4], ncomp = 2), B = 19, seed = 1))
  plot(cca_permutation(cca_classic(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2), B = 19, seed = 1))
}

## plot.cca_comparison --------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(cca_compare(d$X, d$Y, methods = c("classic", "ridge", "shrinkage"),
                   folds = 3, repeats = 2, seed = 1))
  plot(cca_compare(iris[, 1:2], iris[, 3:4], methods = c("classic", "kernel", "rank"),
                   folds = 3, repeats = 2, seed = 1))
  plot(cca_compare(mtcars[, 1:3], mtcars[, 4:6], methods = c("classic", "robust", "copula"),
                   folds = 3, repeats = 2, seed = 1))
}

## plot.cca_tuning ------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(cca_tune(hd$X, hd$Y, method = "ridge", grid = grid, folds = 3, repeats = 2, seed = 1))
  plot(cca_tune(d$X, d$Y, method = "ridge", grid = grid, folds = 3, repeats = 2, seed = 1))
  plot(cca_tune(mtcars[, 1:3], mtcars[, 4:6], method = "ridge", grid = grid,
                folds = 3, repeats = 2, seed = 1))
}

## plot.cca_stability ---------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(cca_stability(cca_sparse(d$X, d$Y, ncomp = 2), B = 19, seed = 1))
  plot(cca_stability(cca_sparse(iris[, 1:2], iris[, 3:4], ncomp = 2), B = 19, seed = 1))
  plot(cca_stability(cca_sparse(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2), B = 19, seed = 1))
}

## plot.cca_tour --------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(cca_tour(cca_classic(d$X, d$Y, ncomp = 2), frames = 12))
  plot(cca_tour(cca_classic(iris[, 1:2], iris[, 3:4], ncomp = 2), frames = 12))
  plot(cca_tour(cca_classic(lin$X, lin$Y, ncomp = 2), frames = 12))
}
}
