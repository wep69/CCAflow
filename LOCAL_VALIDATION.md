# CCAflow 0.1.0 local validation protocol

This protocol is intentionally strict. Run it on a clean local R installation before scientific use, distribution, or CRAN submission.

## 1. Prerequisites

Use R 4.2 or newer. For the full validation matrix, install a current R release plus Rtools on Windows when source packages require compilation.

Create a clean project or temporary library. Do not validate against an uncontrolled user library if reproducibility is important.

```r
R.version.string
.libPaths()
```

## 2. Install core development tools

```r
install.packages(c(
  "devtools", "testthat", "knitr", "rmarkdown",
  "ggplot2", "withr"
))
```

## 3. Optional engine matrix

Install only the engines you intend to test. For a full capability check:

```r
install.packages(c(
  "plotly", "tourr", "spinifex", "robustbase", "MASS",
  "CCA", "CCP", "PMA", "RGCCA", "caret", "mixedCCA", "multiCCA",
  "corpcor", "Matrix", "reticulate", "rstan", "posterior",
  "flextable", "officer", "vdiffr"
))
```

`cca_version_info()` reports which optional engines are visible to CCAflow.

## 4. Unpack and inspect

```r
pkg <- "path/to/CCAflow"
list.files(pkg)
read.dcf(file.path(pkg, "DESCRIPTION"))
```

Confirm that `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, `tests/`, `vignettes/`, `inst/stan/`, and the validation documentation are present.

## 5. Replace publication metadata placeholders

Before CRAN submission, replace the placeholder maintainer identity and email in `DESCRIPTION`, and add repository URL/BugReports fields only after persistent project addresses exist. Do not submit with placeholder metadata.

## 6. Inspect the shipped manual documentation

CCAflow 0.1.0 ships a manually maintained `NAMESPACE` and `.Rd` reference manual. Do **not** run `devtools::document()` as part of routine validation, because that command would regenerate files from roxygen comments and may overwrite the curated grouped documentation.

```r
ns <- readLines(file.path(pkg, "NAMESPACE"), warn = FALSE)
rd <- list.files(file.path(pkg, "man"), pattern = "\\.Rd$", full.names = TRUE)
stopifnot(length(ns) > 0L, length(rd) > 0L)
```

During development, change `NAMESPACE` and the corresponding `.Rd` files together. Verify that every exported symbol is documented before release.

## 7. Install the package

```r
devtools::install(pkg, dependencies = FALSE, upgrade = "never")
library(CCAflow)
cca_version_info()
```

## 8. Core smoke test

```r
d <- cca_agri_data(seed = 1)
fit <- cca_classic(d$X, d$Y, ncomp = 3)
stopifnot(inherits(fit, "cca_fit"))
stopifnot(all(fit$correlations >= 0 & fit$correlations <= 1))
stopifnot(ncol(fit$weights_x) == 3)

cca_infer(fit)
cca_table(fit, "loadings")
cca_table(fit, "redundancy")
```

## 9. Numerical cross-check against base R

For a well-conditioned classical dataset:

```r
x <- scale(as.matrix(d$X))
y <- scale(as.matrix(d$Y))
base <- cancor(x, y)
ours <- cca_classic(x, y, center = FALSE, scale = FALSE)

all.equal(
  abs(base$cor[seq_along(ours$correlations)]),
  abs(ours$correlations),
  tolerance = 1e-6
)
```

The signs of canonical vectors are arbitrary. Compare absolute correlations or align signs before comparing coefficients.

## 10. High-dimensional tests

```r
hd <- cca_simulate(n = 80, scenario = "highdim", seed = 2)
fr <- cca_ridge(hd$X, hd$Y, lambda_x = 0.1, lambda_y = 0.1)
fs <- cca_sparse(hd$X, hd$Y, engine = "internal")

stopifnot(all(is.finite(fr$correlations)))
stopifnot(all(is.finite(fs$correlations)))
```

## 11. Resistant-analysis tests

```r
ro <- cca_simulate(scenario = "outliers", seed = 3)
f0 <- cca_classic(ro$X, ro$Y)
frb <- cca_robust(ro$X, ro$Y)
cca_diagnose(frb)
cca_outlier_tour(ro$X, ro$Y)
```

Do not require the resistant and classical correlations to differ in a particular direction. Confirm instead that the resistant-distance and flagged-observation outputs are coherent and finite.

## 12. Nonlinear tests

```r
nl <- cca_simulate(scenario = "nonmonotone", seed = 4)
fl <- cca_classic(nl$X, nl$Y, ncomp = 1)
fk <- cca_kernel(nl$X, nl$Y, kernel = "rbf", ncomp = 1)

fl$correlations
fk$correlations
```

Use cross-validation when drawing performance conclusions. A higher training correlation alone is not a validation criterion.

## 13. Resampling tests

Use small values first for speed:

```r
bt <- cca_bootstrap(fit, B = 29, seed = 5)
pm <- cca_permutation(fit, B = 29, seed = 5)
cv <- cca_validate(d$X, d$Y, method = "classic", folds = 5, repeats = 2, seed = 5)

stopifnot(nrow(bt$correlation_ci) >= 1)
stopifnot(all(pm$p_value >= 0 & pm$p_value <= 1))
stopifnot(nrow(cv$results) > 0)
```

Then repeat with research-grade replication counts appropriate for the application.

## 14. Nested tuning test

```r
grid <- expand.grid(
  lambda_x = c(0.01, 0.1),
  lambda_y = c(0.01, 0.1)
)

nc <- cca_nested_cv(
  d$X, d$Y,
  method = "ridge",
  grid = grid,
  outer_folds = 3,
  inner_folds = 3,
  seed = 6
)

stopifnot(all(is.finite(nc$test_correlation)))
```

## 15. Graphics tests

```r
p1 <- cca_plot(fit, "scree")
p2 <- cca_plot(fit, "scores")
p3 <- cca_plot(fit, "circle")
p4 <- cca_plot(fit, "heatmap")

cca_export(p1, tempfile(fileext = ".pdf"), preset = "single")
cca_export(p1, tempfile(fileext = ".tiff"), preset = "single", dpi = 600)
```

Visually inspect all plots for clipped labels, unreadable text, overlapping annotations, and distorted scales.

## 16. Table export tests

```r
cca_table(fit, "correlations")
cca_table(fit, "loadings", format = "markdown")
cca_table(fit, "redundancy", format = "latex")
cca_table(fit, "loadings", format = "word", file = tempfile(fileext = ".docx"))
```

Open exported Word/HTML/PDF products and inspect them visually.

## 17. Touring tests

```r
t1 <- cca_tour(fit, "x", frames = 20)
t2 <- cca_guided_tour(fit, "x", frames = 20)
stopifnot(length(t1$bases) == 20)
stopifnot(length(t2$objective) == 20)

# Optional dynamic rendering
# cca_grand_tour(fit, "x")
# cca_tour(fit, "x", renderer = "spinifex", animate = TRUE)
```

## 18. Optional engine tests

### PMA

```r
if (requireNamespace("PMA", quietly = TRUE))
  cca_sparse(d$X, d$Y, engine = "PMA")
```

### mixedCCA

```r
if (requireNamespace("mixedCCA", quietly = TRUE)) {
  md <- cca_simulate(n = 120, scenario = "binary", seed = 7)
  cca_mixed(md$X, md$Y, type_x = "continuous", type_y = "continuous", BICtype = 1)
}
```

Construct dedicated binary/truncated blocks for formal mixed-type validation.

### RGCCA

```r
if (requireNamespace("RGCCA", quietly = TRUE)) {
  mb <- cca_simulate(scenario = "three_blocks", seed = 8)$blocks
  cca_rgcca(mb, ncomp = 2)
  cca_sgcca(mb, ncomp = 2, sparsity = 0.8)
}
```

### multiCCA

Use the repeated-measures block structure required by `multiCCA`:

```r
if (requireNamespace("multiCCA", quietly = TRUE)) {
  set.seed(9)
  n <- 20; Tn <- 10
  blocks <- list(
    lapply(seq_len(n), function(i) matrix(rnorm(Tn * 3), Tn, 3)),
    lapply(seq_len(n), function(i) matrix(rnorm(Tn * 2), Tn, 2)),
    lapply(seq_len(n), function(i) matrix(rnorm(Tn * 4), Tn, 4))
  )
  cca_kernel_gcca(blocks, method = "kernel", ncomp = 2)
  cca_kernel_gcca(blocks, method = "functional", ncomp = 2)
}
```

### Bayesian CCA

Run only after confirming a working Stan toolchain:

```r
if (requireNamespace("rstan", quietly = TRUE)) {
  small <- cca_agri_data(n = 30, seed = 10)
  bf <- cca_bayes(small$X[,1:2], small$Y[,1:2], ncomp = 1, iter = 500, chains = 2)
  print(bf$stanfit)
}
```

Inspect R-hat, effective sample size, divergent transitions, posterior predictive adequacy, and prior sensitivity before scientific interpretation.

### Deep CCA

CCAflow does not install or modify Python environments. First verify the active environment with `cca_version_info()`. For scikit-learn-style CCA-Zoo estimators, validate the direct bridge explicitly:

```r
zfit <- cca_zoo(d$X, d$Y, model = "CCA", module = "linear", ncomp = 2)
predict(zfit, d$X[1:10, , drop = FALSE], d$Y[1:10, , drop = FALSE])
```

CCAflow intentionally does not invent neural-network architectures. Configure a suitable deep CCA architecture in Python and pass a fresh-model factory or explicit fit/transform callbacks to `cca_deep()`.

Confirm Python, PyTorch/Lightning when applicable, CCA-Zoo, device selection, random seeds, and held-out evaluation explicitly.

## 19. Dependence-aware validation

For time-ordered data, validate with rolling origins rather than random folds:

```r
tv <- cca_validate_time(
  d$X, d$Y, method = "ridge",
  params = list(lambda_x = 0.1, lambda_y = 0.1),
  initial = 60, assess = 10, cumulative = TRUE
)
stopifnot(all(tv$results$train_end < tv$results$test_start))
```

For repeated or clustered data, keep complete independent units together:

```r
groups <- rep(seq_len(nrow(d$X) / 3), each = 3)
groups <- groups[seq_len(nrow(d$X))]
gv <- cca_validate_grouped(d$X, d$Y, groups = groups, method = "classic", folds = 4, repeats = 2)
for (sp in gv$splits) stopifnot(length(intersect(sp$train_groups, sp$test_groups)) == 0)
```

Where within-group rows are dependent, validate the grouped bootstrap and compare it with the ordinary row bootstrap. For time-series null inference, assess whether circular-shift permutation is scientifically defensible before using `cca_permutation_time()`.

## 20. Bootstrap forest plots

```r
bt <- cca_bootstrap(fit, B = 49, seed = 12)
plot(bt, quantity = "loading", block = "x", component = 1, top = 15)
plot(bt, quantity = "cross_loading", block = "y", component = 1, top = 15)
```

Visually confirm that interval direction, ordering, labels, and zero reference lines are correct.

## 21. Automated tests

```r
devtools::test(pkg)
```

All tests must pass with no unexpected warnings.

## 22. Full package check

```r
devtools::check(pkg, document = FALSE, cran = TRUE)
```

Also run from a terminal:

```text
R CMD build CCAflow
R CMD check --as-cran CCAflow_0.1.0.tar.gz
```

Resolve every ERROR and WARNING. Investigate every NOTE and document only those that are unavoidable and legitimate.

## 23. Multi-platform validation

Before CRAN submission, validate at minimum on:

- current Windows R release;
- current Linux R release;
- current macOS R release if available;
- R-devel through a service such as win-builder or equivalent CRAN-compatible infrastructure.

## 24. Final scientific validation

Before using CCAflow in a publication:

1. compare classical results to at least one independent implementation;
2. verify regularized/sparse results with simulation scenarios of known truth;
3. inspect bootstrap stability;
4. report out-of-sample association for tuned models;
5. verify that interpretation of weights agrees with loadings/cross-loadings and redundancy;
6. visually inspect every final figure and table;
7. save package, R, optional-engine, and operating-system versions;
8. archive analysis code and random seeds.
