# CCAflow

**CCAflow** is an integrated R workflow for canonical correlation analysis (CCA) and modern multiview association methods. It is designed for teaching, applied research, advanced modeling, and publication-ready reporting.

The package deliberately separates four questions that are often conflated in routine CCA:

1. **Association fit:** how strongly are the fitted canonical variates associated in the training data?
2. **Association inference:** is there evidence that canonical dimensions carry cross-block association?
3. **Association stability:** are weights, selected variables, loadings, and correlations stable under resampling?
4. **Association generalization:** how much canonical association remains in observations not used for model fitting?

A large in-sample canonical correlation is therefore never treated as sufficient evidence of a useful model.

## Installation

The development version can be installed from GitHub with `devtools`:

```r
# install.packages("devtools")  # if not already installed
devtools::install_github("wep69/CCAflow")
```

To install the package together with its built vignettes:

```r
devtools::install_github("wep69/CCAflow", build_vignettes = TRUE)
```

The `build_vignettes = TRUE` form additionally requires `knitr`, `rmarkdown`, and the optional engine packages used by the vignettes. After installation, list the available vignettes with:

```r
vignette(package = "CCAflow")
```

## Cheatsheet

A 10-page cheatsheet (`CCAflow_cheatsheet_10pages_en.pdf`) is distributed separately from the package. It is attached to the GitHub release and can be downloaded here:

[Download CCAflow_cheatsheet_10pages_en.pdf](https://github.com/wep69/CCAflow/releases/download/v0.1.1/CCAflow_cheatsheet_10pages_en.pdf)

## Implemented families

### Core and inference

- Classical CCA using a numerically stable whitened SVD formulation.
- Standardized and raw canonical coefficients.
- Structure coefficients and cross-loadings.
- Variance extracted, redundancy, and communalities.
- Sequential Wilks lambda and Bartlett chi-square approximation.
- Optional Rao approximation through `CCP`.
- Permutation tests, including time-respecting permutation.
- Bootstrap confidence intervals and stability summaries.
- Cross-validation and nested cross-validation, with grouped and rolling-origin designs.
- Training versus out-of-sample association, preserving the training orientation so sign reversals in held-out data remain visible.

### High-dimensional and sparse CCA

- Ridge CCA.
- Shrinkage CCA, with an optional `corpcor` engine.
- Sparse CCA, with an optional `PMA` engine and an internal soft-threshold implementation.
- Elastic-net-inspired CCA.
- Group-structured sparse CCA.
- Stability selection.

### Resistant and distribution-aware CCA

- MCD-based resistant CCA through `robustbase` or `MASS`, plus an internal winsorized engine.
- Resistant sparse CCA and resistant kernel CCA workflows.
- Rank CCA and Gaussian-copula normal-score CCA.
- Mixed-data sparse CCA adapter through `mixedCCA` for continuous, binary, and truncated data.

### Nonlinear and advanced latent models

- Kernel CCA with linear, polynomial, and radial-basis kernels.
- Gaussian probabilistic CCA fitted by latent-factor EM, with convergence and likelihood diagnostics.
- Bayesian probabilistic CCA through the bundled Stan model and optional `rstan`.
- Deep CCA adapter for user-configured Python models through `reticulate`, with `model_factory` support so each resampling fold starts from a fresh model.
- A direct CCA-Zoo bridge. CCAflow never installs or modifies Python environments.

### Functional, temporal, repeated, and multiblock data

- Spline-basis functional CCA.
- Lagged/dynamic CCA.
- GCCA, RGCCA, SGCCA, and tensor TGCCA adapters through `RGCCA`.
- Multiple kernel and multiple functional CCA through `multiCCA`.
- Repeated-measures long-to-wide and multiway/tensor workflows.

### Teaching and visualization

- `cca_teach()` with introductory, intermediate, and advanced learning paths.
- Canonical, guided, and grand projection tours with `tourr` and `spinifex` rendering.
- Classical-versus-resistant basis comparison and outlier tours.
- Agronomy teaching data and controlled simulations (`cca_simulate()`, `cca_agri_data()`).

### Publication-ready output

`cca_plot()` produces editable `ggplot2` objects (canonical-correlation profiles, score plots, dual correlation circles, biplots, redundancy plots, cross-block heatmaps, association networks, sparsity maps, weight-versus-loading suppressor maps, bipolar lollipops, tornado plots, synoptic multi-dimensional heatmaps, polar heliographs, bipolar balance charts, bootstrap intervals with zero-crossing highlighting, validation distributions, permutation nulls, influence plots, method comparisons, and hyperparameter-search plots). `cca_weight_diagnostics()` and `cca_table(fit, "participation")` summarize standardized weights, loadings, cross-loadings, variance inflation factors, and suppressor flags, while `cca_report(..., extra_figures = TRUE)` exports these extra diagnostic figures alongside the standard report. `cca_export()` exports vector PDF/SVG or high-resolution raster files, and `cca_table()` produces standardized tables exportable as Markdown, LaTeX, HTML, CSV, or Word. `cca_report()` generates a reproducible analysis report.

## Minimal example

```r
library(CCAflow)

d <- cca_agri_data(seed = 11)
fit <- cca_classic(d$X, d$Y, ncomp = 3)

fit
cca_infer(fit)
cca_table(fit, "loadings")
cca_table(fit, "redundancy")

p1 <- cca_plot(fit, "scree")
p2 <- cca_plot(fit, "scores", component = 1)
p3 <- cca_plot(fit, "circle", component = 1)

cca_export(p1, "cca_scree.pdf", preset = "single")
```

## Stability and generalization

```r
set.seed(22)
bt <- cca_bootstrap(fit, B = 499)
perm <- cca_permutation(fit, B = 499)
cv <- cca_validate(d$X, d$Y, method = "classic", folds = 5, repeats = 5)

cca_plot(bt, "bootstrap")
cca_plot(cv, "validation")
cca_plot(perm, "permutation")
```

## Optional engines

Run:

```r
cca_version_info()
```

to inspect which optional engines are available locally. The core package does not install or invoke Python, Stan, or specialized R engines automatically. Core sparse, shrinkage, and resistant workflows use deterministic internal defaults; external engines such as `PMA`, `corpcor`, `robustbase`, and `MASS` are selected only when requested explicitly (or with `engine = "auto"`).

## Validation status

CCAflow 0.1.0 passed the complete local validation protocol on Windows with R 4.6.0: full unit-test suite, numerical cross-checks against `stats::cancor()`, method-family tests (high-dimensional, resistant, nonlinear, multiblock, Bayesian, deep), resampling and dependence-aware validation, publication outputs, tours, all vignettes, `R CMD build`, and `R CMD check --as-cran` with 0 ERRORs and 0 WARNINGs. Cross-platform checks on the official remote builders passed as well: win-builder on current R-release and R-devel (0 ERRORs, 0 WARNINGs, 2 NOTEs documented as new-submission and example-timing) and mac-builder on macOS arm64 (Status OK: 0 errors, 0 warnings, 0 notes).

CCAflow 0.1.1 carries the corrections of the nine defects reported by an independent, execution-verified audit of the package API (see `NEWS.md`). Every correction has a regression test; on this version the local `R CMD check --as-cran` reports 0 ERRORs, 0 WARNINGs and the same 2 NOTEs, and the test suite reports 105 successes with 0 failures. Remote-builder checks for 0.1.1 are pending; the 0.1.0 results above remain the reference for the shared code base.

## Methodological safeguards

- Do not select a CCA model solely because it maximizes training correlation.
- Do not interpret canonical weights alone. Examine loadings, cross-loadings, redundancy, stability, and out-of-sample performance.
- Do not delete outliers automatically. Compare classical and resistant solutions and document any data-quality decision.
- Use nested resampling when tuning regularization, sparsity, or kernel hyperparameters.
- For temporal data, use time-respecting validation rather than random folds.
- Classical Wilks/Bartlett inference is not automatically valid for regularized, sparse, kernel, or deep solutions; use permutation/bootstrap or method-specific inference there.

## License

MIT. See the `LICENSE` file.
