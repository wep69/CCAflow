# CCAflow technical architecture

## Design objective

CCAflow provides one coherent workflow across classical CCA, high-dimensional regularization, sparsity, resistant estimation, nonlinear association, latent models, functional/temporal data, multiblock integration, teaching, validation, and publication.

## Object model

The core `cca_fit` object stores:

- method and engine;
- original X and Y blocks;
- preprocessing centers and scales;
- standardized and raw canonical weights;
- canonical scores;
- canonical correlations;
- loadings and cross-loadings;
- variance extracted, redundancy, and communalities;
- tuning parameters;
- optional covariance/influence information.

Specialized classes inherit from `cca_fit`, for example `cca_kernel` and `cca_robust`, so generic methods can reuse the common interface.

## Modules

1. `000-utils.R`: validation, whitening, SVD, metric construction, dispatch.
2. `core.R`: classical, ridge, shrinkage, prediction, print and summary.
3. `sparse-robust.R`: sparse, elastic, group sparse, resistant, rank and copula methods.
4. `kernel.R`: kernels, nonlinear dispatcher, resistant kernel preprocessing.
5. `extended.R`: mixed, probabilistic, Bayesian, deep, functional, dynamic, multiblock, repeated, and multiway methods.
6. `inference-diagnostics.R`: audit, diagnostic, Wilks/Bartlett, permutation, bootstrap and stability.
7. `validation.R`: repeated k-fold, rolling-origin and group-preserving validation, nested tuning, model comparison and recommendation.
8. `plots-tables.R`: publication-ready visualization, interactive conversion, table output and export.
9. `teaching-tours.R`: simulations, agronomy teaching data, learning lab and projection tours.
10. `report-version.R`: Markdown reports and engine-version audit.

## Engine policy

Core functionality is dependency-light. Optional engines are called only after explicit availability checks. No optional dependency is installed automatically.

- `PMA`: penalized matrix decomposition sparse CCA.
- `robustbase`/`MASS`: resistant covariance estimation.
- `mixedCCA`: semiparametric sparse mixed-data CCA.
- `RGCCA`: GCCA, RGCCA, SGCCA, TGCCA, permutation tuning, supervised cross-validation, bootstrap and stability.
- `multiCCA`: repeated-measures multiple kernel and functional CCA.
- `rstan`: Bayesian latent shared-factor model.
- `reticulate`: bridge to already configured CCA-Zoo linear/nonparametric/probabilistic estimators through `cca_zoo()`, and to explicitly configured deep models through `cca_deep()`.
- `tourr` and `spinifex`: dynamic projection visualization.
- `ggplot2`, `plotly`, `flextable`, `officer`: publication output.

## Statistical separation of tasks

CCAflow intentionally distinguishes:

- fit;
- inference;
- stability;
- generalization.

This prevents model selection based only on training correlation and makes regularized/nonlinear models directly comparable on out-of-sample criteria.

## Extension points

New CCA engines should return either a full `cca_fit` object or an external-engine wrapper with a reproducible reference to the engine fit. New core methods should implement:

1. preprocessing;
2. canonical directions;
3. scores;
4. canonical correlations;
5. `.cca_make_fit()` integration;
6. prediction when out-of-sample evaluation is meaningful;
7. tests against known invariants and simulation truth.
