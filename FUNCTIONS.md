# CCAflow public API inventory

Version 0.1.0 development source. The API is grouped by scientific purpose rather than alphabetically.

## Data preparation and audit

- `cca_prepare()` explicit missing-data and preprocessing recipe.
- `cca_audit()` dimensionality, conditioning, finite-value, variance, and collinearity audit.
- `cca_recommend()` method-family recommendations with rationale.
- `cca_set_metadata()` identifiers, groups, and supplementary observation metadata.
- `cca_version_info()` optional-engine availability and versions.

## Core two-block fitting

- `cca_fit()` unified two-block dispatcher.
- `cca_classic()` classical CCA.
- `cca_ridge()` ridge-regularized CCA.
- `cca_shrinkage()` shrinkage CCA.
- `cca_sparse()` sparse CCA.
- `cca_elastic()` elastic-net-inspired CCA.
- `cca_group_sparse()` group-structured sparse CCA.
- `cca_rank()` rank-transformed CCA.
- `cca_copula()` empirical Gaussian-copula CCA.
- `cca_nonlinear()` nonlinear-method dispatcher.

## Resistance and influence

- `cca_robust()` resistant covariance CCA with deterministic winsorized internal default and explicit `robustbase`/`MASS` options.
- `cca_mcd()` MCD CCA convenience wrapper.
- `cca_robust_sparse()` resistant sparse CCA workflow.
- `cca_robust_kernel()` resistant preprocessing plus kernel CCA.
- `cca_robust_compare()` direct classical-versus-resistant comparison.
- `cca_diagnose()` observation-level distance and influence diagnostics.

## Kernel and advanced latent models

- `cca_kernel()` linear, polynomial, or RBF kernel CCA.
- `cca_probabilistic()` Gaussian probabilistic CCA fitted by latent-factor EM, with block-specific noise covariance and convergence diagnostics.
- `cca_bayes()` optional Stan Bayesian probabilistic CCA with posterior canonical-correlation summaries and a posterior-mean `cca_fit` representation.
- `cca_deep()` optional Python/deep CCA adapter with `model_factory` support for leakage-safe resampling.
- `cca_zoo()` - optional CCA-Zoo bridge for linear, nonparametric/kernel, and probabilistic scikit-learn-style Python estimators via `reticulate`.
- `cca_mixed()` optional semiparametric mixed-data sparse CCA adapter.

## Functional, temporal, repeated, and multiway

- `cca_functional()` spline-basis functional CCA.
- `cca_dynamic()` lag-aligned CCA.
- `cca_repeated()` long-to-wide repeated-measures helper.
- `cca_multiway()` explicit tensor matricization workflow.

## Multiblock and multiview

- `cca_gcca()` classical-limit GCCA through RGCCA.
- `cca_rgcca()` regularized generalized CCA adapter.
- `cca_sgcca()` sparse generalized CCA adapter.
- `cca_tgcca()` tensor generalized CCA using the current TGCCA implementation in `RGCCA`.
- `cca_kernel_gcca()` multiple kernel or functional CCA through multiCCA.
- `cca_multiblock_tune()` RGCCA permutation tuning or `multiCCA` kernel/functional grid-search adapter.
- `cca_multiblock_cv()` supervised RGCCA tuning by k-fold or leave-one-out cross-validation.
- `cca_multiblock_stability()` RGCCA stability/bootstrap adapter.

## Inference, uncertainty, stability, and component selection

- `cca_infer()` sequential classical inference with explicit scope safeguards.
- `cca_permutation()` row-exchange permutation inference.
- `cca_permutation_time()` circular-shift time-series null.
- `cca_bootstrap()` row bootstrap for canonical correlations, weights, loadings, cross-loadings, and redundancy.
- `cca_bootstrap_grouped()` cluster bootstrap.
- `cca_stability()` bootstrap variable-selection stability.
- `cca_select_components()` Wilks, permutation, or CV-based dimension selection.

## Validation, tuning, and comparison

- `cca_validate()` repeated random-fold out-of-sample validation.
- `cca_validate_time()` rolling-origin validation.
- `cca_validate_grouped()` group-preserving cross-validation.
- `cca_nested_cv()` nested cross-validation.
- `cca_tune()` hyperparameter grid search based on held-out canonical association.
- `cca_compare()` common-criteria method comparison.

## Teaching, simulation, and touring

- `cca_simulate()` controlled CCA scenarios.
- `cca_agri_data()` synthetic soil-by-plant teaching data.
- `cca_teach()` introductory/intermediate/advanced Canonical Correlation Learning Lab.
- `cca_tour()` canonical projection path.
- `cca_grand_tour()` grand-tour wrapper.
- `cca_guided_tour()` canonical-association guided tour.
- `cca_basis_tour()` basis interpolation.
- `cca_compare_tour()` comparison between two fitted bases.
- `cca_outlier_tour()` classical-versus-resistant geometry under influence.
- `cca_tour_keyframes()` reproducible informative static tour frames.
- `cca_tour_render()` static or optional `spinifex` rendering of stored basis paths.
- `cca_explorer()` four-panel canonical association explorer.

## Publication and reporting

- `cca_plot()` publication-ready plotting dispatcher.
- `theme_cca_publication()` editable journal-style ggplot theme.
- `cca_palette()` accessible color/grayscale palette.
- `cca_export()` vector/high-resolution figure export.
- `cca_table()` standardized publication tables and export.
- `cca_report()` reproducible Markdown/HTML/Word analysis report.

The package also registers prediction, print, summary, plot, tuning, multiblock, validation, and bootstrap S3 methods.

