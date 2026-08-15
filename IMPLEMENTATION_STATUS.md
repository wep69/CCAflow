# CCAflow 0.1.0 implementation status

Date: 2026-08-14

This matrix maps the 44 planned development blocks to concrete package components. "Internal" means the estimator/workflow is implemented in CCAflow without a specialized statistical engine. "Adapter" means CCAflow validates inputs and preserves the external fitted object, but the estimator itself is supplied by the cited optional package. "Workflow" means a cross-cutting analysis, teaching, validation, or publication layer.

| # | Planned capability | Status | Main implementation |
|---:|---|---|---|
| 1 | State of the art | Complete documentation | `STATE_OF_ART.md`, `inst/REFERENCES.bib` |
| 2 | Gap analysis | Complete documentation | `STATE_OF_ART.md`, `ARCHITECTURE.md` |
| 3 | Unified mathematical foundation | Implemented + documented | whitened-SVD covariance solver, vignettes 01–04 |
| 4 | Uniform data architecture | Internal | `cca_prepare()`, common `cca_fit` object, metadata |
| 5 | Automatic audit | Internal | `cca_audit()`, `cca_diagnose()` |
| 6 | Reproducible preprocessing | Internal | stored centers/scales/transformation maps and refit recipes |
| 7 | Classical CCA | Internal | `cca_classic()` |
| 8 | Classical sequential inference | Internal + optional adapter | Wilks/Bartlett in `cca_infer()`, optional `CCP` cross-check |
| 9 | Bootstrap inference | Internal | `cca_bootstrap()`, grouped bootstrap, component order/sign alignment |
| 10 | Canonical coefficients | Internal | standardized and raw weights in `cca_fit` |
| 11 | Structure coefficients | Internal | X/Y loadings |
| 12 | Cross-loadings | Internal | X-to-Y and Y-to-X cross-loadings |
| 13 | Redundancy analysis | Internal | per-component X/Y redundancy |
| 14 | Variance extracted / communality | Internal | stored in every compatible two-block fit |
| 15 | Ridge CCA | Internal | `cca_ridge()` |
| 16 | Shrinkage CCA | Internal + optional adapter | `cca_shrinkage()`, optional `corpcor` |
| 17 | Sparse CCA | Internal + optional adapter | `cca_sparse()`, optional `PMA` |
| 18 | Elastic-net-style CCA | Internal | `cca_elastic()`; explicitly labelled elastic-net-inspired |
| 19 | Group/structured sparse CCA | Internal | `cca_group_sparse()` |
| 20 | Stability selection | Internal | `cca_stability()` |
| 21 | Resistant CCA | Internal + optional adapters | winsorized internal engine; `robustbase`/`MASS` optional |
| 22 | MCD CCA | Adapter/fallback | `cca_mcd()` through `robustbase` or `MASS` |
| 23 | Resistant sparse CCA | Internal workflow | `cca_robust_sparse()` |
| 24 | Rank / copula CCA | Internal | `cca_rank()`, `cca_copula()` with learned out-of-sample maps |
| 25 | Mixed-data CCA | Adapter | `cca_mixed()` through `mixedCCA` |
| 26 | Nonlinear CCA interface | Workflow | `cca_nonlinear()` |
| 27 | Kernel CCA | Internal | linear/polynomial/RBF kernels, kernel-PCA projection + ridge CCA |
| 28 | Resistant kernel CCA | Internal workflow | resistant preprocessing + kernel CCA; not overclaimed as a new operator-level estimator |
| 29 | Probabilistic CCA | Internal | Gaussian latent-factor pCCA fitted by EM |
| 30 | Bayesian CCA | Internal model + adapter | bundled Stan pCCA model through optional `rstan`, posterior rho summaries |
| 31 | Deep CCA | Adapter | `cca_deep()` for user-configured models; fresh-model factory for resampling |
| 32 | Functional CCA | Internal + adapter | spline-basis two-block CCA; multiple functional CCA through `multiCCA` |
| 33 | Dynamic/time-series CCA | Internal workflow | lag alignment, rolling-origin CV, circular-shift null |
| 34 | Generalized CCA | Adapter | `cca_gcca()` using current RGCCA `method="gcca"` |
| 35 | RGCCA / SGCCA | Adapter | `cca_rgcca()`, `cca_sgcca()` |
| 36 | Kernel/functional multiple CCA | Adapter | `cca_kernel_gcca()` and `multiCCA` grid search |
| 37 | Repeated/multiway/tensor structures | Internal + adapter | repeated long-to-wide, explicit matricization, and `cca_tgcca()` without flattening |
| 38 | Model selection | Workflow | `cca_tune()`, `cca_select_components()`, multiblock tuning |
| 39 | Nested cross-validation | Internal | `cca_nested_cv()` |
| 40 | Out-of-sample canonical association | Internal | signed held-out correlations; negative sign reversals are retained |
| 41 | Bootstrap stability | Internal | order/sign component alignment, variable-level draws and intervals |
| 42 | Model comparison | Internal | `cca_compare()` with training, signed validation and sparsity metrics |
| 43 | Didactic + visualization system | Workflow | `cca_teach()`, 15 simulations, agronomy data, tours, explorer, static/interactive plots |
| 44 | Publication and reproducibility layer | Workflow | `cca_plot()`, `cca_table()`, `cca_export()`, `cca_report()`, version audit |

## Additional capabilities added during implementation

- Group-preserving cross-validation and cluster bootstrap for repeated/clustered units.
- Rolling-origin validation with estimable multi-observation assessment windows.
- Tensor GCCA (`TGCCA`) through the current `RGCCA` engine.
- Supervised RGCCA cross-validation via `cca_multiblock_cv()`.
- Optional groupwise normal-theory confidence ellipses on paired canonical-score plots.
- Publication keyframes from stored canonical tour paths.
- Explicit safeguards against categorical-to-integer coercion.
- Fresh-model construction for deep-learning resampling to reduce state leakage.

- CCA-Zoo bridge (`cca_zoo()`) for already configured linear, nonparametric/kernel, and probabilistic Python estimators through `reticulate`; no automatic Python installation.

## Validation status

The source tree has been subjected to static integration checks in the build environment. R and Rscript are not installed in this environment, so numerical tests, vignette execution, Stan compilation, optional-engine runtime tests, `R CMD build`, and `R CMD check --as-cran` must be run locally. See `LOCAL_VALIDATION.md` and the generated static-validation report included with the release bundle.
