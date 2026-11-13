# CCAflow function map

## 1. Audit and choose
`cca_audit()` -> `cca_recommend()` -> `cca_prepare()`

## 2. Core linear CCA
`cca_classic()` | `cca_ridge()` | `cca_shrinkage()`

## 3. Variable selection
`cca_sparse()` | `cca_elastic()` | `cca_group_sparse()` -> `cca_stability()`

## 4. Resistance and non-Gaussian association
`cca_robust()` | `cca_mcd()` | `cca_robust_sparse()` | `cca_rank()` | `cca_copula()`

## 5. Nonlinear association
`cca_kernel()` | `cca_robust_kernel()` | `cca_deep()`

## 6. Latent and mixed data
`cca_probabilistic()` | `cca_bayes()` | `cca_mixed()`

## 7. Structured observations
`cca_functional()` | `cca_dynamic()` | `cca_repeated()` | `cca_multiway()`

## 8. Multiblock
`cca_gcca()` | `cca_rgcca()` | `cca_sgcca()` | `cca_kernel_gcca()`

## 9. Infer, validate and compare
`cca_infer()` | `cca_permutation()` | `cca_bootstrap()` | `cca_weight_diagnostics()` | `cca_validate()` | `cca_nested_cv()` | `cca_compare()`

## 10. Teach and publish
`cca_teach()` | `cca_tour()` | `cca_tour_render()` | `cca_explorer()` | `cca_plot()` (`suppressor`, `bipolar_lollipop`, `tornado`, `synoptic_heatmap`, `heliograph`, `bipolar_balance`, `heliograph_panel`) | `cca_table()` (`"participation"`) | `cca_export()` | `cca_report(..., extra_figures = TRUE)`

**Interpretation rule:** report canonical correlation together with loadings, cross-loadings, redundancy, uncertainty/stability and out-of-sample association. Do not select an estimator by training correlation alone.
