# CCAflow 0.2.0

- Added `cca_weight_diagnostics()` (with `print()` and `plot()` S3 methods) to compare standardized weights against canonical loadings and cross-loadings, quantify block VIF (with automatic ridge regularization when collinear or high-dimensional), and flag sign inversions and suppressor variables.
- Added seven bipolar weight-participation visualizations to `cca_plot()` (`"suppressor"`, `"bipolar_lollipop"`, `"weight_vs_loading"` / `"tornado"`, `"synoptic_heatmap"`, `"heliograph"`, `"bipolar_balance"`, and `"heliograph_panel"` via optional `patchwork`), plus optional label collision avoidance via `ggrepel`.
- Extended `cca_plot()` for `"coefficients"`, `"sparsity"`, and `cca_bootstrap` objects with `components`, `block = c("x", "y", "both")`, `style = c("dot", "lollipop")`, `polarity_color`, `sparsity_mode`, `highlight_zero`, `variable_labels`, `vif_cap`, `squish`, `order_by_component`, and `metric` arguments.
- Added `cca_table(fit, "participation")` and `cca_report(..., extra_figures = FALSE)` (`extra_figures = TRUE` exports `bipolar_lollipop`, `suppressor`, and `synoptic_heatmap`).
- Enforced explicit guardrails rejecting kernel models for original-variable weight diagnostics and rejecting `interactive = TRUE` for polar heliograph plots.

# CCAflow 0.1.1

Corrections reported by the audited tutorial study (see the study's
RELATORIO-AO-AUTOR.md). Each item is covered by a regression test.

## Silent failures

* `cca_validate()`, `cca_validate_time()`, `cca_validate_grouped()`, `cca_tune()`,
  `cca_nested_cv()` and `cca_compare()` now reject an unknown `method` name instead
  of returning an empty but legitimate-looking result, and warn when no fold
  produced a usable fit. `cca_compare()` names the method it recorded as `NA`.
* `cca_kernel()` warns when a canonical correlation lands within `1e-4` of the
  ceiling, which with a wide radial kernel indicates memorisation; the message names
  the observed value. `gamma_x`/`gamma_y` are documented accordingly.
* `cca_kernel_gcca()` requires the documented repeated-measures block form (a list
  with one matrix per subject) instead of silently reading `n` as rows x columns.

## Reproducibility and audit trail

* `cca_robust()`, `cca_mcd()`, `cca_robust_sparse()` and `cca_robust_kernel()`
  accept `seed`, because the MCD engines draw random subsets. The seed is recorded
  in `tuning` and in `refit_params`.
* `cca_sparse(engine = "auto")` reports which engine it selected and warns that the
  penalty scale differs between engines.
* `cca_rgcca()` records the effective `tau` and the effective `ncomp`, and warns when
  the requested number of components is reduced by the block rank.

## Outputs that did not happen

* `plot(<cca_bootstrap>, quantity = "loading" | "cross_loading" | "weight")` now
  draws: the interval layer used column names produced by named `quantile()` values
  (`lower.2.5.`), which no longer happens.
* `cca_plot(fit, "scree", bootstrap = bt)` and `cca_report(bootstrap = bt,
  include_figures = TRUE)` now produce their figure and their file: the interval
  layer no longer inherits the base plot's `correlation` aesthetic.
* `cca_plot()` serves inference and validation objects (bootstrap, validation,
  permutation, comparison, tuning, stability, robust comparison) through the same
  figure helpers as their `plot()` methods, so the documented figure kinds are
  reachable from one entry point without a second round trip.
* `cca_plot()` gains `type = "communality"` for fitted objects, the one figure kind
  that had no route, and an inference key applied to a fit now points at the object
  to build first. A companion object (`bootstrap =`, `permutation =`) that the
  requested type cannot use warns instead of being ignored in silence.
* `cca_export()` refuses plotly widgets with an explicit message instead of an
  obscure device error.
* `cca_grand_tour()` gains `print()` and `plot()` methods; the plot method points to
  `cca_tour()`/`cca_tour_render()` because a grand tour has no fixed frame set.

## Interface and messages

* `predict()` methods require paired blocks with the same number of rows and name the
  expected number of variables when the block does not match.
* `cca_robust_sparse()` and the sparse family report the penalties when sparsity
  removes every weight, instead of a generic 'no finite component' error.
* `cca_repeated()` names missing columns; `cca_tour_keyframes()` warns when
  `strategy = "objective"` has no objective trace to use.
* `cca_multiblock_cv()` no longer leaks the engine's `extra argument 'connection'`
  warning, and `cca_multiblock_stability(type = "stability")` explains that a
  sparse fit from `cca_sgcca()` is required.
* Deprecated `geom_errorbarh()` replaced by `geom_errorbar(orientation = "y")`.

# CCAflow 0.1.0

Initial comprehensive development release.

- Added a common object model for two-block CCA, interpretation quantities, diagnostics, inference, resampling, validation, teaching, visualization, and reporting.
- Added classical, ridge, shrinkage, sparse, elastic, group-sparse, resistant, rank, copula, kernel, mixed-data adapter, probabilistic, Bayesian adapter, deep adapter, functional, dynamic, repeated, multiway, GCCA/RGCCA/SGCCA/TGCCA, and multiple kernel/functional workflows.
- Added an optional CCA-Zoo bridge through `reticulate` for already configured linear, nonparametric/kernel, and probabilistic Python estimators; Python environments are never installed or modified automatically.
- Added training-versus-generalization safeguards, nested cross-validation, time-respecting validation, group-preserving validation, cluster bootstrap, and time-aware circular-shift permutation.
- Added canonical tours, guided tours, keyframe extraction, interactive exploration, and publication-ready figures/tables.
- Added six vignettes, automated tests, a local validation protocol, references, reporting template, and cheatsheet blueprint.
