# CCAflow 0.1.0

Initial comprehensive development release.

- Added a common object model for two-block CCA, interpretation quantities, diagnostics, inference, resampling, validation, teaching, visualization, and reporting.
- Added classical, ridge, shrinkage, sparse, elastic, group-sparse, resistant, rank, copula, kernel, mixed-data adapter, probabilistic, Bayesian adapter, deep adapter, functional, dynamic, repeated, multiway, GCCA/RGCCA/SGCCA/TGCCA, and multiple kernel/functional workflows.
- Added an optional CCA-Zoo bridge through `reticulate` for already configured linear, nonparametric/kernel, and probabilistic Python estimators; Python environments are never installed or modified automatically.
- Added training-versus-generalization safeguards, nested cross-validation, time-respecting validation, group-preserving validation, cluster bootstrap, and time-aware circular-shift permutation.
- Added canonical tours, guided tours, keyframe extraction, interactive exploration, and publication-ready figures/tables.
- Added six vignettes, automated tests, a local validation protocol, references, reporting template, and cheatsheet blueprint.
