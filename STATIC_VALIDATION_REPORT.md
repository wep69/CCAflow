# CCAflow 0.1.0 static validation report

**Validation date:** 2026-08-14  
**Package:** CCAflow  
**Version:** 0.1.0  
**Validation scope:** source-code and package-structure validation in the current build environment.

## Executive result

The static integration audit completed with **no detected structural failures**. All public exports resolve to function definitions, all registered S3 methods have definitions, all exported functions are covered by the shipped manual aliases, all R and Stan source files passed delimiter-balance checks, and every package referenced with `pkg::fun` syntax is declared in `Imports` or `Suggests`.

This is **not** a substitute for `R CMD check`. The current build environment does not contain `R` or `Rscript`, so R parsing, package loading, numerical tests, vignettes, Stan compilation, optional-engine execution, and `R CMD build/check` could not be executed here.

## Package inventory

| Item | Count |
|---|---:|
| Public exported functions | 70 |
| Registered S3 methods | 28 |
| R source modules | 11 |
| Manual `.Rd` pages | 9 |
| Test files | 6 |
| Vignettes | 6 |
| Bibliographic entries | 14 |
| Stan models | 1 |

## Static checks completed

- Balanced parentheses, brackets, and braces in all 11 R source files after excluding comments and quoted strings.
- Balanced delimiters in `inst/stan/pcca.stan`.
- Every symbol exported in `NAMESPACE` has a detected function definition.
- Every registered S3 method has a detected implementation.
- Every exported function has a corresponding manual alias.
- All nine `.Rd` files have balanced braces.
- Required `DESCRIPTION` fields are present.
- All packages referenced with explicit namespace syntax are declared in `Imports` or `Suggests`.
- No duplicate public/internal function definitions were detected by the static scanner.
- Required package paths are present, including `R/`, `man/`, `tests/`, `vignettes/`, `inst/stan/`, report template, examples, and citation files.
- No `TODO` or `FIXME` markers remain in the source tree.
- `NAMESPACE` and `.Rd` documentation are intentionally maintained manually in this development release; the local protocol therefore uses `devtools::check(..., document = FALSE)`.

## Scientific integration checks represented in the source

The source implements or provides explicit adapters for:

- classical CCA, coefficients, loadings, cross-loadings, variance extracted, redundancy, and classical sequential inference;
- permutation and bootstrap inference, including component order/sign alignment across resamples;
- ridge, shrinkage, sparse, elastic-net-inspired, grouped sparse, rank/copula, and resistant variants;
- nonlinear/kernel CCA and a resistant-preprocessed kernel workflow;
- latent Gaussian probabilistic CCA by EM and Bayesian pCCA through an included Stan model;
- user-configured Deep CCA and a direct CCA-Zoo adapter through `reticulate`;
- mixed-data CCA, functional CCA, lagged/dynamic workflows, repeated-measures reshaping, multiway workflows, GCCA/RGCCA/SGCCA/TGCCA, kernel/multiple CCA, and multiblock tuning/validation adapters;
- ordinary repeated k-fold validation, grouped validation, rolling-origin temporal validation, nested cross-validation, stability selection, tuning, and method comparison;
- a teaching laboratory with simulations, agronomic examples, canonical/guided/grand tours, tour keyframes, and an exploratory four-panel display;
- publication-oriented plots, vector/raster export, formatted tables, and automated report generation.

## Checks that still require a local R installation

The following must be completed before scientific release, submission to CRAN, or use as the sole computational basis of a manuscript:

1. Parse and load all R source files with the target R version.
2. Run `devtools::test()` and confirm the complete test suite passes.
3. Render all six vignettes and the report template.
4. Run the numerical cross-check against `stats::cancor()` included in `LOCAL_VALIDATION.md`.
5. Exercise optional engines individually: CCA/CCP, PMA, RGCCA, mixedCCA, multiCCA, spinifex/tourr, rstan/posterior, and reticulate/CCA-Zoo.
6. Compile and sample the included Stan pCCA model.
7. Inspect publication figures and exported tables visually.
8. Run `R CMD build CCAflow` and `R CMD check --as-cran CCAflow_0.1.0.tar.gz`.
9. Repeat checks on current Windows and Linux R releases; add macOS/R-devel validation when preparing a public release.
10. Replace the placeholder maintainer name/email in `DESCRIPTION` before public submission.

## Runtime detected in the current build environment

- `R`: not available
- `Rscript`: not available
- `pandoc`: `/usr/bin/pandoc`

## Release interpretation

CCAflow 0.1.0 should therefore be treated as a **comprehensive development source package that passed static integration checks**, not as a CRAN-validated binary or an `R CMD check`-certified release. Follow `LOCAL_VALIDATION.md` on the intended workstation before publication-grade use.
