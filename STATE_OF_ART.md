# State of the art and gap addressed by CCAflow

Date of package-design review: 2026-08-14.

## Current ecosystem

The R package `CCA` extends classical CCA and includes regularized CCA. `RGCCA` provides a broad multiblock framework with regularized, sparse, generalized, and tensor GCCA methods, permutation or supervised cross-validation for selected hyperparameters, stability/bootstrap assessment, and prediction. `mixedCCA` implements semiparametric sparse CCA for continuous, binary, and truncated continuous variables using latent Gaussian-copula correlations. `multiCCA`, published on CRAN in 2026, provides multiple kernel CCA and multiple functional CCA for repeated-measures/functional multiblock data, along with tuning and diagnostic tools. Python CCA-Zoo provides reference implementations spanning classical, regularized, sparse, kernel/nonparametric, probabilistic, and deep CCA families, together with model-selection utilities.

Dynamic high-dimensional projections are separately developed in `tourr`, while `spinifex` adds manual tours and a layered rendering system that can consume tour paths.

## Integration gap

The methodological capabilities are distributed across packages with different object systems, tuning conventions, validation practices, output structures, and teaching philosophies. CCAflow addresses the integration problem by placing each engine in the same conceptual workflow:

**audit -> method choice -> fit -> inference -> interpretation -> stability -> validation -> visualization -> report**.

The intended novelty is therefore workflow integration, methodological safeguards, teaching, and publication support rather than claiming novelty for established individual estimators.

## Principles incorporated from modern CCA practice

- High-dimensional CCA requires regularization or structural assumptions.
- Sparse CCA improves interpretability but selected variables should be assessed for stability.
- Resistant CCA can reveal whether a small subset of observations drives the classical solution.
- Kernel and deep methods can model nonlinear dependence but increase the need for rigorous tuning and out-of-sample evaluation.
- Mixed-data methods should model the latent association structure rather than treating binary/zero-inflated measurements as ordinary Gaussian variables by default.
- Multiblock methods generalize the two-view objective and need explicit block-connection and regularization choices.
- Dynamic projection tours add interpretive information unavailable from a single static 2D view.

## Core references

See `inst/REFERENCES.bib`. Metadata and DOI information for the central package/method references were checked against CRAN, journal, PubMed/PMLR/JMLR, or publisher records during package design.
