## Test environments

- local: Windows 11 x64, R 4.6.0 (ucrt), Rtools45
- win-builder: R-release (R 4.6.1, Windows Server 2022)
- win-builder: R-devel (r90407 ucrt, Windows Server 2022)
- mac-builder: r-release-macosx-arm64 (Apple M1, macOS Tahoe 26.6, R 4.6.1 Patched, clang 17)

## R CMD check results

- local `R CMD check --as-cran`: 0 ERRORs | 0 WARNINGs | 2 NOTEs
- win-builder R-release: 0 ERRORs | 0 WARNINGs | 2 NOTEs
- win-builder R-devel: 0 ERRORs | 0 WARNINGs | 2 NOTEs
- mac-builder macOS arm64: Status OK (0 errors, 0 warnings, 0 notes)

## Notes

1. New submission. The incoming-feasibility check also flags "multiview",
   "multiblock" and "multiway" as possibly misspelled words; these are
   deliberate, standard terms of the multivariate-analysis literature and are
   used consistently throughout the package documentation.
2. Example timings: the `advanced-methods` help page exceeds the 5 s example
   budget because its Bayesian examples compile the bundled Stan model on
   first use (subsequent calls reuse the compiled model), and the
   `S3-methods` page exercises several multiblock adapters with small
   resampling counts. All heavy examples are wrapped in
   `requireNamespace()` guards and are skipped cleanly when the
   corresponding optional engines are unavailable.

## Downstream dependencies

None. This is a new submission.
