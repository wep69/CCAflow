\name{cca_audit}
\alias{cca_audit}\alias{cca_diagnose}\alias{cca_weight_diagnostics}\alias{cca_recommend}\alias{cca_robust_compare}
\alias{cca_set_metadata}\alias{cca_version_info}\alias{cca_prepare}
\title{Data audit, diagnostics, method recommendation and provenance}
\description{Inspect dimensionality, numerical conditioning, missingness, intra-block VIF, weight-vs-loading suppressor inversions and potential influence before interpreting canonical correlations.}
\usage{
cca_audit(X, Y, near_cor = 0.95)
cca_diagnose(object, cutoff = 0.975)
cca_weight_diagnostics(object, components = NULL, min_abs = 0.05,
  weight_inflation = 0.85, labels = NULL, ridge_eps = 1e-4)
cca_recommend(X, Y, nonlinear = FALSE, mixed = FALSE,
  functional = FALSE, multiblock = FALSE, outliers = NULL)
cca_robust_compare(X, Y, ncomp = NULL, robust_engine = "winsor", ...)
cca_set_metadata(object, ids = NULL, groups = NULL,
  supplementary = NULL)
cca_version_info()
cca_prepare(X, Y, center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"))
}
\arguments{
\item{X,Y}{Paired data blocks.}
\item{object}{A fitted object.}
\item{near_cor}{Correlation threshold above which pairs of variables are reported as nearly collinear.}
\item{cutoff}{Quantile of the chi-square reference distribution used to flag observations.}
\item{components}{Integer vector of canonical components to inspect (defaults to up to 4).}
\item{min_abs}{Minimum absolute coefficient magnitude before flagging sign reversal between weight and loading.}
\item{weight_inflation}{Threshold on absolute standardized canonical weight above which collinear inflation is flagged.}
\item{labels}{Optional named character vector mapping variable names to display labels.}
\item{ridge_eps}{Minimal diagonal regularization added when computing intra-block VIF on singular or rank-deficient correlation matrices.}
\item{nonlinear,mixed,functional,multiblock,outliers}{Scientifically motivated characteristics used to generate recommendations.}
\item{ncomp}{Number of canonical components retained in the compared fits.}
\item{robust_engine}{Resistant covariance engine forwarded to \code{cca_robust()}.}
\item{ids,groups,supplementary}{Optional metadata attached to the fitted object.}
\item{center,scale,na_action}{Preprocessing controls for \code{cca_prepare()}.}
\item{...}{Additional method-specific arguments.}
}
\value{A diagnostic/audit object, recommendation list, comparison object, prepared data object, or version table.}
\details{Recommendations are advisory and never select a method solely because it maximizes training correlation.}
\examples{
d  <- cca_agri_data(seed = 1)
out <- cca_simulate(scenario = "outliers", seed = 2)
nm <- cca_simulate(scenario = "nonmonotone", seed = 3)
hd <- cca_simulate(scenario = "highdim", seed = 4)
ht <- cca_simulate(scenario = "heavy_tail", seed = 5)

## cca_audit ----------------------------------------------------------------
print(cca_audit(d$X, d$Y))
print(cca_audit(iris[, 1:2], iris[, 3:4]))
print(cca_audit(mtcars[, 1:3], mtcars[, 4:6]))

## cca_diagnose -------------------------------------------------------------
dg <- cca_diagnose(cca_classic(d$X, d$Y, ncomp = 2))
head(dg$data)
dg <- cca_diagnose(cca_robust(out$X, out$Y))
head(dg$data)
dg <- cca_diagnose(cca_kernel(nm$X, nm$Y, ncomp = 2))
head(dg$data)

## cca_weight_diagnostics ---------------------------------------------------
wd <- cca_weight_diagnostics(cca_classic(d$X, d$Y, ncomp = 2))
print(wd)
head(wd$table)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(wd)
}

fit_hd_wd <- cca_ridge(hd$X, hd$Y, lambda_x = 0.2, lambda_y = 0.2, ncomp = 2)
wd_hd <- cca_weight_diagnostics(fit_hd_wd, components = 1:2)
head(wd_hd$table)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(wd_hd)
}

wd_mt <- cca_weight_diagnostics(
  cca_classic(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2),
  weight_inflation = 0.7
)
print(wd_mt)

## cca_recommend ------------------------------------------------------------
cca_recommend(d$X, d$Y, nonlinear = TRUE)$recommendations
cca_recommend(iris[, 1:2], iris[, 3:4], mixed = TRUE)$recommendations
cca_recommend(hd$X, hd$Y, outliers = TRUE)$recommendations

## cca_robust_compare -------------------------------------------------------
cca_robust_compare(out$X, out$Y)$comparison
cca_robust_compare(ht$X, ht$Y)$comparison
cca_robust_compare(d$X, d$Y)$comparison

## cca_set_metadata ---------------------------------------------------------
fit <- cca_classic(d$X, d$Y, ncomp = 2)
fit <- cca_set_metadata(fit, groups = d$treatment)
head(fit$metadata$groups)

fit <- cca_classic(iris[, 1:2], iris[, 3:4], ncomp = 2)
fit <- cca_set_metadata(fit, ids = paste0("obs", seq_len(nrow(iris))))
head(fit$metadata$ids)

fit <- cca_classic(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2)
fit <- cca_set_metadata(fit, supplementary = data.frame(cyl = mtcars$cyl, hp = mtcars$hp))
head(fit$metadata$supplementary)

## cca_version_info ---------------------------------------------------------
cca_version_info()
vi <- cca_version_info()
vi[vi$available, ]
vi[!is.na(vi$version), c("package", "version")]
n_engines <- sum(cca_version_info()$available)
n_engines

## cca_prepare --------------------------------------------------------------
prep <- cca_prepare(d$X, d$Y)
names(prep)
prep <- cca_prepare(iris[, 1:2], iris[, 3:4], center = FALSE, scale = TRUE)
head(prep$Xs)
prep <- cca_prepare(mtcars[, 1:3], mtcars[, 4:6], center = TRUE, scale = FALSE)
head(prep$Ys)
}
