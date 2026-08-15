\name{cca_infer}
\alias{cca_infer}\alias{cca_permutation}\alias{cca_permutation_time}\alias{cca_bootstrap}\alias{cca_bootstrap_grouped}
\alias{cca_stability}\alias{cca_select_components}\alias{cca_validate}\alias{cca_validate_time}\alias{cca_validate_grouped}
\alias{cca_nested_cv}\alias{cca_tune}\alias{cca_compare}
\title{Inference, stability, tuning and out-of-sample validation}
\description{Separate inferential evidence, resampling stability, hyperparameter selection and held-out generalization for canonical association models.}
\usage{
cca_infer(object, rao = FALSE)
cca_permutation(object, B = 999, seed = NULL, stepwise = TRUE, max_stat = TRUE)
cca_permutation_time(object, B = 999, min_shift = 1, seed = NULL)
cca_bootstrap(object, B = 999, conf = 0.95, seed = NULL)
cca_bootstrap_grouped(object, groups, B = 999, conf = 0.95, seed = NULL)
cca_stability(object, B = 500, threshold = 1e-8, seed = NULL)
cca_select_components(object, method = c("wilks", "permutation", "cv"), alpha = 0.05,
  B = 499, X = NULL, Y = NULL, folds = 5, repeats = 3, seed = NULL,
  cv_conf = 0.95)
cca_validate(X, Y, method = "classic", params = list(), folds = 5, repeats = 5, seed = 1)
cca_validate_time(X, Y, method = "ridge", params = list(), initial, assess = 5,
  skip = 0, cumulative = TRUE)
cca_validate_grouped(X, Y, groups, method = "ridge", params = list(), folds = 5,
  repeats = 3, seed = 1)
cca_tune(X, Y, method, grid, fixed_params = list(), folds = 5, repeats = 3,
  seed = 1, objective = c("first", "mean"))
cca_nested_cv(X, Y, method, grid, fixed_params = list(), outer_folds = 5,
  inner_folds = 4, inner_repeats = 2, seed = 1, objective = c("first", "mean"))
cca_compare(X, Y,
  methods = c("classic", "ridge", "shrinkage", "sparse", "copula",
    "kernel"),
  specs = list(), validate = TRUE, folds = 5, repeats = 3, seed = 1)
}
\arguments{
\item{object}{A fitted CCA object.}
\item{X,Y}{Paired numeric blocks.}
\item{B}{Number of resamples or permutations.}
\item{conf,cv_conf}{Confidence level for bootstrap intervals or the central held-out interval used by the CV component-selection heuristic.}
\item{alpha}{Significance level for inferential component-selection rules.}
\item{rao}{If \code{TRUE}, optionally request a CCP-based Wilks cross-check when \pkg{CCP} is installed.}
\item{stepwise,max_stat}{Permutation controls for sequential classical dimensions and max-statistic multiplicity adjustment.}
\item{min_shift}{Minimum circular shift for time-respecting permutation.}
\item{threshold}{Weight threshold used to define selection in stability summaries.}
\item{method}{Method name.}
\item{params}{Named method parameters.}
\item{folds,repeats,seed}{Cross-validation controls.}
\item{initial,assess}{Rolling-origin training and assessment sizes.}
\item{skip,cumulative}{Spacing between rolling origins and whether the training window expands.}
\item{groups}{Experimental-unit or subject grouping variable that must stay within folds.}
\item{grid}{Data frame of hyperparameter combinations.}
\item{fixed_params}{Method arguments held constant during tuning.}
\item{outer_folds,inner_folds,inner_repeats}{Nested cross-validation controls.}
\item{objective}{Whether tuning scores the first canonical dimension or the mean over available dimensions.}
\item{specs,validate}{Per-method comparison specifications and whether to run held-out validation.}
\item{methods}{Methods to compare.}
}
\value{An inference, permutation, bootstrap, stability, validation, tuning, nested-CV or comparison object.}
\details{Held-out correlations preserve the orientation learned in the training data; a negative validation correlation is therefore retained as evidence of poor directional generalization rather than replaced by its absolute value. Classical Wilks/Bartlett inference is confined to the classical setting. Regularized, sparse and nonlinear methods should preferentially use design-compatible permutation and bootstrap procedures. Random row permutations require exchangeability; time series and grouped designs require their dedicated resampling functions.}
\examples{
## ---- cca_infer: classical asymptotic (Wilks/Bartlett) inference ----
d_agri <- cca_agri_data(seed = 1)
print(cca_infer(cca_classic(d_agri$X, d_agri$Y)))

print(cca_infer(cca_classic(iris[, 1:2], iris[, 3:4])))

print(cca_infer(cca_classic(mtcars[, 1:4], mtcars[, 5:8])))

## ---- cca_permutation: exchangeable row permutations ----
d_agri <- cca_agri_data(seed = 2)
perm_classic <- cca_permutation(cca_classic(d_agri$X, d_agri$Y), B = 19, seed = 1)
print(perm_classic$p_value)

hd <- cca_simulate(n = 60, scenario = "highdim", seed = 3)
perm_ridge <- cca_permutation(cca_ridge(hd$X, hd$Y, lambda_x = 0.1, lambda_y = 0.1,
                                        ncomp = 3), B = 19, seed = 2)
print(perm_ridge$p_value)

nm <- cca_simulate(n = 80, scenario = "nonmonotone", seed = 4)
perm_kernel <- cca_permutation(cca_kernel(nm$X, nm$Y, ncomp = 2), B = 19, seed = 3)
print(perm_kernel$p_value)

## ---- cca_permutation_time: circular-shift permutations for time series ----
X_long <- longley[, 1:3]
Y_long <- longley[, 4:6]
pt_long <- cca_permutation_time(cca_classic(X_long, Y_long), B = 19, min_shift = 1, seed = 1)
print(pt_long$p_value)

aq <- airquality[complete.cases(airquality), ]
pt_aq <- cca_permutation_time(cca_ridge(aq[, c("Ozone", "Solar.R", "Wind")],
                                         aq[, c("Temp", "Month", "Day")], ncomp = 3),
                              B = 19, min_shift = 1, seed = 2)
print(pt_aq$p_value)

set.seed(1)
n_ar <- 60
latent <- as.numeric(arima.sim(model = list(ar = 0.7), n = n_ar))
X_ar <- cbind(x1 = latent, x2 = 0.5 * latent + rnorm(n_ar))
Y_ar <- cbind(y1 = 0.8 * latent + rnorm(n_ar), y2 = 0.4 * latent + rnorm(n_ar))
pt_ar <- cca_permutation_time(cca_ridge(X_ar, Y_ar, ncomp = 2), B = 19, min_shift = 1, seed = 3)
print(pt_ar$p_value)

## ---- cca_bootstrap: case-resampling bootstrap intervals ----
d_agri <- cca_agri_data(seed = 5)
boot_classic <- cca_bootstrap(cca_classic(d_agri$X, d_agri$Y), B = 19, seed = 1)
print(boot_classic$correlation_ci)

hd <- cca_simulate(n = 60, scenario = "highdim", seed = 6)
boot_ridge <- cca_bootstrap(cca_ridge(hd$X, hd$Y, lambda_x = 0.1, lambda_y = 0.1,
                                      ncomp = 3), B = 19, seed = 2)
print(boot_ridge$correlation_ci)

nm <- cca_simulate(n = 80, scenario = "nonmonotone", seed = 7)
boot_kernel <- cca_bootstrap(cca_kernel(nm$X, nm$Y, ncomp = 2), B = 19, seed = 3)
print(boot_kernel$correlation_ci)

## ---- cca_bootstrap_grouped: cluster-resampling bootstrap intervals ----
os <- OrchardSprays
boot_os <- cca_bootstrap_grouped(cca_classic(os[, c("rowpos", "colpos")],
                                             os[, "decrease", drop = FALSE]),
                                 groups = os$treatment, B = 19, seed = 1)
print(boot_os$correlation_ci)

d_agri <- cca_agri_data(seed = 8)
boot_agri <- cca_bootstrap_grouped(cca_classic(d_agri$X, d_agri$Y),
                                   groups = rep(1:40, each = 3), B = 19, seed = 2)
print(boot_agri$correlation_ci)

boot_iris <- cca_bootstrap_grouped(cca_classic(iris[, 1:2], iris[, 3:4]),
                                   groups = iris$Species, B = 19, seed = 3)
print(boot_iris$correlation_ci)

## ---- cca_stability: bootstrap variable-selection stability ----
d_agri <- cca_agri_data(seed = 9)
stab_classic <- cca_stability(cca_classic(d_agri$X, d_agri$Y), B = 19, seed = 1)
print(stab_classic$table)

hd <- cca_simulate(n = 60, scenario = "highdim", seed = 10)
stab_sparse <- cca_stability(cca_sparse(hd$X, hd$Y, penalty_x = 0.3, penalty_y = 0.3),
                             B = 19, seed = 2)
print(stab_sparse$table)

stab_ridge <- cca_stability(cca_ridge(hd$X, hd$Y, lambda_x = 0.1, lambda_y = 0.1,
                                      ncomp = 3), B = 19, seed = 3)
print(stab_ridge$table)

## ---- cca_select_components: three selection strategies ----
d_agri <- cca_agri_data(seed = 11)
sel_wilks <- cca_select_components(cca_classic(d_agri$X, d_agri$Y),
                                   method = "wilks", alpha = 0.05)
print(sel_wilks)

sel_perm <- cca_select_components(cca_classic(mtcars[, 1:4], mtcars[, 5:8]),
                                  method = "permutation", B = 19, seed = 1)
print(sel_perm)

sel_cv <- cca_select_components(cca_classic(iris[, 1:2], iris[, 3:4]),
                                method = "cv", folds = 3, repeats = 2, seed = 1)
print(sel_cv)

## ---- cca_validate: repeated k-fold out-of-sample validation ----
d_agri <- cca_agri_data(seed = 12)
val_classic <- cca_validate(d_agri$X, d_agri$Y, method = "classic",
                            folds = 3, repeats = 2, seed = 1)
print(val_classic)

hd <- cca_simulate(n = 60, scenario = "highdim", seed = 13)
val_ridge <- cca_validate(hd$X, hd$Y, method = "ridge",
                          params = list(lambda_x = 0.1, lambda_y = 0.1),
                          folds = 3, repeats = 2, seed = 1)
print(val_ridge)

nm <- cca_simulate(n = 80, scenario = "nonmonotone", seed = 14)
val_kernel <- cca_validate(nm$X, nm$Y, method = "kernel", params = list(ncomp = 2),
                           folds = 3, repeats = 2, seed = 1)
print(val_kernel)

## ---- cca_validate_time: rolling-origin validation for time-ordered data ----
X_long <- longley[, 1:3]
Y_long <- longley[, 4:6]
vt_long <- cca_validate_time(X_long, Y_long, method = "ridge", initial = 8, assess = 3)
print(vt_long)

aq <- airquality[complete.cases(airquality), ]
vt_aq <- cca_validate_time(aq[, c("Ozone", "Solar.R", "Wind")],
                           aq[, c("Temp", "Month", "Day")],
                           method = "ridge", initial = 80, assess = 10)
print(vt_aq)

set.seed(1)
n_ar <- 60
latent <- as.numeric(arima.sim(model = list(ar = 0.7), n = n_ar))
X_ar <- cbind(x1 = latent, x2 = 0.5 * latent + rnorm(n_ar))
Y_ar <- cbind(y1 = 0.8 * latent + rnorm(n_ar), y2 = 0.4 * latent + rnorm(n_ar))
vt_ar <- cca_validate_time(X_ar, Y_ar, method = "ridge", initial = 20, assess = 5)
print(vt_ar)

## ---- cca_validate_grouped: group-preserving cross-validation ----
os <- OrchardSprays
vg_os <- cca_validate_grouped(os[, c("rowpos", "colpos")], os[, "decrease", drop = FALSE],
                              groups = os$treatment, method = "classic",
                              folds = 3, repeats = 2, seed = 1)
print(vg_os)

d_agri <- cca_agri_data(seed = 15)
vg_agri <- cca_validate_grouped(d_agri$X, d_agri$Y, groups = rep(1:40, each = 3),
                                method = "classic", folds = 3, repeats = 2, seed = 1)
print(vg_agri)

vg_iris <- cca_validate_grouped(iris[, 1:2], iris[, 3:4], groups = iris$Species,
                                method = "classic", folds = 3, repeats = 2, seed = 1)
print(vg_iris)

## ---- cca_nested_cv: nested cross-validation for ridge tuning ----
ridge_grid <- expand.grid(lambda_x = c(0.01, 0.1), lambda_y = c(0.01, 0.1))

hd <- cca_simulate(n = 60, scenario = "highdim", seed = 16)
ncv_hd <- cca_nested_cv(hd$X, hd$Y, method = "ridge", grid = ridge_grid,
                        outer_folds = 3, inner_folds = 3, seed = 1)
print(ncv_hd)

d_agri <- cca_agri_data(seed = 17)
ncv_agri <- cca_nested_cv(d_agri$X, d_agri$Y, method = "ridge", grid = ridge_grid,
                          outer_folds = 3, inner_folds = 3, seed = 1)
print(ncv_agri)

ncv_mtcars <- cca_nested_cv(mtcars[, 1:4], mtcars[, 5:8], method = "ridge",
                            grid = ridge_grid, outer_folds = 3, inner_folds = 3, seed = 1)
print(ncv_mtcars)

## ---- cca_tune: hyperparameter grid search ----
ridge_grid <- expand.grid(lambda_x = c(0.01, 0.1), lambda_y = c(0.01, 0.1))

hd <- cca_simulate(n = 60, scenario = "highdim", seed = 18)
tune_hd <- cca_tune(hd$X, hd$Y, method = "ridge", grid = ridge_grid,
                    folds = 3, repeats = 2, seed = 1)
print(tune_hd)

d_agri <- cca_agri_data(seed = 19)
tune_agri <- cca_tune(d_agri$X, d_agri$Y, method = "ridge", grid = ridge_grid,
                      folds = 3, repeats = 2, seed = 1)
print(tune_agri)

tune_mtcars <- cca_tune(mtcars[, 1:4], mtcars[, 5:8], method = "ridge",
                        grid = ridge_grid, folds = 3, repeats = 2, seed = 1)
print(tune_mtcars)

## ---- cca_compare: compare estimators on common criteria ----
hd <- cca_simulate(n = 60, scenario = "highdim", seed = 20)
cmp_hd <- cca_compare(hd$X, hd$Y, methods = c("classic", "ridge", "sparse"),
                      specs = list(sparse = list(penalty_x = 0.3, penalty_y = 0.3)),
                      validate = TRUE, folds = 3, repeats = 2, seed = 1)
print(cmp_hd)

nm <- cca_simulate(n = 80, scenario = "nonmonotone", seed = 21)
cmp_nm <- cca_compare(nm$X, nm$Y, methods = c("classic", "kernel", "copula"),
                      specs = list(kernel = list(ncomp = 2)),
                      validate = TRUE, folds = 3, repeats = 2, seed = 1)
print(cmp_nm)

d_agri <- cca_agri_data(seed = 22)
cmp_agri <- cca_compare(d_agri$X, d_agri$Y, methods = c("classic", "ridge", "shrinkage"),
                        validate = TRUE, folds = 3, repeats = 2, seed = 1)
print(cmp_agri)
}
