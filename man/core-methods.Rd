\name{cca_fit}
\alias{cca_fit}\alias{cca_classic}\alias{cca_ridge}\alias{cca_shrinkage}
\alias{cca_sparse}\alias{cca_elastic}\alias{cca_group_sparse}\alias{cca_rank}\alias{cca_copula}\alias{cca_nonlinear}
\alias{cca_robust}\alias{cca_mcd}\alias{cca_robust_sparse}\alias{cca_robust_kernel}\alias{cca_kernel}
\title{Fit canonical correlation models}
\description{Fit classical and modern two-block canonical correlation estimators using a common object interface.}
\usage{
cca_fit(X, Y,
  method = c("classic", "ridge", "shrinkage", "sparse", "elastic",
    "group_sparse", "robust", "mcd", "rank", "copula", "kernel",
    "robust_sparse", "robust_kernel", "probabilistic", "bayes",
    "mixed", "deep"),
  ...)
cca_classic(X, Y, ncomp = NULL, center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"), tol = 1e-10)
cca_ridge(X, Y, lambda_x = 0.1, lambda_y = 0.1, ncomp = NULL,
  center = TRUE, scale = TRUE, na_action = c("fail", "complete"),
  tol = 1e-10)
cca_shrinkage(X, Y, alpha_x = NULL, alpha_y = NULL, ncomp = NULL,
  engine = c("internal", "corpcor", "auto"), center = TRUE,
  scale = TRUE, na_action = c("fail", "complete"), tol = 1e-10)
cca_sparse(X, Y, ncomp = 1, penalty_x = 0.1, penalty_y = 0.1,
  engine = c("internal", "PMA", "auto"), center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"), max_iter = 200, tol = 1e-7)
cca_elastic(X, Y, ncomp = 1, l1_x = 0.1, l1_y = 0.1,
  l2_x = 0.1, l2_y = 0.1, center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"), max_iter = 200, tol = 1e-7)
cca_group_sparse(X, Y, groups_x, groups_y, ncomp = 1,
  penalty_x = 0.1, penalty_y = 0.1, center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"), max_iter = 200, tol = 1e-7)
cca_nonlinear(X, Y, method = c("kernel", "rank", "copula", "deep"), ...)
cca_rank(X, Y, ncomp = NULL, center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"), tol = 1e-10)
cca_copula(X, Y, ncomp = NULL, center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"), tol = 1e-10)
cca_kernel(X, Y, kernel = c("rbf", "linear", "polynomial"),
  gamma_x = NULL, gamma_y = NULL, degree = 2, coef0 = 1,
  lambda_x = 0.1, lambda_y = 0.1, ncomp = 2, center = TRUE,
  scale = TRUE, na_action = c("fail", "complete"), tol = 1e-10)
cca_robust(X, Y, ncomp = NULL,
  engine = c("winsor", "robustbase", "MASS", "auto"),
  winsor = 4.5, tol = 1e-10, seed = NULL)
cca_mcd(X, Y, ncomp = NULL, engine = c("auto", "robustbase", "MASS"),
  tol = 1e-10, seed = NULL)
cca_robust_sparse(X, Y, ncomp = 1, penalty_x = 0.1, penalty_y = 0.1,
  robust_engine = c("winsor", "robustbase", "MASS", "auto"),
  winsor = 4.5, max_iter = 200, tol = 1e-7, seed = NULL)
cca_robust_kernel(X, Y, kernel = c("rbf", "linear", "polynomial"),
  gamma_x = NULL, gamma_y = NULL, degree = 2, coef0 = 1,
  lambda_x = 0.1, lambda_y = 0.1, ncomp = 2,
  robust_engine = c("winsor", "robustbase", "MASS", "auto"),
  winsor = 4.5, tol = 1e-10, seed = NULL)
}
\arguments{
\item{X,Y}{Paired numeric matrices or data frames. Categorical variables must be encoded explicitly or handled by a method designed for mixed data.}
\item{method}{Estimator name for the high-level dispatcher.}
\item{ncomp}{Requested number of canonical components.}
\item{center,scale}{Logical preprocessing controls.}
\item{na_action}{Missing-data handling: fail with an error or restrict to complete observations.}
\item{lambda_x,lambda_y,alpha_x,alpha_y,penalty_x,penalty_y,l1_x,l1_y,l2_x,l2_y}{Method-specific regularization parameters.}
\item{kernel}{Kernel family for kernel CCA.}
\item{engine}{Numerical engine: internal implementations or optional external engines.}
\item{groups_x,groups_y}{Group identifiers, one per variable in each block, for group-structured sparse CCA.}
\item{gamma_x,gamma_y,degree,coef0}{Kernel hyperparameters (RBF scale, polynomial degree and intercept). A large RBF \code{gamma} drives the Gram matrix towards the identity, so the fit interpolates: the training correlation approaches one even when the blocks are independent. Large values therefore require a larger \code{lambda} and an out-of-sample check with \code{\link{cca_validate}}; \code{cca_kernel()} warns when a training correlation lands within \code{1e-4} of the ceiling.}
\item{robust_engine}{Resistant covariance engine: internal winsorization or optional \pkg{robustbase}/\pkg{MASS}.}
\item{winsor}{Winsorization bound for the internal resistant engine.}
\item{max_iter}{Maximum number of iterations for iterative estimators.}
\item{seed}{Optional integer seed. The MCD engines (\pkg{robustbase}, \pkg{MASS}) draw random subsets; supplying a seed makes the resistant fit reproducible and records it in \code{tuning}. The internal winsor engine is deterministic.}
\item{tol}{Numerical tolerance.}
\item{...}{Additional method-specific arguments.}
}
\value{An object inheriting from \code{cca_fit}, containing canonical scores, correlations, weights when defined, loadings, cross-loadings, variance extracted, redundancy and preprocessing information.}
\details{Training canonical correlation is not a model-selection criterion by itself. Use \code{cca_validate}, \code{cca_stability} and, when tuning, \code{cca_nested_cv}.}
\examples{
## Synthetic and teaching datasets (fast, reproducible)
d_linear   <- cca_simulate(n = 80, scenario = "linear",      seed = 1)
d_agri     <- cca_agri_data(seed = 2)
d_highdim  <- cca_simulate(n = 70, p = 50, q = 40, scenario = "highdim",    seed = 3)
d_collin   <- cca_simulate(n = 80, scenario = "collinear",   seed = 4)
d_sparse   <- cca_simulate(n = 80, scenario = "sparse",      seed = 5)
d_outliers <- cca_simulate(n = 80, scenario = "outliers",    seed = 6)
d_heavy    <- cca_simulate(n = 80, scenario = "heavy_tail",  seed = 7)
d_nonmono  <- cca_simulate(n = 80, scenario = "nonmonotone", seed = 8)

## Base R datasets (numeric columns only)
iris_X   <- iris[c("Sepal.Length", "Sepal.Width")]
iris_Y   <- iris[c("Petal.Length", "Petal.Width")]
mtcars_X <- mtcars[c("mpg", "disp", "hp", "wt")]
mtcars_Y <- mtcars[c("drat", "qsec", "gear")]

## cca_fit: dispatcher over three estimators on three datasets
fit_dispatch_classic <- cca_fit(d_linear$X, d_linear$Y, method = "classic", ncomp = 2)
fit_dispatch_classic$correlations
fit_dispatch_ridge <- cca_fit(d_highdim$X, d_highdim$Y, method = "ridge",
                              lambda_x = 0.2, lambda_y = 0.2, ncomp = 2)
fit_dispatch_ridge$correlations
fit_dispatch_kernel <- cca_fit(d_nonmono$X, d_nonmono$Y, method = "kernel",
                               kernel = "rbf", ncomp = 2)
fit_dispatch_kernel$correlations

## cca_classic: low-dimensional linear associations
fit_classic_linear <- cca_classic(d_linear$X, d_linear$Y, ncomp = 2)
fit_classic_linear$correlations
fit_classic_agri <- cca_classic(d_agri$X, d_agri$Y, ncomp = 2)
fit_classic_agri$correlations
fit_classic_iris <- cca_classic(iris_X, iris_Y, ncomp = 2)
fit_classic_iris$correlations
fit_classic_mtcars <- cca_classic(mtcars_X, mtcars_Y, ncomp = 2)
fit_classic_mtcars$correlations

## cca_ridge: high-dimensional and ill-conditioned blocks
fit_ridge_highdim <- cca_ridge(d_highdim$X, d_highdim$Y,
                               lambda_x = 0.2, lambda_y = 0.2, ncomp = 2)
fit_ridge_highdim$correlations
fit_ridge_collin <- cca_ridge(d_collin$X, d_collin$Y,
                              lambda_x = 0.2, lambda_y = 0.2, ncomp = 2)
fit_ridge_collin$correlations
fit_ridge_mtcars <- cca_ridge(mtcars_X, mtcars_Y,
                              lambda_x = 0.2, lambda_y = 0.2, ncomp = 2)
fit_ridge_mtcars$correlations

## cca_shrinkage: internal shrinkage engine (no external package required)
fit_shrink_highdim <- cca_shrinkage(d_highdim$X, d_highdim$Y, ncomp = 2)
fit_shrink_highdim$correlations
fit_shrink_collin <- cca_shrinkage(d_collin$X, d_collin$Y,
                                   alpha_x = 0.2, alpha_y = 0.2, ncomp = 2)
fit_shrink_collin$correlations
fit_shrink_mtcars <- cca_shrinkage(mtcars_X, mtcars_Y, ncomp = 2)
fit_shrink_mtcars$correlations

## cca_sparse: small penalties with the internal engine
fit_sparse_highdim <- cca_sparse(d_highdim$X, d_highdim$Y,
                                 ncomp = 1, penalty_x = 0.1, penalty_y = 0.1)
fit_sparse_highdim$correlations
fit_sparse_sparse <- cca_sparse(d_sparse$X, d_sparse$Y,
                                ncomp = 1, penalty_x = 0.1, penalty_y = 0.1)
fit_sparse_sparse$correlations
fit_sparse_collin <- cca_sparse(d_collin$X, d_collin$Y,
                                ncomp = 1, penalty_x = 0.1, penalty_y = 0.1,
                                engine = "internal")
fit_sparse_collin$correlations

## cca_elastic: l1 + l2 penalties, internal engine
fit_elastic_highdim <- cca_elastic(d_highdim$X, d_highdim$Y,
                                   ncomp = 1, l1_x = 0.1, l1_y = 0.1)
fit_elastic_highdim$correlations
fit_elastic_sparse <- cca_elastic(d_sparse$X, d_sparse$Y, ncomp = 1,
                                  l1_x = 0.1, l1_y = 0.1, l2_x = 0.05, l2_y = 0.05)
fit_elastic_sparse$correlations
fit_elastic_collin <- cca_elastic(d_collin$X, d_collin$Y,
                                  ncomp = 1, l1_x = 0.05, l1_y = 0.05)
fit_elastic_collin$correlations

## cca_group_sparse: grouped variable selection
groups_highdim_x <- rep(1:4, length.out = ncol(d_highdim$X))
groups_highdim_y <- rep(1:4, length.out = ncol(d_highdim$Y))
fit_gs_highdim <- cca_group_sparse(d_highdim$X, d_highdim$Y,
                                   groups_x = groups_highdim_x,
                                   groups_y = groups_highdim_y,
                                   ncomp = 1, penalty_x = 0.1, penalty_y = 0.1)
fit_gs_highdim$correlations
groups_sparse_x <- rep(1:3, length.out = ncol(d_sparse$X))
groups_sparse_y <- rep(1:3, length.out = ncol(d_sparse$Y))
fit_gs_sparse <- cca_group_sparse(d_sparse$X, d_sparse$Y,
                                  groups_x = groups_sparse_x,
                                  groups_y = groups_sparse_y, ncomp = 1)
fit_gs_sparse$correlations
groups_collin_x <- rep(1:2, length.out = ncol(d_collin$X))
groups_collin_y <- rep(1:2, length.out = ncol(d_collin$Y))
fit_gs_collin <- cca_group_sparse(d_collin$X, d_collin$Y,
                                  groups_x = groups_collin_x,
                                  groups_y = groups_collin_y, ncomp = 1)
fit_gs_collin$correlations

## cca_rank: rank-transformed CCA
fit_rank_linear <- cca_rank(d_linear$X, d_linear$Y, ncomp = 2)
fit_rank_linear$correlations
fit_rank_agri <- cca_rank(d_agri$X, d_agri$Y, ncomp = 2)
fit_rank_agri$correlations
fit_rank_iris <- cca_rank(iris_X, iris_Y, ncomp = 2)
fit_rank_iris$correlations
fit_rank_mtcars <- cca_rank(mtcars_X, mtcars_Y, ncomp = 2)
fit_rank_mtcars$correlations

## cca_copula: Gaussian-copula normal-score CCA
fit_copula_linear <- cca_copula(d_linear$X, d_linear$Y, ncomp = 2)
fit_copula_linear$correlations
fit_copula_agri <- cca_copula(d_agri$X, d_agri$Y, ncomp = 2)
fit_copula_agri$correlations
fit_copula_mtcars <- cca_copula(mtcars_X, mtcars_Y, ncomp = 2)
fit_copula_mtcars$correlations

## cca_nonlinear: dispatcher over nonlinear estimators
fit_nonlinear_kernel <- cca_nonlinear(d_linear$X, d_linear$Y,
                                      method = "kernel", kernel = "rbf", ncomp = 2)
fit_nonlinear_kernel$correlations
fit_nonlinear_rank <- cca_nonlinear(d_agri$X, d_agri$Y,
                                    method = "rank", ncomp = 2)
fit_nonlinear_rank$correlations
fit_nonlinear_copula <- cca_nonlinear(mtcars_X, mtcars_Y,
                                      method = "copula", ncomp = 2)
fit_nonlinear_copula$correlations

## cca_robust: winsorized (internal) resistant CCA
fit_robust_outliers <- cca_robust(d_outliers$X, d_outliers$Y, ncomp = 2)
fit_robust_outliers$correlations
fit_robust_heavy <- cca_robust(d_heavy$X, d_heavy$Y, ncomp = 2)
fit_robust_heavy$correlations
fit_robust_agri <- cca_robust(d_agri$X, d_agri$Y, ncomp = 2)
fit_robust_agri$correlations

## cca_mcd: MCD covariance via the robustbase engine (optional package)
if (requireNamespace("robustbase", quietly = TRUE)) {
  fit_mcd_outliers <- cca_mcd(d_outliers$X, d_outliers$Y, ncomp = 2)
  fit_mcd_outliers$correlations
  fit_mcd_heavy <- cca_mcd(d_heavy$X, d_heavy$Y, ncomp = 2)
  fit_mcd_heavy$correlations
  fit_mcd_agri <- cca_mcd(d_agri$X, d_agri$Y, ncomp = 2)
  fit_mcd_agri$correlations
}

## cca_robust_sparse: resistant preprocessing + internal sparse CCA
fit_rs_outliers <- cca_robust_sparse(d_outliers$X, d_outliers$Y,
                                     ncomp = 1, penalty_x = 0.1, penalty_y = 0.1)
fit_rs_outliers$correlations
fit_rs_heavy <- cca_robust_sparse(d_heavy$X, d_heavy$Y,
                                  ncomp = 1, penalty_x = 0.1, penalty_y = 0.1)
fit_rs_heavy$correlations
fit_rs_agri <- cca_robust_sparse(d_agri$X, d_agri$Y,
                                 ncomp = 1, penalty_x = 0.1, penalty_y = 0.1)
fit_rs_agri$correlations

## cca_robust_kernel: resistant preprocessing + kernel CCA
fit_rk_outliers <- cca_robust_kernel(d_outliers$X, d_outliers$Y,
                                     kernel = "rbf", ncomp = 2)
fit_rk_outliers$correlations
fit_rk_heavy <- cca_robust_kernel(d_heavy$X, d_heavy$Y,
                                  kernel = "rbf", ncomp = 2)
fit_rk_heavy$correlations
fit_rk_agri <- cca_robust_kernel(d_agri$X, d_agri$Y,
                                 kernel = "rbf", ncomp = 2)
fit_rk_agri$correlations

## cca_kernel: linear, RBF and polynomial kernels
fit_kernel_nonmono <- cca_kernel(d_nonmono$X, d_nonmono$Y, kernel = "rbf", ncomp = 2)
fit_kernel_nonmono$correlations
fit_kernel_linear <- cca_kernel(d_linear$X, d_linear$Y, kernel = "linear", ncomp = 2)
fit_kernel_linear$correlations
fit_kernel_agri <- cca_kernel(d_agri$X, d_agri$Y,
                              kernel = "polynomial", degree = 2, ncomp = 2)
fit_kernel_agri$correlations
}
