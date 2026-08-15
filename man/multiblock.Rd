\name{cca_gcca}
\alias{cca_gcca}\alias{cca_rgcca}\alias{cca_sgcca}\alias{cca_tgcca}\alias{cca_kernel_gcca}
\alias{cca_multiblock_tune}\alias{cca_multiblock_cv}\alias{cca_multiblock_stability}
\title{Generalized and multiblock canonical association}
\description{Adapters for generalized, regularized, sparse, kernel and functional multiblock canonical component models.}
\usage{
cca_gcca(blocks, ncomp = 2, connection = NULL, scheme = "factorial",
  scale = TRUE, ...)
cca_rgcca(blocks, ncomp = 2, connection = NULL, tau = "optimal",
  sparsity = 1, scheme = "factorial", scale = TRUE,
  method = "rgcca", ...)
cca_sgcca(blocks, ncomp = 2, sparsity = 0.8, connection = NULL,
  scheme = "factorial", scale = TRUE, ...)
cca_tgcca(blocks, ncomp = 2, connection = NULL, tau = 1, rank = 1,
  mode_orth = 1, separable = TRUE, scheme = "factorial",
  scale = TRUE, ...)
cca_kernel_gcca(blocks, method = c("kernel", "functional"),
  ncomp = 2, eps = 0.001, gamma = NULL, nbasis = 5,
  basis_type = "fourier", argvals = NULL, ...)
cca_multiblock_tune(blocks, engine = c("RGCCA", "multiCCA"),
  par_type = c("tau", "sparsity", "ncomp"), par_value = NULL,
  par_length = 10, n_perms = 50, ncomp = 1, connection = NULL,
  method = "rgcca", multi_method = c("kernel", "functional"),
  eps_grid = c(0.0001, 0.001, 0.01), gamma_grid = NULL,
  nbasis_grid = NULL, ncomp_eval = 2, criterion = c("first", "sumk"),
  basis_type = "fourier", argvals = NULL, ...)
cca_multiblock_cv(blocks, response, method = "rgcca",
  par_type = c("tau", "sparsity", "ncomp"), par_value = NULL,
  par_length = 10, validation = c("kfold", "loo"),
  prediction_model = "lm", metric = NULL, k = 5, n_run = 1,
  n_cores = 1, connection = NULL, ...)
cca_multiblock_stability(object, type = c("bootstrap", "stability"),
  n_boot = 200, ...)
}
\arguments{
\item{blocks}{A list of numeric blocks, or the repeated-measures block representation required by the selected engine.}
\item{ncomp}{Number of multiblock components.}
\item{connection}{Connection design matrix forwarded to \pkg{RGCCA}.}
\item{tau}{Regularization values for RGCCA (scalar, vector, or \code{"optimal"}).}
\item{sparsity}{Sparsity parameter for SGCCA.}
\item{scheme}{RGCCA scheme function (e.g., \code{"factorial"}, \code{"centroid"}, \code{"horst"}).}
\item{scale}{Whether blocks are scaled before the multiblock analysis.}
\item{method}{Selected RGCCA or multiple kernel/functional CCA method.}
\item{rank,mode_orth,separable}{Tensor-GCCA decomposition and covariance-structure controls.}
\item{eps,gamma,nbasis,basis_type,argvals}{Controls forwarded to \pkg{multiCCA} kernel/functional fits and grid searches.}
\item{object}{A fitted multiblock object when required.}
\item{engine}{Tuning engine: RGCCA permutation or multiCCA grid search.}
\item{par_type}{Hyperparameter family tuned by permutation or supervised cross-validation.}
\item{par_value,par_length,n_perms}{Permutation-tuning grid controls.}
\item{multi_method,eps_grid,gamma_grid,nbasis_grid,ncomp_eval,criterion}{multiCCA grid-search controls.}
\item{response}{Index of the supervised response block for RGCCA cross-validation.}
\item{validation,prediction_model,metric,k,n_run,n_cores}{RGCCA cross-validation controls.}
\item{type,n_boot}{Bootstrap or stability selection type and number of resamples.}
\item{...}{Arguments forwarded to \pkg{RGCCA} or \pkg{multiCCA}.}
}
\value{A multiblock wrapper retaining the external-engine fit and method metadata.}
\details{These functions require optional packages and preserve the external fit rather than silently translating unsupported inferential quantities into the two-block interface.}
\examples{
## Multiblock datasets used across the examples
tb <- cca_simulate(n = 60, scenario = "three_blocks", rho = 0.7, seed = 1)$blocks
ir_bl <- list(sepal = as.matrix(iris[, 1:2]), petal = as.matrix(iris[, 3:4]))
mt_bl <- list(engine = as.matrix(mtcars[, c("mpg", "cyl", "disp")]),
              drivetrain = as.matrix(mtcars[, c("hp", "drat", "wt")]),
              performance = as.matrix(mtcars[, c("qsec", "vs", "am", "gear", "carb")]))
us_bl <- list(violent = as.matrix(USArrests[, c("Murder", "Assault")]),
              social = as.matrix(USArrests[, c("UrbanPop", "Rape")]))
ad <- cca_agri_data(n = 60, seed = 1)
ag_bl <- list(soil = as.matrix(ad$X[, 1:3]), plant = as.matrix(ad$Y[, 1:3]))

## Generalized, regularized and sparse GCCA (all require RGCCA)
if (requireNamespace("RGCCA", quietly = TRUE)) {
  fit <- cca_gcca(tb, ncomp = 2)
  print(fit)

  fit <- cca_gcca(ir_bl, ncomp = 2)
  print(fit)

  fit <- cca_gcca(mt_bl, ncomp = 2)
  print(fit)

  fit <- cca_rgcca(tb, ncomp = 2)
  print(fit)

  fit <- cca_rgcca(us_bl, ncomp = 2)
  print(fit)

  fit <- cca_rgcca(ag_bl, ncomp = 2)
  print(fit)

  fit <- cca_sgcca(tb, ncomp = 2, sparsity = 0.8)
  print(fit)

  fit <- cca_sgcca(mt_bl, ncomp = 2, sparsity = 0.8)
  print(fit)

  fit <- cca_sgcca(ag_bl, ncomp = 2, sparsity = 0.8)
  print(fit)
}

## Tensor GCCA: present but skipped when the installed RGCCA lacks the 'tgcca' method
tgcca_ok <- requireNamespace("RGCCA", quietly = TRUE) &&
            exists("tgcca", envir = asNamespace("RGCCA"), inherits = FALSE)
if (tgcca_ok) {
  fit <- cca_tgcca(tb, ncomp = 2)
  print(fit)

  fit <- cca_tgcca(ir_bl, ncomp = 2)
  print(fit)

  fit <- cca_tgcca(us_bl, ncomp = 2)
  print(fit)
}

## Kernel and functional multiblock CCA on repeated measures (requires multiCCA)
if (requireNamespace("multiCCA", quietly = TRUE)) {
  set.seed(1); n <- 20; T_len <- 10
  X <- list(lapply(seq_len(n), function(i) matrix(rnorm(T_len * 3), T_len, 3)),
            lapply(seq_len(n), function(i) matrix(rnorm(T_len * 2), T_len, 2)))
  fit <- cca_kernel_gcca(X, method = "kernel", ncomp = 2, eps = 0.01)
  print(fit)

  set.seed(2); n <- 20; T_len <- 12
  X <- list(lapply(seq_len(n), function(i)
             matrix(sin(seq_len(T_len) / 2) + rnorm(T_len), T_len, 2)),
            lapply(seq_len(n), function(i)
             matrix(cos(seq_len(T_len) / 2) + rnorm(T_len), T_len, 2)))
  fit <- cca_kernel_gcca(X, method = "functional", ncomp = 2, nbasis = 5, basis_type = "fourier")
  print(fit)

  set.seed(3); n <- 20; T_len <- 10
  tt <- seq(0, 1, length.out = T_len)
  X <- list(lapply(seq_len(n), function(i) matrix(sin(2 * pi * tt) + rnorm(T_len), T_len, 3)),
            lapply(seq_len(n), function(i) matrix(cos(2 * pi * tt) + rnorm(T_len), T_len, 2)))
  fit <- cca_kernel_gcca(X, method = "functional", ncomp = 2, nbasis = 5,
                         basis_type = "bspline", argvals = tt)
  print(fit)
}

## Permutation and grid tuning
if (requireNamespace("RGCCA", quietly = TRUE)) {
  perm <- cca_multiblock_tune(tb, engine = "RGCCA", par_type = "tau",
                              par_length = 3, n_perms = 10, ncomp = 1, method = "rgcca")
  perm$stats

  perm <- cca_multiblock_tune(ir_bl, engine = "RGCCA", par_type = "ncomp",
                              par_length = 3, n_perms = 10, ncomp = 1, method = "rgcca")
  perm$stats
}
if (requireNamespace("multiCCA", quietly = TRUE)) {
  set.seed(1); n <- 12; T_len <- 8
  X <- list(lapply(seq_len(n), function(i) matrix(rnorm(T_len * 3), T_len, 3)),
            lapply(seq_len(n), function(i) matrix(rnorm(T_len * 2), T_len, 2)))
  gs <- cca_multiblock_tune(X, engine = "multiCCA", multi_method = "kernel",
                            eps_grid = c(0.01, 0.05), gamma_grid = c(0.5, 1),
                            ncomp_eval = 2, criterion = "first")
  print(gs)
}

## Supervised tuning by cross-validation (response = index of the block to predict)
if (requireNamespace("RGCCA", quietly = TRUE)) {
  cv <- suppressWarnings(cca_multiblock_cv(tb, response = 3, method = "rgcca",
        par_type = "tau", par_length = 3, k = 3, n_run = 1))
  cv$stats

  cv <- suppressWarnings(cca_multiblock_cv(mt_bl, response = 1, method = "rgcca",
        par_type = "tau", par_length = 3, k = 3, n_run = 1))
  cv$stats

  cv <- suppressWarnings(cca_multiblock_cv(ag_bl, response = 2, method = "rgcca",
        par_type = "tau", par_length = 3, k = 3, n_run = 1))
  cv$stats
}

## Bootstrap and sparse stability assessment
if (requireNamespace("RGCCA", quietly = TRUE)) {
  fit <- cca_rgcca(tb, ncomp = 2)
  bs <- cca_multiblock_stability(fit, type = "bootstrap", n_boot = 10)
  bs$stats

  fit <- cca_sgcca(tb, ncomp = 2, sparsity = 0.8)
  st <- cca_multiblock_stability(fit, type = "stability", n_boot = 10)
  print(st)

  fit <- cca_rgcca(ir_bl, ncomp = 2)
  bs <- cca_multiblock_stability(fit, type = "bootstrap", n_boot = 10)
  bs$stats
}
}
