\name{cca_probabilistic}
\alias{cca_probabilistic}\alias{cca_bayes}\alias{cca_deep}\alias{cca_zoo}\alias{cca_mixed}
\alias{cca_functional}\alias{cca_dynamic}\alias{cca_repeated}\alias{cca_multiway}
\title{Latent, mixed, functional, temporal and deep canonical association methods}
\description{Advanced canonical association methods and adapters for data structures beyond ordinary two-block Gaussian CCA.}
\usage{
cca_probabilistic(X, Y, ncomp = 2, center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"), max_iter = 500, tol = 1e-6,
  psi_floor = 1e-6, seed = 1)
cca_bayes(X, Y, ncomp = 2, center = TRUE, scale = TRUE,
  na_action = c("fail", "complete"), chains = 4, iter = 2000,
  warmup = floor(iter / 2), seed = 1, refresh = 0,
  max_posterior_draws = 1000, ...)
cca_deep(X, Y, model = NULL, model_factory = NULL,
  fit_fun = NULL, transform_fun = NULL, ncomp = 2, center = TRUE,
  scale = TRUE, na_action = c("fail", "complete"), ...)
cca_zoo(X, Y, model = "CCA",
  module = c("linear", "nonparametric", "probabilistic"),
  ncomp = 2, model_args = list(), python = NULL, center = TRUE,
  scale = TRUE, na_action = c("fail", "complete"))
cca_mixed(X, Y, type_x = "continuous", type_y = "continuous",
  BICtype = 1, nlamseq = 20, lam.eps = 0.01, ...)
cca_functional(X, Y,
  time_x = seq_len(ncol(as.matrix(if (is.list(X)) X[[1]] else X))),
  time_y = seq_len(ncol(as.matrix(if (is.list(Y)) Y[[1]] else Y))),
  basis_df_x = 6, basis_df_y = 6, degree = 3,
  method = c("ridge", "classic"), ncomp = 2,
  lambda_x = 0.1, lambda_y = 0.1)
cca_dynamic(X, Y, lag = 0L, method = "ridge", ncomp = 2, ...)
cca_repeated(data, id, time, x_vars, y_vars, method = "ridge",
  ncomp = 2, ...)
cca_multiway(X, Y, method = "ridge", ncomp = 2, ...)
}
\arguments{
\item{X,Y}{Paired input blocks, with method-specific structure.}
\item{ncomp}{Number of latent/canonical dimensions.}
\item{center,scale,na_action}{Preprocessing controls shared with the core estimators.}
\item{model}{For \code{cca_deep()}, a configured deep CCA model; for \code{cca_zoo()}, the CCA-Zoo estimator class name.}
\item{module,model_args,python}{CCA-Zoo Python namespace, estimator-constructor arguments, and optional Python executable path.}
\item{model_factory}{Function that constructs a fresh deep model. This is the preferred interface for cross-validation.}
\item{fit_fun,transform_fun}{Custom R callbacks for fitting and out-of-sample transformation. If \code{fit_fun} declares a \code{model} argument, the object created by \code{model_factory} is passed explicitly; the callback must return \code{scores_x} and \code{scores_y} and may return the fitted \code{model}.}
\item{max_iter,tol,psi_floor,seed}{EM controls and seed for \code{cca_probabilistic()}.}
\item{chains,iter,warmup,refresh,max_posterior_draws}{Stan controls for \code{cca_bayes()}.}
\item{type_x,type_y}{Data-type specifications accepted by the optional \pkg{mixedCCA} engine.}
\item{BICtype,nlamseq,lam.eps}{\pkg{mixedCCA} tuning controls.}
\item{time_x,time_y}{Observation-grid positions for the functional variables of each block.}
\item{basis_df_x,basis_df_y,degree}{Spline basis dimension and polynomial degree for functional CCA.}
\item{lambda_x,lambda_y}{Ridge penalties used by spline-basis functional CCA.}
\item{lag}{Temporal lag used to align the blocks.}
\item{method}{Underlying CCA method used after structural representation.}
\item{data,id,time,x_vars,y_vars}{Long-format repeated-measures data specification.}
\item{...}{Additional method-specific arguments.}
}
\value{A \code{cca_fit} object or a compatible advanced wrapper.}
\details{Bayesian fitting requires \pkg{rstan}. Mixed-data fitting requires \pkg{mixedCCA}. The CCA-Zoo bridge requires an already configured Python \code{cca-zoo} environment through \pkg{reticulate} and never installs Python packages automatically. Deep models are user-configured; use \code{model_factory} rather than a reused mutable model when resampling. Functional CCA uses spline-basis coefficients and must be validated at the experimental-unit level.}
\examples{
## Probabilistic CCA via EM (internal engine, no optional packages needed)
d <- cca_agri_data(n = 80, seed = 1)
fit <- cca_probabilistic(d$X, d$Y, ncomp = 2)
fit
fit$converged

fit <- cca_probabilistic(as.matrix(iris[, 1:2]), as.matrix(iris[, 3:4]), ncomp = 2)
fit

fit <- cca_probabilistic(as.matrix(mtcars[, 1:4]), as.matrix(mtcars[, 5:7]), ncomp = 2)
fit

## Bayesian probabilistic CCA (requires rstan; tiny subsets keep the demo fast).
## Compiling a Stan model needs a working C++ toolchain and a StanHeaders version
## matching rstan; when the local installation cannot build one, the example stays
## quiet instead of failing the check.
if (requireNamespace("rstan", quietly = TRUE)) {
  d <- cca_agri_data(n = 40, seed = 1)
  bf <- try(cca_bayes(as.matrix(d$X[1:40, c("soil_N", "soil_P")]),
                      as.matrix(d$Y[1:40, c("photosynthesis", "biomass")]),
                      ncomp = 2, iter = 120, chains = 1, refresh = 0, seed = 1),
                silent = TRUE)
  if (!inherits(bf, "try-error")) {
    print(bf)
    print(bf$stanfit, pars = c("Wx", "Wy", "sigma_x", "sigma_y"))
    print(cca_bayes(as.matrix(iris[1:40, 1:2]), as.matrix(iris[1:40, 3:4]),
                    ncomp = 2, iter = 120, chains = 1, refresh = 0, seed = 1))
    print(cca_bayes(as.matrix(mtcars[, c("mpg", "cyl")]),
                    as.matrix(mtcars[, c("disp", "hp")]),
                    ncomp = 2, iter = 120, chains = 1, refresh = 0, seed = 1))
  }
}

## Deep CCA through pure-R callbacks (no Python required)
deep_fit <- function(X, Y, ncomp, ...) {
  Sxx <- crossprod(X); Syy <- crossprod(Y); Sxy <- crossprod(X, Y)
  isqrt <- function(M) {
    e <- eigen((M + t(M)) / 2, symmetric = TRUE)
    e$vectors \%*\% diag(1 / sqrt(pmax(e$values, 1e-8)), length(e$values)) \%*\% t(e$vectors)
  }
  Wx <- isqrt(Sxx); Wy <- isqrt(Syy)
  s <- svd(Wx \%*\% Sxy \%*\% Wy, nu = ncomp, nv = ncomp)
  list(scores_x = X \%*\% Wx \%*\% s$u,
       scores_y = Y \%*\% Wy \%*\% s$v,
       model = list(Wx = Wx \%*\% s$u, Wy = Wy \%*\% s$v))
}
deep_transform <- function(model, X, Y, ...) {
  list(x = X \%*\% model$Wx, y = Y \%*\% model$Wy)
}

d <- cca_agri_data(n = 60, seed = 1)
fit <- cca_deep(d$X, d$Y, fit_fun = deep_fit, transform_fun = deep_transform, ncomp = 2)
fit
predict(fit, newdata_x = d$X[1:5, ], newdata_y = d$Y[1:5, ])

fit <- cca_deep(as.matrix(iris[, 1:2]), as.matrix(iris[, 3:4]),
                fit_fun = deep_fit, ncomp = 2)
fit

fit <- cca_deep(as.matrix(mtcars[, 1:4]), as.matrix(mtcars[, 5:7]),
                fit_fun = deep_fit, ncomp = 2)
fit

## CCA-Zoo bridge (requires reticulate and a configured 'cca-zoo' environment)
pyok <- requireNamespace("reticulate", quietly = TRUE) &&
        reticulate::py_available(initialize = FALSE) &&
        reticulate::py_module_available("cca_zoo")
if (pyok) {
  d <- cca_agri_data(n = 60, seed = 1)
  zfit <- cca_zoo(d$X, d$Y, model = "CCA", module = "linear", ncomp = 2)
  print(zfit)
  print(predict(zfit, newdata_x = d$X[1:5, ], newdata_y = d$Y[1:5, ]))

  zfit <- cca_zoo(as.matrix(iris[, 1:2]), as.matrix(iris[, 3:4]),
                  model = "CCA", module = "linear", ncomp = 2)
  print(zfit)

  zfit <- cca_zoo(as.matrix(mtcars[, 1:4]), as.matrix(mtcars[, 5:7]),
                  model = "CCA", module = "linear", ncomp = 2)
  print(zfit)
}

## Semiparametric sparse CCA for mixed data types (requires mixedCCA)
if (requireNamespace("mixedCCA", quietly = TRUE)) {
  b <- cca_simulate(n = 60, scenario = "binary", seed = 1)
  fit <- cca_mixed(b$X, b$Y, type_x = "binary", type_y = "binary")
  fit

  z <- cca_simulate(n = 60, scenario = "zero_inflated", seed = 1)
  fit <- cca_mixed(z$X, z$Y, type_x = "continuous", type_y = "continuous")
  fit

  set.seed(1); n <- 80; lat <- rnorm(n)
  X <- cbind(x1 = lat + rnorm(n), x2 = -0.5 * lat + rnorm(n), x3 = rnorm(n))
  Y <- cbind(y1 = rbinom(n, 1, plogis(lat)),
             y2 = rbinom(n, 1, plogis(-lat)),
             y3 = rbinom(n, 1, plogis(0.5 * lat)))
  fit <- cca_mixed(X, Y, type_x = "continuous", type_y = "binary")
  fit
}

## Functional CCA on spline-basis coefficients
f <- cca_simulate(n = 40, scenario = "functional", rho = 0.7, seed = 1)
fit <- cca_functional(f$X, f$Y, time_x = f$time_x, time_y = f$time_y, ncomp = 2)
fit

set.seed(1); n <- 40; tt <- seq(0, 1, length.out = 30)
a <- rnorm(n); b <- 0.7 * a + sqrt(1 - 0.7^2) * rnorm(n)
X <- outer(a, sin(2 * pi * tt)) + matrix(rnorm(n * length(tt), sd = 0.3), n)
Y <- outer(b, cos(2 * pi * tt)) + matrix(rnorm(n * length(tt), sd = 0.3), n)
fit <- cca_functional(X, Y, time_x = tt, time_y = tt, ncomp = 2)
fit

set.seed(2); n <- 40; tt <- seq(0, 1, length.out = 24)
u <- rnorm(n); v <- 0.6 * u + sqrt(1 - 0.6^2) * rnorm(n)
X <- outer(u, tt^2) + outer(u, tt) + matrix(rnorm(n * length(tt), sd = 0.3), n)
Y <- outer(v, tt^2) - outer(v, tt) + matrix(rnorm(n * length(tt), sd = 0.3), n)
fit <- cca_functional(X, Y, time_x = tt, time_y = tt, method = "ridge", ncomp = 2)
fit

## Dynamic lag-aligned CCA
X <- as.matrix(longley[, c("GNP", "Unemployed", "Armed.Forces")])
Y <- as.matrix(longley[, c("Employed", "Population")])
fit <- cca_dynamic(X, Y, lag = 1, method = "ridge", ncomp = 2)
fit

aq <- airquality[complete.cases(airquality), ]
fit <- cca_dynamic(as.matrix(aq[, c("Ozone", "Solar.R")]),
                   as.matrix(aq[, c("Wind", "Temp")]),
                   lag = 0, method = "ridge", ncomp = 2)
fit

set.seed(1); n <- 60
ar1 <- function(n, phi, sd) {
  x <- numeric(n); x[1] <- rnorm(1)
  for (i in 2:n) x[i] <- phi * x[i - 1] + rnorm(1, sd = sd)
  x
}
z <- rnorm(n)
X <- cbind(x1 = ar1(n, 0.6, 0.5) + z, x2 = ar1(n, 0.4, 0.5) - z, x3 = ar1(n, 0.2, 0.5))
Y <- cbind(y1 = ar1(n, 0.6, 0.5) + 0.8 * z, y2 = ar1(n, 0.4, 0.5))
fit <- cca_dynamic(X, Y, lag = 1, method = "ridge", ncomp = 2)
fit

## Repeated-measures long-to-wide CCA
d <- cca_agri_data(n = 40, seed = 1)
ids <- seq_len(nrow(d$X))
long <- do.call(rbind, lapply(1:3, function(t) {
  data.frame(id = ids, time = t,
             soil_N = d$X$soil_N + 0.5 * t * rnorm(length(ids)),
             soil_P = d$X$soil_P + 0.3 * t * rnorm(length(ids)),
             photosynthesis = d$Y$photosynthesis + 0.6 * t * rnorm(length(ids)),
             biomass = d$Y$biomass + 0.4 * t * rnorm(length(ids)))
}))
fit <- cca_repeated(long, id = "id", time = "time",
                    x_vars = c("soil_N", "soil_P"),
                    y_vars = c("photosynthesis", "biomass"),
                    method = "ridge", ncomp = 2)
fit

set.seed(1); n <- 50; Tn <- 4
latent <- rnorm(n)
long <- do.call(rbind, lapply(seq_len(Tn), function(t) {
  data.frame(id = seq_len(n), time = t,
             x1 = latent + rnorm(n) + 0.2 * t,
             x2 = -0.5 * latent + rnorm(n),
             y1 = 0.8 * latent + rnorm(n) + 0.1 * t,
             y2 = 0.6 * latent + rnorm(n))
}))
fit <- cca_repeated(long, "id", "time", c("x1", "x2"), c("y1", "y2"),
                    method = "ridge", ncomp = 2)
fit

ir <- iris[1:90, ]
long <- do.call(rbind, lapply(1:3, function(t) {
  data.frame(id = seq_len(nrow(ir)), time = t,
             Sepal.Length = ir$Sepal.Length + 0.1 * t * rnorm(nrow(ir)),
             Sepal.Width = ir$Sepal.Width + 0.1 * t * rnorm(nrow(ir)),
             Petal.Length = ir$Petal.Length + 0.1 * t * rnorm(nrow(ir)),
             Petal.Width = ir$Petal.Width + 0.1 * t * rnorm(nrow(ir)))
}))
fit <- cca_repeated(long, "id", "time",
                    c("Sepal.Length", "Sepal.Width"),
                    c("Petal.Length", "Petal.Width"),
                    method = "ridge", ncomp = 2)
fit

## Multiway CCA via explicit tensor matricization
set.seed(1)
X <- array(rnorm(20 * 4 * 3), c(20, 4, 3))
Y <- array(rnorm(20 * 3 * 2), c(20, 3, 2))
fit <- cca_multiway(X, Y, ncomp = 2)
fit

set.seed(2); n <- 30
lat <- rnorm(n)
X <- array(outer(lat, 1:8) + rnorm(n * 8, sd = 0.5), c(n, 4, 2))
Y <- array(outer(0.8 * lat, 1:6) + rnorm(n * 6, sd = 0.5), c(n, 3, 2))
fit <- cca_multiway(X, Y, ncomp = 2)
fit

set.seed(3); n <- 40
u <- rnorm(n); v <- 0.7 * u + sqrt(1 - 0.7^2) * rnorm(n)
X <- array(outer(u, sin(seq(0, pi, length.out = 6))) + rnorm(n * 6, sd = 0.4), c(n, 3, 2))
Y <- array(outer(v, cos(seq(0, pi, length.out = 4))) + rnorm(n * 4, sd = 0.4), c(n, 2, 2))
fit <- cca_multiway(X, Y, ncomp = 2)
fit
}
