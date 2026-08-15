\name{cca_teach}
\alias{cca_teach}\alias{cca_simulate}\alias{cca_agri_data}\alias{cca_tour}\alias{cca_grand_tour}
\alias{cca_guided_tour}\alias{cca_basis_tour}\alias{cca_compare_tour}\alias{cca_outlier_tour}\alias{cca_tour_keyframes}
\alias{cca_tour_render}\alias{cca_explorer}
\title{Canonical Correlation Learning Lab and dynamic projection tours}
\description{Generate controlled simulations, agronomy teaching data, canonical projection paths, guided paths, publication keyframes and layered tour renderings.}
\usage{
cca_simulate(n = 150, p = 6, q = 5,
  scenario = c("linear", "weak", "null", "collinear", "highdim",
    "sparse", "outliers", "heavy_tail", "monotone", "nonmonotone",
    "binary", "zero_inflated", "missing", "functional",
    "three_blocks"),
  rho = 0.8, seed = 1)
cca_agri_data(n = 120, seed = 1)
cca_teach(X = NULL, Y = NULL,
  level = c("intro", "intermediate", "advanced"), seed = 1)
cca_tour(object, block = c("x", "y"), frames = 40,
  start = c("pca", "identity"), renderer = c("none", "spinifex"),
  animate = FALSE)
cca_grand_tour(object, block = c("x", "y"), ...)
cca_guided_tour(object, block = c("x", "y"), frames = 40,
  start = "pca")
cca_basis_tour(data, start, target, frames = 40)
cca_compare_tour(object1, object2, block = c("x", "y"), frames = 40)
cca_outlier_tour(X, Y, block = c("x", "y"), frames = 40,
  robust_engine = "winsor")
cca_tour_keyframes(tour, n = 6, strategy = c("even", "objective"),
  plots = FALSE, group = NULL)
cca_tour_render(tour, engine = c("spinifex", "static"), frame = NULL,
  animate = c("none", "plotly", "gganimate"), group = NULL,
  point_size = 1.7, alpha = 0.75, ...)
cca_explorer(object, component = 1, top = 15, interactive = FALSE)
}
\arguments{
\item{n,p,q,rho,seed,scenario}{Simulation controls.}
\item{X,Y}{Optional paired blocks.}
\item{level}{Teaching depth.}
\item{object,object1,object2}{Fitted CCA objects.}
\item{block}{Block to tour.}
\item{frames,renderer,animate}{Tour construction and rendering controls.}
\item{data,start,target}{Raw data and start/target bases for \code{cca_basis_tour()}.}
\item{robust_engine}{Resistant engine used by \code{cca_outlier_tour()}.}
\item{tour}{A \code{cca_tour} object.}
\item{strategy,plots,group}{Keyframe and grouping controls.}
\item{engine,frame}{Rendering engine and selected static frame.}
\item{point_size,alpha}{Point size and transparency for tour renderings.}
\item{component,top,interactive}{Explorer controls.}
\item{...}{Additional rendering arguments.}
}
\value{Simulation data, a teaching-lab object, a tour object, a list of keyframes, a plot/widget, or a four-panel explorer.}
\details{The canonical tour interpolates from an initial projection to a canonical basis. \pkg{spinifex} can render the retained basis path as a layered tour and animate it via plotly or gganimate. Static keyframes remain reproducible and suitable for manuscript panels.}
\examples{
## Reusable fits and discrepant datasets (fast, reproducible)
ag <- cca_agri_data(n = 120, seed = 1)
fit_ag <- cca_classic(ag$X, ag$Y)
iris_X <- iris[c("Sepal.Length", "Sepal.Width")]
iris_Y <- iris[c("Petal.Length", "Petal.Width")]
fit_iris <- cca_classic(iris_X, iris_Y)
lin <- cca_simulate(n = 80, scenario = "linear", rho = 0.8, seed = 1)
fit_lin <- cca_classic(lin$X, lin$Y)
mtcars_X <- mtcars[c("mpg", "disp", "hp", "wt")]
mtcars_Y <- mtcars[c("drat", "qsec", "gear")]
fit_mtcars <- cca_classic(mtcars_X, mtcars_Y)

## ---- cca_teach: three teaching depths ----
lab_intro <- cca_teach(level = "intro", seed = 1)
print(lab_intro$fit$correlations)
print(lab_intro$narrative[1])

lab_inter <- cca_teach(level = "intermediate", seed = 2)
print(lab_inter$fit$correlations)

lab_adv <- cca_teach(level = "advanced", seed = 3)
print(lab_adv$fit$correlations)

## ---- cca_simulate: controlled scenarios (linear, high-dimensional, multiblock) ----
sim_linear <- cca_simulate(n = 80, scenario = "linear", rho = 0.8, seed = 1)
print(sim_linear$truth)
dim(sim_linear$X)

sim_highdim <- cca_simulate(n = 60, scenario = "highdim", seed = 2)
dim(sim_highdim$X)
dim(sim_highdim$Y)

sim_blocks <- cca_simulate(n = 60, scenario = "three_blocks", seed = 3)
names(sim_blocks$blocks)
lapply(sim_blocks$blocks, dim)

## Outlier and nonlinear scenarios feed the tour functions below
sim_outliers <- cca_simulate(n = 80, scenario = "outliers", seed = 4)
sim_nonmono  <- cca_simulate(n = 80, scenario = "nonmonotone", seed = 5)

## ---- cca_agri_data: three sample sizes with distinct seeds ----
ag60 <- cca_agri_data(n = 60, seed = 10)
nrow(ag60$X)
names(ag60$X)
names(ag60$Y)
levels(ag60$treatment)

ag120 <- cca_agri_data(n = 120, seed = 20)
nrow(ag120$X)

ag200 <- cca_agri_data(n = 200, seed = 30)
nrow(ag200$X)

## ---- cca_tour: canonical projection paths (renderer = "none") ----
tour_agri <- cca_tour(fit_ag, block = "x", frames = 12, renderer = "none")
print(tour_agri)

tour_iris <- cca_tour(fit_iris, block = "y", frames = 12, renderer = "none")
length(tour_iris$bases)

tour_lin <- cca_tour(fit_lin, block = "x", frames = 12, renderer = "none")
length(tour_lin$bases)

## ---- cca_grand_tour: tourr grand-tour wrapper (optional package) ----
if (requireNamespace("tourr", quietly = TRUE)) {
  gt_agri <- cca_grand_tour(fit_ag, block = "x")
  gt_iris <- cca_grand_tour(fit_iris, block = "x")
  gt_mtcars <- cca_grand_tour(fit_mtcars, block = "y")
  class(gt_agri)
  class(gt_iris)
  class(gt_mtcars)
}

## ---- cca_guided_tour: association-guided paths ----
guided_agri <- cca_guided_tour(fit_ag, block = "x", frames = 12)
length(guided_agri$objective)

guided_iris <- cca_guided_tour(fit_iris, block = "y", frames = 12)
length(guided_iris$objective)

guided_lin <- cca_guided_tour(fit_lin, block = "x", frames = 12)
round(head(guided_lin$objective), 3)

## ---- cca_basis_tour: interpolate between two projection bases ----
tr_agri <- cca_tour(fit_ag, block = "x", frames = 12, renderer = "none")
bt_agri <- cca_basis_tour(tr_agri$data, tr_agri$bases[[1]], tr_agri$bases[[12]],
                          frames = 12)
length(bt_agri$bases)

cmp_out <- cca_robust_compare(sim_outliers$X, sim_outliers$Y)
tr_classic <- cca_tour(cmp_out$classical, block = "x", frames = 12, renderer = "none")
tr_robust  <- cca_tour(cmp_out$resistant, block = "x", frames = 12, renderer = "none")
bt_classic_robust <- cca_basis_tour(cmp_out$classical$Xs,
                                    tr_classic$bases[[12]], tr_robust$bases[[12]],
                                    frames = 12)
length(bt_classic_robust$bases)

tr_iris <- cca_tour(fit_iris, block = "x", frames = 12, renderer = "none")
bt_iris <- cca_basis_tour(tr_iris$data, tr_iris$bases[[1]], tr_iris$bases[[12]],
                          frames = 12)
length(bt_iris$bases)

## ---- cca_compare_tour: two fitted bases on a shared variable space ----
fit_ag_ridge <- cca_ridge(ag$X, ag$Y, lambda_x = 0.2, lambda_y = 0.2)
cmp_tour_agri <- cca_compare_tour(fit_ag, fit_ag_ridge, block = "x", frames = 12)
length(cmp_tour_agri$bases)

fit_lin_ridge <- cca_ridge(lin$X, lin$Y, lambda_x = 0.2, lambda_y = 0.2)
cmp_tour_lin <- cca_compare_tour(fit_lin, fit_lin_ridge, block = "x", frames = 12)
length(cmp_tour_lin$bases)

fit_mtcars_ridge <- cca_ridge(mtcars_X, mtcars_Y, lambda_x = 0.2, lambda_y = 0.2)
cmp_tour_mtcars <- cca_compare_tour(fit_mtcars, fit_mtcars_ridge,
                                    block = "x", frames = 12)
length(cmp_tour_mtcars$bases)

## ---- cca_outlier_tour: classical-versus-resistant paths ----
ot_outliers <- cca_outlier_tour(sim_outliers$X, sim_outliers$Y,
                                block = "x", frames = 12)
print(ot_outliers$type)

ht <- cca_simulate(n = 80, scenario = "heavy_tail", seed = 6)
ot_heavy <- cca_outlier_tour(ht$X, ht$Y, block = "y", frames = 12)
print(ot_heavy$type)

ot_agri <- cca_outlier_tour(ag$X, ag$Y, block = "x", frames = 12)
print(ot_agri$type)

## ---- cca_tour_keyframes: static keyframes for manuscript panels ----
kf_agri <- cca_tour_keyframes(tour_agri, n = 3, strategy = "even")
print(kf_agri$index)

kf_guided <- cca_tour_keyframes(guided_lin, n = 3, strategy = "objective")
print(kf_guided$index)

kf_iris <- cca_tour_keyframes(tour_iris, n = 5, strategy = "even", plots = TRUE)
length(kf_iris$plots)

## ---- cca_tour_render: static keyframe and layered spinifex rendering ----
p_frame <- cca_tour_render(tour_agri, engine = "static", frame = 6)
class(p_frame)

p_group <- cca_tour_render(tour_agri, engine = "static", frame = 12,
                           group = ag$treatment)
class(p_group)

## spinifex renders the retained path; angle = 0 keeps the stored frames
if (requireNamespace("spinifex", quietly = TRUE)) {
  g_spinifex <- cca_tour_render(tour_lin, engine = "spinifex", angle = 0)
  class(g_spinifex)
}

## ---- cca_explorer: four-panel association explorer ----
ex_agri <- cca_explorer(fit_ag, component = 1, interactive = FALSE)
names(ex_agri)

ex_iris <- cca_explorer(fit_iris, component = 1, interactive = FALSE)
names(ex_iris)

ex_mtcars <- cca_explorer(fit_mtcars, component = 1, interactive = FALSE)
names(ex_mtcars)
}
