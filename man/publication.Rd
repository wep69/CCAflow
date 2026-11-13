\name{cca_plot}
\alias{cca_plot}\alias{theme_cca_publication}\alias{cca_palette}\alias{cca_export}\alias{cca_table}\alias{cca_report}
\title{Publication-ready plots, tables and reports}
\description{Create editable plots, standardized tables, vector/raster exports and reproducible canonical-correlation reports.}
\usage{
cca_plot(object, type = "scree", ..., interactive = FALSE)
theme_cca_publication(base_size = 10, base_family = "", grid = FALSE)
cca_palette(n = 8, grayscale = FALSE)
cca_export(plot, file, preset = c("single", "double", "custom"),
  width = NULL, height = NULL, units = c("in", "cm", "mm"),
  dpi = 600, bg = "white", ...)
cca_table(object, type = "correlations",
  format = c("data.frame", "markdown", "latex", "html", "csv",
  "word"), file = NULL)
cca_report(object, file = "cca_report.md",
  format = c("markdown", "html", "word"), bootstrap = NULL,
  validation = NULL, permutation = NULL,
  title = "Canonical Correlation Analysis Report",
  include_figures = TRUE, quiet = TRUE, extra_figures = FALSE)
}
\arguments{
\item{object}{CCAflow object.}
\item{type}{Plot or table type. For \code{cca_plot()} on a fitted object the accepted keys are \code{"scree"}, \code{"scores"}, \code{"circle"}, \code{"dual_circle"}, \code{"loadings"}, \code{"cross_loadings"}, \code{"coefficients"}, \code{"biplot"}, \code{"redundancy"}, \code{"heatmap"}, \code{"network"}, \code{"sparsity"}, \code{"influence"}, \code{"communality"}, \code{"suppressor"}, \code{"bipolar_lollipop"}, \code{"weight_vs_loading"} (alias \code{"tornado"}), \code{"synoptic_heatmap"}, \code{"heliograph"}, \code{"bipolar_balance"} and \code{"heliograph_panel"}. For \code{cca_table()} on a fitted object, \code{"participation"} returns the unified variable-level diagnostic table. Figures of inference and validation are not drawn from a fit: build the companion object first and plot it (\code{plot(cca_bootstrap(fit))}, \code{plot(cca_validate(...))}, \code{plot(cca_permutation(...))}, \code{plot(cca_compare(...))}, \code{plot(cca_tune(...))}, \code{plot(cca_stability(...))}). \code{cca_plot()} also accepts those objects and \code{cca_weight_diagnostics} directly and maps \code{type} to the same route.}
\item{interactive}{Convert supported ggplot outputs to plotly (not supported for \code{"heliograph"} and \code{"heliograph_panel"}).}
\item{base_size,base_family,grid}{Publication theme settings.}
\item{n,grayscale}{Palette controls.}
\item{plot,file,preset,width,height,units,dpi,bg}{Figure output controls.}
\item{format}{Requested output format.}
\item{bootstrap,validation,permutation}{Optional companion objects appended to the report.}
\item{title}{Report title.}
\item{include_figures,quiet,extra_figures}{Report figure inclusion, rendering verbosity and optional inclusion of suppressor/bipolar diagnostic figures.}
\item{...}{Additional plot, table or report controls, such as \code{components}, \code{block} (\code{"x"}, \code{"y"} or \code{"both"}), \code{style} (\code{"dot"} or \code{"lollipop"}), \code{polarity_color}, \code{thresholds}, \code{vif_cap}, \code{squish}, \code{variable_labels}, \code{order_by_component}, \code{sparsity_mode} (\code{"binary"}, \code{"signed"} or \code{"magnitude"}) and \code{highlight_zero}.}
}
\value{A ggplot/plotly object, palette, output path, data frame/formatted table, or report path.}
\details{Static figures are editable ggplot objects whenever possible. PDF/SVG are preferred for vector output and TIFF/PNG can be exported at journal-ready resolution.}
\examples{
## Reusable fits on disparate datasets (fast, reproducible)
ag <- cca_agri_data(n = 120, seed = 1)
fit_ag <- cca_classic(ag$X, ag$Y, ncomp = 2)
fit_sp <- cca_sparse(ag$X, ag$Y, ncomp = 2, penalty_x = 0.1, penalty_y = 0.1,
                     engine = "internal")

hd <- cca_simulate(n = 60, scenario = "highdim", seed = 2)
fit_ridge <- cca_ridge(hd$X, hd$Y, lambda_x = 0.2, lambda_y = 0.2, ncomp = 3)

nm <- cca_simulate(n = 80, scenario = "nonmonotone", seed = 3)
fit_kernel <- cca_kernel(nm$X, nm$Y, kernel = "rbf", ncomp = 2)

fit_mt <- cca_classic(mtcars[, 1:3], mtcars[, 4:6], ncomp = 2)

## ---- cca_plot: core types and new bipolar/suppressor diagnostics ----
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_scree <- cca_plot(fit_ag, "scree")
  p_scores <- cca_plot(fit_ag, "scores")
  p_circle <- cca_plot(fit_ridge, "circle")
  p_heatmap <- cca_plot(fit_kernel, "heatmap")
  p_biplot <- cca_plot(fit_ag, "biplot")
  p_redundancy <- cca_plot(fit_ridge, "redundancy")
  p_network <- cca_plot(fit_ag, "network", threshold = 0.3)

  ## Extended coefficients & sparsity styles (>= 2 self-contained examples each)
  p_coef_1 <- cca_plot(fit_ag, "coefficients", components = 1:2,
                       block = "both", style = "lollipop", polarity_color = TRUE)
  p_coef_2 <- cca_plot(fit_mt, "coefficients", components = 1:2,
                       block = "both", style = "dot", polarity_color = TRUE)

  p_sp_1 <- cca_plot(fit_sp, "sparsity", block = "both", sparsity_mode = "signed")
  p_sp_2 <- cca_plot(fit_sp, "sparsity", block = "both", sparsity_mode = "magnitude")

  ## Suppressor & VIF diagnostics
  p_supp_1 <- cca_plot(fit_ag, "suppressor", components = 1:2)
  p_supp_2 <- cca_plot(fit_ridge, "suppressor", components = 1:2, vif_cap = 20)

  ## Bipolar lollipop participation
  p_lol_1 <- cca_plot(fit_ag, "bipolar_lollipop", components = 1:2)
  p_lol_2 <- cca_plot(fit_ridge, "bipolar_lollipop", components = 1:2,
                      metric = "cross_loadings", order_by_component = 1)

  ## Tornado / weight vs. loading duality
  p_tor_1 <- cca_plot(fit_ag, "weight_vs_loading", components = 1:2)
  p_tor_2 <- cca_plot(fit_ridge, "tornado", components = 1:2,
                      order_by_component = 1, top = 10)

  ## Multi-dimensional synoptic heatmap
  p_syn_1 <- cca_plot(fit_ag, "synoptic_heatmap", components = 1:2)
  p_syn_2 <- cca_plot(fit_ridge, "synoptic_heatmap", components = 1:3, squish = 0.9)

  ## Polar heliograph
  p_hel_1 <- cca_plot(fit_ag, "heliograph", component = 1)
  p_hel_2 <- cca_plot(fit_ridge, "heliograph", component = 2, top = 8)

  ## Bipolar variance balance
  p_bal_1 <- cca_plot(fit_ag, "bipolar_balance", components = 1:2)
  p_bal_2 <- cca_plot(fit_ridge, "bipolar_balance", components = 1:3)

  ## Bootstrap bipolar weights with zero-crossing highlight
  boot_ex1 <- cca_bootstrap(fit_ag, B = 15, seed = 10)
  boot_ex2 <- cca_bootstrap(fit_mt, B = 15, seed = 11)
  p_bt_1 <- cca_plot(boot_ex1, quantity = "weight", block = "both", highlight_zero = TRUE)
  p_bt_2 <- cca_plot(boot_ex2, quantity = "weight", block = "both",
                     highlight_zero = TRUE, grayscale = TRUE)

  ## Direct weight-diagnostics plotting via cca_plot()
  p_wd_1 <- cca_plot(cca_weight_diagnostics(fit_ag))
  p_wd_2 <- cca_plot(cca_weight_diagnostics(fit_ridge, components = 1:2))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_hpan_1 <- cca_plot(fit_ag, "heliograph_panel", components = 1:2)
    p_hpan_2 <- cca_plot(fit_ridge, "heliograph_panel", components = 1:2, top = 8)
  }
  class(p_scree)
}

## Interactive plotly conversion (optional package)
if (requireNamespace("plotly", quietly = TRUE)) {
  p_plotly <- cca_plot(fit_ag, "scores", interactive = TRUE)
  class(p_plotly)
}

## ---- theme_cca_publication: size and grid variants ----
if (requireNamespace("ggplot2", quietly = TRUE)) {
  base_plot <- ggplot2::ggplot(ag$X, ggplot2::aes(soil_N, soil_P)) +
    ggplot2::geom_point()

  p_default <- base_plot + theme_cca_publication()
  p_large   <- base_plot + theme_cca_publication(base_size = 12)
  p_grid    <- base_plot + theme_cca_publication(base_size = 8, grid = TRUE)
  class(p_grid)
}

## ---- cca_palette: color, grayscale and short palettes ----
pal8 <- cca_palette(8)
length(pal8)

pal_gray <- cca_palette(5, grayscale = TRUE)
pal_gray

pal3 <- cca_palette(3)
pal3

## ---- cca_export: vector and raster journal-ready exports ----
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p1 <- cca_plot(fit_ag, "scores")

  f_pdf <- tempfile(fileext = ".pdf")
  cca_export(p1, f_pdf, preset = "single")

  f_tiff <- tempfile(fileext = ".tiff")
  cca_export(p1, f_tiff, preset = "double", dpi = 600)

  f_svg <- tempfile(fileext = ".svg")
  cca_export(p1, f_svg, preset = "single")

  file.exists(c(f_pdf, f_tiff, f_svg))
}

## ---- cca_table: types and formats across three datasets ----
head(cca_table(fit_ag, "correlations"))

head(cca_table(fit_ag, "coefficients"))

head(cca_table(fit_ag, "loadings"))

cat(cca_table(fit_ag, "redundancy", format = "markdown"))

head(cca_table(fit_ridge, "cross_loadings"))

cat(cca_table(fit_ridge, "correlations", format = "latex"))

f_csv <- tempfile(fileext = ".csv")
cca_table(fit_ridge, "diagnostics", format = "csv", file = f_csv)
head(read.csv(f_csv))

head(cca_table(fit_kernel, "communality"))

cat(cca_table(fit_kernel, "correlations", format = "markdown"))

cat(cca_table(fit_kernel, "correlations", format = "html"))

head(cca_table(fit_ag, "participation"))

cat(cca_table(fit_mt, "participation", format = "markdown"))

## Word export requires flextable and officer (optional packages)
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  ft <- cca_table(fit_ag, "correlations", format = "word")
  class(ft)
}

## ---- cca_report: reproducible markdown reports ----
boot_ag <- cca_bootstrap(fit_ag, B = 19, seed = 1)
f1 <- tempfile(fileext = ".md")
cca_report(fit_ag, file = f1, format = "markdown", bootstrap = boot_ag,
           include_figures = FALSE)
file.exists(f1)

val_ridge <- cca_validate(hd$X, hd$Y, method = "ridge",
                          params = list(lambda_x = 0.2, lambda_y = 0.2, ncomp = 3),
                          folds = 3, repeats = 2, seed = 1)
perm_ridge <- suppressWarnings(cca_permutation(fit_ridge, B = 19, seed = 2))
f2 <- tempfile(fileext = ".md")
cca_report(fit_ridge, file = f2, format = "markdown",
           validation = val_ridge, permutation = perm_ridge,
           include_figures = FALSE)
file.exists(f2)

f3 <- tempfile(fileext = ".md")
cca_report(fit_kernel, file = f3, format = "markdown", include_figures = FALSE)
file.exists(f3)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  f4 <- tempfile(fileext = ".md")
  cca_report(fit_ag, file = f4, format = "markdown",
             include_figures = TRUE, extra_figures = TRUE)
  f5 <- tempfile(fileext = ".md")
  cca_report(fit_mt, file = f5, format = "markdown",
             include_figures = TRUE, extra_figures = TRUE)
  file.exists(c(f4, f5))
}

## HTML and Word rendering require rmarkdown with a working pandoc
if (requireNamespace("rmarkdown", quietly = TRUE) &&
    rmarkdown::pandoc_available()) {
  f_html <- tempfile(fileext = ".html")
  cca_report(fit_ag, file = f_html, format = "html", include_figures = FALSE)
  f_word <- tempfile(fileext = ".docx")
  cca_report(fit_ag, file = f_word, format = "word", include_figures = FALSE)
}
}
