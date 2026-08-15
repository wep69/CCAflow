# Publication-ready graphics and tables -------------------------------------

#' Accessible publication palette
#' @export
cca_palette <- function(n = 8, grayscale = FALSE) {
  if (grayscale) {
    pal <- grDevices::gray.colors(max(n, 2), start = 0.15, end = 0.80)
  } else {
    pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
    if (n > length(pal)) pal <- grDevices::colorRampPalette(pal)(n)
  }
  pal[seq_len(n)]
}

#' Publication theme for CCAflow plots
#' @export
theme_cca_publication <- function(base_size = 10, base_family = "", grid = FALSE) {
  .cca_require("ggplot2", "theme_cca_publication")
  th <- ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(hjust = 0),
      axis.title = ggplot2::element_text(face = "plain"),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(linewidth = 0.4),
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text = ggplot2::element_text(face = "bold")
    )
  if (!grid) th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  th
}

.cca_component_name <- function(k) paste0("CV", k)

.cca_plot_fit <- function(object, type, component = 1, component_y = 2,
                          block = c("x", "y"), top = 20, threshold = 0.3,
                          groups = NULL, labels = NULL, bootstrap = NULL,
                          permutation = NULL, grayscale = FALSE,
                          group_ellipse = FALSE, ellipse_level = 0.95, ...) {
  .cca_require("ggplot2", "cca_plot")
  block <- match.arg(block)
  pal <- cca_palette(8, grayscale)
  type <- match.arg(type, c("scree", "scores", "circle", "dual_circle", "loadings",
                            "cross_loadings", "coefficients", "biplot", "redundancy",
                            "heatmap", "network", "sparsity", "influence"))
  if (type == "scree") {
    d <- data.frame(component = seq_along(object$correlations), correlation = object$correlations)
    p <- ggplot2::ggplot(d, ggplot2::aes(component, correlation)) +
      ggplot2::geom_line(linewidth = 0.5) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_x_continuous(breaks = d$component) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = "Canonical dimension", y = "Canonical correlation", title = "Canonical correlation profile")
    if (!is.null(bootstrap) && inherits(bootstrap, "cca_bootstrap")) {
      ci <- bootstrap$correlation_ci
      p <- p + ggplot2::geom_errorbar(data = ci, ggplot2::aes(x = component, ymin = lower, ymax = upper), width = 0.12)
    }
    if (!is.null(permutation) && inherits(permutation, "cca_permutation")) {
      env <- apply(permutation$null, 2, stats::quantile, probs = 0.95, na.rm = TRUE)
      e <- data.frame(component = seq_along(env), null95 = env)
      p <- p + ggplot2::geom_line(data = e, ggplot2::aes(component, null95), linetype = 2)
    }
    return(p + theme_cca_publication())
  }
  if (type == "scores") {
    k <- min(component, ncol(object$scores_x), ncol(object$scores_y))
    d <- data.frame(U = object$scores_x[, k], V = object$scores_y[, k])
    if (is.null(groups)) groups <- object$metadata$groups
    if (!is.null(groups)) d$group <- as.factor(groups)
    if (!is.null(labels)) d$label <- labels
    aes <- if (!is.null(groups)) ggplot2::aes(U, V, shape = group, group = group) else ggplot2::aes(U, V)
    p <- ggplot2::ggplot(d, aes) + ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.5) +
      ggplot2::labs(x = paste0("X canonical score, CV", k), y = paste0("Y canonical score, CV", k),
                    title = sprintf("Paired canonical scores (r = %.3f)", object$correlations[k]))
    if (isTRUE(group_ellipse) && !is.null(groups)) {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes(linetype = group), type = "norm",
                                     level = ellipse_level, linewidth = 0.45, show.legend = TRUE)
    }
    if (!is.null(labels)) p <- p + ggplot2::geom_text(ggplot2::aes(label = label), check_overlap = TRUE, vjust = -0.6, size = 3)
    return(p + theme_cca_publication())
  }
  if (type %in% c("circle", "dual_circle")) {
    kx <- min(component, ncol(object$loadings_x)); ky <- min(component_y, ncol(object$loadings_x))
    make <- function(M, b) data.frame(variable = rownames(M), x = M[, kx], y = if (ncol(M) >= ky && ky != kx) M[, ky] else 0, block = b)
    d <- if (type == "dual_circle") rbind(make(object$loadings_x, "X"), make(object$loadings_y, "Y"))
         else make(if (block == "x") object$loadings_x else object$loadings_y, toupper(block))
    theta <- seq(0, 2*pi, length.out = 300); circle <- data.frame(x = cos(theta), y = sin(theta))
    p <- ggplot2::ggplot(d, ggplot2::aes(0, 0)) +
      ggplot2::geom_path(data = circle, ggplot2::aes(x, y), inherit.aes = FALSE, linewidth = 0.35) +
      ggplot2::geom_segment(ggplot2::aes(xend = x, yend = y, linetype = block), arrow = ggplot2::arrow(length = grid::unit(0.12, "cm")), linewidth = 0.45) +
      ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = variable), size = 3, check_overlap = TRUE, vjust = -0.4) +
      ggplot2::coord_equal(xlim = c(-1.08, 1.08), ylim = c(-1.08, 1.08)) +
      ggplot2::labs(x = .cca_component_name(kx), y = if (ky != kx) .cca_component_name(ky) else "0",
                    title = if (type == "dual_circle") "Dual canonical correlation circle" else paste("Canonical correlation circle:", toupper(block)))
    return(p + theme_cca_publication())
  }
  if (type %in% c("loadings", "cross_loadings", "coefficients")) {
    if (type == "loadings") M <- if (block == "x") object$loadings_x else object$loadings_y
    if (type == "cross_loadings") M <- if (block == "x") object$cross_loadings_x else object$cross_loadings_y
    if (type == "coefficients") {
      M <- if (block == "x") object$weights_x else object$weights_y
      if (is.null(M) || nrow(M) != if (block == "x") object$p else object$q) .cca_stop("Original-variable canonical coefficients are unavailable for this estimator.")
    }
    k <- min(component, ncol(M)); idx <- .cca_top(M[, k], top)
    d <- data.frame(variable = factor(rownames(M)[idx], levels = rownames(M)[idx][order(M[idx,k])]), value = M[idx, k])
    p <- ggplot2::ggplot(d, ggplot2::aes(value, variable)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) + ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = switch(type, loadings = "Structure coefficient", cross_loadings = "Cross-loading", coefficients = "Canonical coefficient"),
                    y = NULL, title = paste(gsub("_", " ", type), toupper(block), .cca_component_name(k)))
    return(p + theme_cca_publication())
  }
  if (type == "redundancy") {
    d <- rbind(data.frame(component = seq_along(object$redundancy_x), redundancy = object$redundancy_x, block = "X"),
               data.frame(component = seq_along(object$redundancy_y), redundancy = object$redundancy_y, block = "Y"))
    return(ggplot2::ggplot(d, ggplot2::aes(component, redundancy, linetype = block, shape = block)) +
      ggplot2::geom_line(linewidth = 0.5) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_x_continuous(breaks = unique(d$component)) +
      ggplot2::labs(x = "Canonical dimension", y = "Redundancy index", title = "Canonical redundancy") + theme_cca_publication())
  }
  if (type == "heatmap") {
    C <- stats::cor(object$Xs, object$Ys, use = "pairwise.complete.obs")
    d <- as.data.frame(as.table(C), stringsAsFactors = FALSE); names(d) <- c("X", "Y", "correlation")
    return(ggplot2::ggplot(d, ggplot2::aes(Y, X, fill = correlation)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = pal[2], mid = "white", high = pal[1], midpoint = 0, limits = c(-1,1)) +
      ggplot2::labs(x = "Y variables", y = "X variables", fill = "r", title = "Cross-block correlation heatmap") +
      theme_cca_publication() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  if (type == "network") {
    C <- stats::cor(object$Xs, object$Ys, use = "pairwise.complete.obs")
    ij <- which(abs(C) >= threshold, arr.ind = TRUE)
    if (!nrow(ij)) .cca_stop("No cross-block associations exceed the requested threshold.")
    dx <- data.frame(name = rownames(C), x = 0, y = seq(1, 0, length.out = nrow(C)), block = "X")
    dy <- data.frame(name = colnames(C), x = 1, y = seq(1, 0, length.out = ncol(C)), block = "Y")
    ed <- data.frame(x=0, y=dx$y[ij[,1]], xend=1, yend=dy$y[ij[,2]], value=C[ij])
    p <- ggplot2::ggplot() + ggplot2::geom_segment(data=ed, ggplot2::aes(x,y,xend=xend,yend=yend, linewidth=abs(value), linetype=value>0), alpha=.55) +
      ggplot2::geom_point(data=rbind(dx,dy), ggplot2::aes(x,y,shape=block), size=2) +
      ggplot2::geom_text(data=dx, ggplot2::aes(x,y,label=name), hjust=1.1, size=3) +
      ggplot2::geom_text(data=dy, ggplot2::aes(x,y,label=name), hjust=-.1, size=3) +
      ggplot2::scale_x_continuous(limits=c(-.35,1.35), breaks=c(0,1), labels=c("X","Y")) +
      ggplot2::guides(linewidth = "none", linetype = "none") +
      ggplot2::labs(x=NULL,y=NULL,title=paste0("Cross-block association network |r| >= ",threshold)) + theme_cca_publication()
    return(p)
  }
  if (type == "sparsity") {
    if (is.null(object$weights_x) || is.null(object$weights_y) || nrow(object$weights_x) != object$p || nrow(object$weights_y) != object$q) .cca_stop("Original-variable sparse weights are unavailable.")
    d <- rbind(.cca_long_matrix(object$weights_x, "X", "weight"), .cca_long_matrix(object$weights_y, "Y", "weight"))
    d$selected <- abs(d$value) > 1e-8
    return(ggplot2::ggplot(d, ggplot2::aes(component, variable, fill = selected)) + ggplot2::geom_tile() +
      ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_fill_manual(values = c("FALSE" = "white", "TRUE" = pal[1])) +
      ggplot2::labs(x = "Component", y = NULL, fill = "Selected", title = "Canonical sparsity map") + theme_cca_publication())
  }
  if (type == "influence") {
    d <- cca_diagnose(object)$data
    return(ggplot2::ggplot(d, ggplot2::aes(classical_distance, score_energy)) + ggplot2::geom_point(size=2, alpha=.75) +
      ggplot2::labs(x="Classical multivariate distance", y="Canonical score energy", title="CCA influence diagnostic") + theme_cca_publication())
  }
  if (type == "biplot") {
    kx <- min(component, ncol(object$scores_x)); ky <- min(component_y, ncol(object$scores_x))
    sc <- data.frame(x=object$scores_x[,kx], y=if(ky==kx) object$scores_y[,kx] else object$scores_x[,ky])
    L <- object$loadings_x[, c(kx, if(ky==kx) kx else ky), drop=FALSE]
    ld <- data.frame(variable=rownames(L), x=L[,1], y=if(ky==kx) 0 else L[,2])
    mult <- 0.7 * min(diff(range(sc$x))/2, diff(range(sc$y))/2)
    return(ggplot2::ggplot(sc, ggplot2::aes(x,y)) + ggplot2::geom_point(alpha=.65) +
      ggplot2::geom_segment(data=ld, ggplot2::aes(x=0,y=0,xend=x*mult,yend=y*mult), inherit.aes=FALSE,
                            arrow=ggplot2::arrow(length=grid::unit(.1,"cm"))) +
      ggplot2::geom_text(data=ld, ggplot2::aes(x=x*mult,y=y*mult,label=variable), inherit.aes=FALSE, size=3, check_overlap=TRUE) +
      ggplot2::labs(x=.cca_component_name(kx), y=.cca_component_name(ky), title="Canonical biplot") + theme_cca_publication())
  }
}

.cca_plot_other <- function(object, type, component = 1, block = "x", top = 20, ...) {
  .cca_require("ggplot2", "cca_plot")
  if (inherits(object, "cca_bootstrap")) {
    if (type %in% c("bootstrap", "scree")) {
      d <- object$correlation_ci
      return(ggplot2::ggplot(d, ggplot2::aes(component, estimate)) + ggplot2::geom_point(size=2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=lower,ymax=upper), width=.12) +
        ggplot2::scale_x_continuous(breaks=d$component) + ggplot2::coord_cartesian(ylim=c(0,1)) +
        ggplot2::labs(x="Canonical dimension",y="Canonical correlation",title="Bootstrap uncertainty") + theme_cca_publication())
    }
    quantity <- if (type %in% c("loading","cross_loading","weight")) type else "loading"
    d <- object$variable_draws
    d <- d[d$quantity == quantity & tolower(d$block) == tolower(block) & d$component == paste0("CV",component),,drop=FALSE]
    if (!nrow(d)) .cca_stop("Requested bootstrap quantity is unavailable.")
    sp <- split(d$value, d$variable)
    sm <- do.call(rbind, lapply(sp, function(z) c(median=stats::median(z,na.rm=TRUE), lower=stats::quantile(z,.025,na.rm=TRUE), upper=stats::quantile(z,.975,na.rm=TRUE))))
    dd <- data.frame(variable=rownames(sm), sm, row.names=NULL); idx <- .cca_top(dd$median, top); dd <- dd[idx,]
    dd$variable <- factor(dd$variable, levels=dd$variable[order(dd$median)])
    return(ggplot2::ggplot(dd, ggplot2::aes(median, variable)) + ggplot2::geom_vline(xintercept=0,linewidth=.3) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=lower,xmax=upper),height=.15) + ggplot2::geom_point(size=2) +
      ggplot2::labs(x=paste("Bootstrap",gsub("_"," ",quantity)),y=NULL,title=paste("Bootstrap",toupper(block),paste0("CV",component))) + theme_cca_publication())
  }
  if (inherits(object, "cca_validation")) {
    d <- object$results; if (!nrow(d)) .cca_stop("No validation results to plot.")
    return(ggplot2::ggplot(d, ggplot2::aes(factor(component), correlation)) + ggplot2::geom_hline(yintercept=0, linewidth=.3) +
      ggplot2::geom_boxplot(outlier.shape=NA) +
      ggplot2::geom_jitter(width=.12,height=0,alpha=.45,size=1.4) + ggplot2::coord_cartesian(ylim=c(-1,1)) +
      ggplot2::labs(x="Canonical dimension",y="Signed out-of-sample correlation",title="Canonical generalization") + theme_cca_publication())
  }
  if (inherits(object, "cca_permutation")) {
    k <- component; d <- data.frame(value=object$null[,k]); obs <- object$observed[k]
    return(ggplot2::ggplot(d, ggplot2::aes(value)) + ggplot2::geom_histogram(bins=35,boundary=0) +
      ggplot2::geom_vline(xintercept=obs,linetype=2,linewidth=.7) +
      ggplot2::labs(x=paste0("Null correlation, CV",k),y="Frequency",title=sprintf("Permutation null (observed = %.3f)",obs)) + theme_cca_publication())
  }
  if (inherits(object, "cca_robust_comparison")) {
    d <- object$comparison
    long <- rbind(data.frame(component=d$component, correlation=d$classical, estimator="Classical"),
                  data.frame(component=d$component, correlation=d$resistant, estimator="Resistant"))
    return(ggplot2::ggplot(long, ggplot2::aes(component, correlation, linetype=estimator, shape=estimator)) +
      ggplot2::geom_line(linewidth=.5) + ggplot2::geom_point(size=2) +
      ggplot2::scale_x_continuous(breaks=unique(long$component)) + ggplot2::coord_cartesian(ylim=c(0,1)) +
      ggplot2::labs(x="Canonical dimension", y="Canonical correlation", title="Classical versus resistant CCA") + theme_cca_publication())
  }
  if (inherits(object, "cca_comparison")) {
    d <- object$table[dplyr_safe_isfinite(object$table$training_rho1),,drop=FALSE]
    return(ggplot2::ggplot(d, ggplot2::aes(training_rho1, validation_rho1, label=method)) + ggplot2::geom_abline(linetype=2) +
      ggplot2::geom_point(size=2) + ggplot2::geom_text(vjust=-.6,check_overlap=TRUE,size=3) +
      ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(-1,1)) + ggplot2::labs(x="Training rho1",y="Signed validation r1",title="CCA method comparison") + theme_cca_publication())
  }
  if (inherits(object, "cca_tuning")) {
    d <- object$results
    pars <- setdiff(names(d), c("mean_validation","sd_validation")); pars <- pars[vapply(d[pars], is.numeric, logical(1))]
    if (length(pars) >= 2) {
      return(ggplot2::ggplot(d, ggplot2::aes(x = .data[[pars[1]]],
                                             y = .data[[pars[2]]],
                                             fill = .data[["mean_validation"]])) +
        ggplot2::geom_tile() + ggplot2::labs(x=pars[1],y=pars[2],fill="Validation",title="CCA tuning surface") + theme_cca_publication())
    }
    xvar <- if (length(pars)) pars[1] else names(d)[1]
    return(ggplot2::ggplot(d, ggplot2::aes(x = .data[[xvar]],
                                           y = .data[["mean_validation"]])) + ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(x=xvar,y="Mean validation association",title="CCA tuning") + theme_cca_publication())
  }
  if (inherits(object, "cca_stability")) {
    d <- object$table; d <- d[order(d$selection_frequency,decreasing=TRUE),]; d <- head(d, top)
    d$variable <- factor(d$variable, levels=rev(unique(d$variable)))
    return(ggplot2::ggplot(d, ggplot2::aes(selection_frequency,variable,shape=block)) + ggplot2::geom_point(size=2) +
      ggplot2::coord_cartesian(xlim=c(0,1)) + ggplot2::labs(x="Selection frequency",y=NULL,title="Sparse CCA stability") + theme_cca_publication())
  }
  .cca_stop("Unsupported object/type combination for cca_plot().")
}

# Helper avoids dplyr dependency.
dplyr_safe_isfinite <- function(x) is.finite(x)

#' Publication-ready plotting dispatcher
#' @export
cca_plot <- function(object, type = "scree", ..., interactive = FALSE) {
  p <- if (inherits(object, "cca_fit")) .cca_plot_fit(object, type = type, ...) else .cca_plot_other(object, type = type, ...)
  if (isTRUE(interactive)) {
    .cca_require("plotly", "cca_plot(interactive=TRUE)")
    return(plotly::ggplotly(p))
  }
  p
}

#' @export
plot.cca_fit <- function(x, type = "scree", ...) cca_plot(x, type = type, ...)
#' @export
plot.cca_bootstrap <- function(x, type = NULL, quantity = NULL, ...) { if (!is.null(quantity)) type <- quantity else if (is.null(type)) type <- "bootstrap"; cca_plot(x, type = type, ...) }
#' @export
plot.cca_validation <- function(x, type = "validation", ...) cca_plot(x, type = type, ...)
#' @export
plot.cca_permutation <- function(x, type = "permutation", ...) cca_plot(x, type = type, ...)
#' @export
plot.cca_comparison <- function(x, type = "comparison", ...) cca_plot(x, type = type, ...)
#' @export
plot.cca_robust_comparison <- function(x, type = "comparison", ...) cca_plot(x, type = type, ...)
#' @export
plot.cca_tuning <- function(x, type = "tuning", ...) cca_plot(x, type = type, ...)
#' @export
plot.cca_stability <- function(x, type = "stability", ...) cca_plot(x, type = type, ...)

#' Export publication-ready figures
#' @export
cca_export <- function(plot, file, preset = c("single", "double", "custom"),
                       width = NULL, height = NULL, units = c("in", "cm", "mm"), dpi = 600,
                       bg = "white", ...) {
  .cca_require("ggplot2", "cca_export")
  preset <- match.arg(preset); units <- match.arg(units)
  ext <- tolower(tools::file_ext(file))
  if (ext == "svg") .cca_require("svglite", "cca_export(file='*.svg')")
  # Journal-oriented reference sizes are defined in inches, then converted.
  wh_in <- switch(preset, single = c(3.35, 3.0), double = c(7.0, 4.5), custom = c(6.0, 4.0))
  convert <- switch(units, `in` = 1, cm = 2.54, mm = 25.4)
  if (is.null(width)) width <- wh_in[1] * convert
  if (is.null(height)) height <- wh_in[2] * convert
  dots <- list(...)
  if (ext %in% c("tif", "tiff") && !"compression" %in% names(dots)) dots$compression <- "lzw"
  do.call(ggplot2::ggsave, c(list(filename = file, plot = plot, width = width, height = height,
                                  units = units, dpi = dpi, bg = bg), dots))
  invisible(normalizePath(file, winslash = "/", mustWork = FALSE))
}

.cca_table_fit <- function(object, type) {
  type <- match.arg(type, c("correlations","coefficients","loadings","cross_loadings","redundancy","communality","diagnostics"))
  if (type == "correlations") {
    inf <- if (object$method == "classic") cca_infer(object) else NULL
    d <- data.frame(component=seq_along(object$correlations), correlation=object$correlations)
    if (!is.null(inf$table)) d <- merge(d, inf$table[,c("component","wilks_lambda","bartlett_chisq","df","p_value")], by="component", all.x=TRUE)
    return(d)
  }
  if (type == "coefficients") {
    if (is.null(object$weights_x) || nrow(object$weights_x) != object$p) .cca_stop("Original-variable coefficients unavailable for this estimator.")
    return(rbind(.cca_long_matrix(object$weights_x,"X","coefficient"), .cca_long_matrix(object$weights_y,"Y","coefficient")))
  }
  if (type == "loadings") return(rbind(.cca_long_matrix(object$loadings_x,"X","loading"), .cca_long_matrix(object$loadings_y,"Y","loading")))
  if (type == "cross_loadings") return(rbind(.cca_long_matrix(object$cross_loadings_x,"X","cross_loading"), .cca_long_matrix(object$cross_loadings_y,"Y","cross_loading")))
  if (type == "redundancy") return(rbind(
    data.frame(component=seq_along(object$redundancy_x),block="X",variance_extracted=object$variance_extracted_x,redundancy=object$redundancy_x),
    data.frame(component=seq_along(object$redundancy_y),block="Y",variance_extracted=object$variance_extracted_y,redundancy=object$redundancy_y)))
  if (type == "communality") return(rbind(data.frame(block="X",variable=names(object$communality_x),communality=object$communality_x),
                                           data.frame(block="Y",variable=names(object$communality_y),communality=object$communality_y)))
  cca_diagnose(object)$data
}

.cca_markdown_table <- function(d) {
  vals <- lapply(d, function(z) if (is.numeric(z)) format(round(z, 5), trim=TRUE) else as.character(z))
  d2 <- as.data.frame(vals, stringsAsFactors=FALSE)
  header <- paste0("| ", paste(names(d2), collapse=" | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(d2)), collapse=" | "), " |")
  rows <- apply(d2,1,function(z) paste0("| ",paste(z,collapse=" | ")," |"))
  paste(c(header,sep,rows),collapse="\n")
}

.cca_latex_table <- function(d) {
  align <- paste(rep("l", ncol(d)), collapse="")
  hdr <- paste(names(d), collapse=" & ")
  rows <- apply(d,1,function(z) paste(z,collapse=" & "))
  paste(c(paste0("\\begin{tabular}{",align,"}"),"\\hline",paste0(hdr," \\\\"),"\\hline",
          paste0(rows," \\\\"),"\\hline","\\end{tabular}"),collapse="\n")
}

#' Standardized publication tables and export
#' @export
cca_table <- function(object, type = "correlations",
                      format = c("data.frame","markdown","latex","html","csv","word"), file = NULL) {
  format <- match.arg(format)
  if (inherits(object,"cca_fit")) d <- .cca_table_fit(object,type)
  else if (inherits(object,"cca_bootstrap")) d <- object$correlation_ci
  else if (inherits(object,"cca_validation")) d <- object$results
  else if (inherits(object,"cca_permutation")) d <- data.frame(component=seq_along(object$observed),observed=object$observed,p_value=object$p_value,p_maxT=object$p_maxT %||% NA_real_)
  else if (inherits(object,"cca_robust_comparison")) d <- object$comparison
  else if (inherits(object,"cca_comparison")) d <- object$table
  else if (inherits(object,"cca_stability")) d <- object$table
  else if (inherits(object,"cca_tuning")) d <- object$results
  else .cca_stop("Unsupported object for cca_table().")
  if (format == "data.frame") return(d)
  if (format == "markdown") out <- .cca_markdown_table(d)
  if (format == "latex") out <- .cca_latex_table(d)
  if (format == "html") out <- paste0("<table><thead><tr>",paste0("<th>",names(d),"</th>",collapse=""),"</tr></thead><tbody>",
                                        paste(apply(d,1,function(z) paste0("<tr>",paste0("<td>",z,"</td>",collapse=""),"</tr>")),collapse=""),"</tbody></table>")
  if (format == "csv") {
    if (is.null(file)) { tf <- tempfile(fileext=".csv"); utils::write.csv(d,tf,row.names=FALSE); out <- paste(readLines(tf,warn=FALSE),collapse="\n"); unlink(tf); return(out) }
    utils::write.csv(d,file,row.names=FALSE); return(invisible(file))
  }
  if (format == "word") {
    .cca_require("flextable","cca_table(format='word')"); .cca_require("officer","cca_table(format='word')")
    ft <- flextable::flextable(d); ft <- flextable::autofit(ft)
    if (is.null(file)) return(ft)
    doc <- officer::read_docx(); doc <- flextable::body_add_flextable(doc,ft); print(doc,target=file); return(invisible(file))
  }
  if (!is.null(file)) { writeLines(out,file); return(invisible(file)) }
  out
}
