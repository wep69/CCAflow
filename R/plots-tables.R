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

.cca_bipolar_palette <- function(pal, grayscale = FALSE) {
  if (grayscale) c("Positive (+)" = "#262626", "Negative (-)" = "#8c8c8c")
  else c("Positive (+)" = pal[1], "Negative (-)" = pal[2])
}

.cca_filter_diag_data <- function(diag_obj, block = "both", top = 20,
                                  order_by_component = NULL, rank_by = "loading") {
  d <- diag_obj$table
  if (block != "both") {
    d <- d[d$block == toupper(block), , drop = FALSE]
  }
  d$block <- droplevels(d$block)
  if (!nrow(d)) .cca_stop("No variables remain to plot after block filtering.")

  ## Keep the top variables overall per block so multi-component displays align rows
  if (is.finite(top) && length(unique(d$variable)) > top) {
    keep_vars <- character()
    for (blk in levels(d$block)) {
      db <- d[d$block == blk, , drop = FALSE]
      sc <- tapply(abs(db[[rank_by]]), db$variable, max, na.rm = TRUE)
      keep_vars <- c(keep_vars, names(sort(sc, decreasing = TRUE))[seq_len(min(length(sc), top))])
    }
    d <- d[d$variable %in% keep_vars, , drop = FALSE]
  }

  comps_avail <- unique(d$component_index)
  ref_k <- order_by_component %||% if (2L %in% comps_avail) 2L else comps_avail[1L]
  if (!ref_k %in% comps_avail) ref_k <- comps_avail[1L]
  d_ref <- d[d$component_index == ref_k, , drop = FALSE]
  ord <- names(sort(tapply(d_ref[[rank_by]], d_ref$label, mean, na.rm = TRUE)))
  d$label_ord <- factor(d$label, levels = unique(c(ord, unique(d$label))))
  d
}

.cca_plot_coeff_multi <- function(object, type, components, block, top,
                                  variable_labels, style, polarity_color, pal, grayscale) {
  comps <- .cca_resolve_components(object, components, max_default = 3L)
  if (type == "coefficients") .cca_assert_orig_weights(object, block, "type = 'coefficients'")
  dg <- cca_weight_diagnostics(object, components = comps, labels = variable_labels)
  col_nm <- switch(type, loadings = "loading", cross_loadings = "cross_loading", coefficients = "weight")
  d <- .cca_filter_diag_data(dg, block = block, top = top, order_by_component = comps[1L], rank_by = col_nm)
  d$val <- d[[col_nm]]
  d$pol <- factor(ifelse(d$val >= 0, "Positive (+)", "Negative (-)"),
                  levels = c("Positive (+)", "Negative (-)"))

  p <- if (isTRUE(polarity_color)) {
    ggplot2::ggplot(d, ggplot2::aes(x = .data[["val"]], y = .data[["label_ord"]], color = .data[["pol"]])) +
      ggplot2::scale_color_manual(values = .cca_bipolar_palette(pal, grayscale), name = "Polarity")
  } else {
    ggplot2::ggplot(d, ggplot2::aes(x = .data[["val"]], y = .data[["label_ord"]]))
  }
  p <- p + ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "#334155")
  if (style == "lollipop") {
    p <- p + ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data[["val"]], y = .data[["label_ord"]], yend = .data[["label_ord"]]),
      linewidth = 0.75, alpha = 0.85
    )
  }
  p <- p + ggplot2::geom_point(size = 2.2)
  if (block == "both" || length(comps) > 1L) {
    p <- p + ggplot2::facet_grid(block ~ component_label, scales = "free_y", space = "free_y")
  }
  x_lab <- switch(type, loadings = "Structure coefficient",
                  cross_loadings = "Cross-loading", coefficients = "Canonical coefficient")
  ttl <- paste(gsub("_", " ", type), if (block == "both") "X and Y" else toupper(block),
               paste(paste0("CV", comps), collapse = ", "))
  p + ggplot2::labs(x = x_lab, y = NULL, title = ttl) + theme_cca_publication()
}

.cca_plot_suppressor <- function(object, components, block, top, variable_labels,
                                 vif_cap, pal, grayscale) {
  .cca_assert_orig_weights(object, block, "type = 'suppressor'")
  comps <- .cca_resolve_components(object, components, max_default = 2L)
  dg <- cca_weight_diagnostics(object, components = comps, labels = variable_labels)
  d <- .cca_filter_diag_data(dg, block = block, top = top, order_by_component = comps[1L], rank_by = "loading")
  d$inv_fac <- factor(d$suppressor, levels = c(FALSE, TRUE))
  d$vif_plot <- pmin(d$vif, vif_cap)
  blk_cols <- if (grayscale) c("X" = "#262626", "Y" = "#737373") else c("X" = pal[1], "Y" = pal[2])
  quad_fill <- if (grayscale) "#999999" else "#dc2626"

  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[["loading"]], y = .data[["weight"]],
    color = .data[["block"]], shape = .data[["inv_fac"]]
  )) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#64748b", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#64748b", linewidth = 0.4) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#94a3b8", linewidth = 0.45) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, alpha = 0.07, fill = quad_fill) +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, alpha = 0.07, fill = quad_fill) +
    ggplot2::geom_point(ggplot2::aes(size = .data[["vif_plot"]]), alpha = 0.9) +
    ggplot2::scale_color_manual(values = blk_cols, name = "Block") +
    ggplot2::scale_shape_manual(
      values = c("FALSE" = 16, "TRUE" = 17), drop = FALSE,
      labels = c("FALSE" = "Aligned sign", "TRUE" = "Inverted sign (suppressor)"),
      name = "Sign relation"
    ) +
    ggplot2::scale_size_continuous(name = "VIF")

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = .data[["label"]]),
      size = 2.6, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
      min.segment.length = 0, seed = 42, show.legend = FALSE
    )
  } else {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data[["label"]]),
      size = 2.6, check_overlap = TRUE, vjust = -0.6, show.legend = FALSE
    )
  }
  if (length(comps) > 1L) p <- p + ggplot2::facet_wrap(~ component_label)
  p + ggplot2::labs(
    title = "Suppressor paradox and intra-block collinearity",
    subtitle = "Shaded quadrants indicate sign reversal between weight and structure loading; size reflects VIF",
    x = "Intra-block structure coefficient (r_s)",
    y = "Standardized canonical weight (w)"
  ) + theme_cca_publication()
}

.cca_plot_bipolar_lollipop <- function(object, components, block, top, variable_labels,
                                       thresholds, order_by_component, metric, pal, grayscale) {
  metric <- match.arg(metric, c("loadings", "cross_loadings", "coefficients"))
  if (metric == "coefficients") .cca_assert_orig_weights(object, block, "metric = 'coefficients'")
  comps <- .cca_resolve_components(object, components, max_default = 3L)
  dg <- cca_weight_diagnostics(object, components = comps, labels = variable_labels)
  col_nm <- switch(metric, loadings = "loading", cross_loadings = "cross_loading", coefficients = "weight")
  d <- .cca_filter_diag_data(dg, block = block, top = top,
                             order_by_component = order_by_component, rank_by = col_nm)
  d$val <- d[[col_nm]]
  d$pol <- factor(ifelse(d$val >= 0, "Positive (+)", "Negative (-)"),
                  levels = c("Positive (+)", "Negative (-)"))
  d$lbl_val <- sprintf("%+.2f", d$val)
  max_abs <- max(1.0, max(abs(d$val), na.rm = TRUE))
  nudge <- 0.055 * max_abs
  d$x_lbl <- d$val + ifelse(d$val >= 0, nudge, -nudge)
  d$hjust_lbl <- ifelse(d$val >= 0, 0, 1)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[["val"]], y = .data[["label_ord"]], color = .data[["pol"]])) +
    ggplot2::geom_vline(xintercept = 0, color = "#1e293b", linewidth = 0.6)
  if (length(thresholds) >= 1L && is.finite(thresholds[1L])) {
    p <- p + ggplot2::geom_vline(xintercept = c(-thresholds[1L], thresholds[1L]),
                                 linetype = "dashed", color = "#64748b", linewidth = 0.45)
  }
  if (length(thresholds) >= 2L && is.finite(thresholds[2L])) {
    p <- p + ggplot2::geom_vline(xintercept = c(-thresholds[2L], thresholds[2L]),
                                 linetype = "dotted", color = "#94a3b8", linewidth = 0.45)
  }
  p <- p +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data[["val"]], y = .data[["label_ord"]], yend = .data[["label_ord"]]),
      linewidth = 0.9, alpha = 0.88
    ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data[["x_lbl"]], label = .data[["lbl_val"]], hjust = .data[["hjust_lbl"]]),
      size = 2.5, fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::facet_grid(block ~ component_label, scales = "free_y", space = "free_y") +
    ggplot2::scale_x_continuous(limits = c(-max_abs * 1.24, max_abs * 1.28)) +
    ggplot2::scale_color_manual(values = .cca_bipolar_palette(pal, grayscale), name = "Polarity") +
    ggplot2::labs(
      title = "Bipolar variable participation across canonical pairs",
      subtitle = sprintf("Dashed/dotted cutoffs at +/-%s", paste(format(thresholds, trim = TRUE), collapse = " and +/-")),
      x = switch(metric, loadings = "Structure coefficient (r_s)",
                 cross_loadings = "Cross-loading (r_c)", coefficients = "Standardized weight (w)"),
      y = NULL
    ) +
    theme_cca_publication() +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.major.y = ggplot2::element_line(color = "#eef2f6", linewidth = 0.35))
  p
}

.cca_plot_tornado <- function(object, components, block, top, variable_labels,
                              order_by_component, pal, grayscale) {
  .cca_assert_orig_weights(object, block, "type = 'weight_vs_loading'")
  comps <- .cca_resolve_components(object, components, max_default = 3L)
  dg <- cca_weight_diagnostics(object, components = comps, labels = variable_labels)
  d <- .cca_filter_diag_data(dg, block = block, top = top,
                             order_by_component = order_by_component, rank_by = "loading")
  d_long <- rbind(
    data.frame(label_ord = d$label_ord, block = d$block, component_label = d$component_label,
               value = d$loading, metric = "Structure loading (r_s)", stringsAsFactors = FALSE),
    data.frame(label_ord = d$label_ord, block = d$block, component_label = d$component_label,
               value = d$weight, metric = "Standardized weight (w)", stringsAsFactors = FALSE)
  )
  d_long$metric <- factor(d_long$metric, levels = c("Standardized weight (w)", "Structure loading (r_s)"))
  d_long$fill_group <- factor(
    paste0(ifelse(d_long$metric == "Structure loading (r_s)", "r_s ", "w "),
           ifelse(d_long$value >= 0, "(+)", "(-)")),
    levels = c("r_s (+)", "r_s (-)", "w (+)", "w (-)")
  )
  fill_cols <- if (grayscale) {
    c("r_s (+)" = "#1f2937", "r_s (-)" = "#6b7280", "w (+)" = "#9ca3af", "w (-)" = "#d1d5db")
  } else {
    c("r_s (+)" = pal[1], "r_s (-)" = pal[2], "w (+)" = pal[6], "w (-)" = pal[5])
  }

  lim <- max(1.05, max(abs(d_long$value), na.rm = TRUE))
  pad <- 0.08 * lim
  d_inv <- d[d$suppressor, , drop = FALSE]
  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, color = "#0f172a", linewidth = 0.6) +
    ggplot2::geom_col(
      data = d_long,
      ggplot2::aes(x = .data[["value"]], y = .data[["label_ord"]],
                   fill = .data[["fill_group"]], group = .data[["metric"]]),
      position = ggplot2::position_dodge(width = 0.72), width = 0.66, alpha = 0.92
    )
  if (nrow(d_inv)) {
    dom <- ifelse(abs(d_inv$weight) > abs(d_inv$loading), d_inv$weight, d_inv$loading)
    sgn <- ifelse(dom >= 0, 1, -1)
    d_inv$mark_x <- dom + sgn * pad
    d_inv$txt_x <- dom + sgn * (pad * 1.75)
    d_inv$txt_hjust <- ifelse(sgn > 0, 0, 1)
    d_inv$sup_lbl <- sprintf("sup: w=%+.2f / rs=%+.2f", d_inv$weight, d_inv$loading)
    pt_col <- if (grayscale) "#000000" else "#dc2626"
    p <- p +
      ggplot2::geom_point(
        data = d_inv,
        ggplot2::aes(x = .data[["mark_x"]], y = .data[["label_ord"]]),
        shape = 18, size = 3.0, color = pt_col, inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = d_inv,
        ggplot2::aes(x = .data[["txt_x"]], y = .data[["label_ord"]],
                     hjust = .data[["txt_hjust"]], label = .data[["sup_lbl"]]),
        size = 2.3, color = pt_col, fontface = "bold", inherit.aes = FALSE
      )
  }
  p +
    ggplot2::facet_grid(block ~ component_label, scales = "free_y", space = "free_y") +
    ggplot2::scale_x_continuous(limits = c(-lim * 1.45, lim * 1.48)) +
    ggplot2::scale_fill_manual(values = fill_cols, name = "Coefficient and polarity") +
    ggplot2::labs(
      title = "Bidirectional mirrored comparison: standardized weight (w) vs structure loading (r_s)",
      subtitle = "Diamonds highlight sign-reversed suppressor variables between weight and structure loading",
      x = "Standardized coefficient magnitude and polarity", y = NULL
    ) +
    theme_cca_publication() +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.major.y = ggplot2::element_line(color = "#eef2f6", linewidth = 0.35))
}

.cca_plot_synoptic_heatmap <- function(object, components, block, top, variable_labels,
                                       squish, pal, grayscale) {
  .cca_assert_orig_weights(object, block, "type = 'synoptic_heatmap'")
  comps <- .cca_resolve_components(object, components, max_default = 4L)
  dg <- cca_weight_diagnostics(object, components = comps, labels = variable_labels)
  d <- .cca_filter_diag_data(dg, block = block, top = top, order_by_component = comps[1L], rank_by = "loading")
  d$label_rev <- factor(d$label, levels = rev(levels(d$label_ord)))
  dh <- rbind(
    data.frame(label_rev = d$label_rev, block = d$block, component = d$component,
               metric = "1. Structure loading (r_s)", value = d$loading, highlight = FALSE, stringsAsFactors = FALSE),
    data.frame(label_rev = d$label_rev, block = d$block, component = d$component,
               metric = "2. Cross-loading (r_c)", value = d$cross_loading, highlight = FALSE, stringsAsFactors = FALSE),
    data.frame(label_rev = d$label_rev, block = d$block, component = d$component,
               metric = "3. Standardized weight (w)", value = d$weight,
               highlight = d$suppressor | d$collinear_inflation, stringsAsFactors = FALSE)
  )
  dh$metric <- factor(dh$metric, levels = c("1. Structure loading (r_s)",
                                            "2. Cross-loading (r_c)",
                                            "3. Standardized weight (w)"))
  dh$fill_squish <- pmax(pmin(dh$value, squish), -squish)
  dh$txt_col <- ifelse(abs(dh$fill_squish) >= 0.55 * squish, "white", "#0f172a")
  dh$tile_lbl <- sprintf("%+.2f%s", dh$value, ifelse(dh$highlight, "*", ""))
  low_col <- if (grayscale) "#595959" else pal[2]
  high_col <- if (grayscale) "#1a1a1a" else pal[1]

  ggplot2::ggplot(dh, ggplot2::aes(x = .data[["component"]], y = .data[["label_rev"]], fill = .data[["fill_squish"]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_tile(
      data = dh[dh$highlight, , drop = FALSE],
      ggplot2::aes(x = .data[["component"]], y = .data[["label_rev"]]),
      fill = NA, color = "#0f172a", linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[["tile_lbl"]], color = .data[["txt_col"]]),
      size = 2.6, fontface = "bold"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(
      low = low_col, mid = "#f8fafc", high = high_col, midpoint = 0,
      limits = c(-squish, squish),
      name = sprintf("Value (color scale clipped at +/-%.2f)", squish)
    ) +
    ggplot2::facet_grid(block ~ metric, scales = "free_y", space = "free_y") +
    ggplot2::labs(
      title = "Synoptic multidimensional participation map",
      subtitle = "Asterisk (*) and dashed outline mark sign reversal or collinear weight inflation",
      x = "Canonical pair", y = NULL
    ) +
    theme_cca_publication() +
    ggplot2::theme(legend.position = "bottom")
}

.cca_plot_heliograph <- function(object, component, block, top, variable_labels,
                                 thresholds, pal, grayscale) {
  k <- .cca_resolve_components(object, component, max_default = 1L)[1L]
  dg <- cca_weight_diagnostics(object, components = k, labels = variable_labels)
  d <- .cca_filter_diag_data(dg, block = block, top = top, order_by_component = k, rank_by = "loading")
  n_v <- nrow(d)
  d$ord_id <- seq_len(n_v)
  angs <- 90 - 360 * (d$ord_id - 0.5) / n_v
  d$ang_text <- ifelse(angs < -90, angs + 180, angs)
  d$hjust_text <- ifelse(angs < -90, 1, 0)
  d$abs_val <- abs(d$loading)
  d$lbl_pos <- pmax(d$abs_val + 0.08, 0.38)
  d$lbl_polar <- sprintf("%s (%+.2f)", d$label, d$loading)
  pol_cols <- .cca_bipolar_palette(pal, grayscale)

  ggplot2::ggplot(d, ggplot2::aes(x = .data[["ord_id"]], y = .data[["abs_val"]],
                                  fill = .data[["polarity_loading"]])) +
    ggplot2::geom_hline(yintercept = thresholds[1L], linetype = "dashed", color = "#475569", linewidth = 0.45) +
    ggplot2::geom_hline(yintercept = c(0.60, 0.90), linetype = "dotted", color = "#94a3b8", linewidth = 0.4) +
    ggplot2::geom_col(width = 0.78, alpha = 0.92, color = "white", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(y = .data[["lbl_pos"]], label = .data[["lbl_polar"]],
                   angle = .data[["ang_text"]], hjust = .data[["hjust_text"]],
                   color = .data[["polarity_loading"]]),
      size = 2.4, fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::coord_polar(start = 0, clip = "off") +
    ggplot2::scale_x_continuous(limits = c(0.5, n_v + 0.5), breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1.48), breaks = c(thresholds[1L], 0.6, 0.9)) +
    ggplot2::scale_fill_manual(values = pol_cols, name = "Polarity (r_s)") +
    ggplot2::scale_color_manual(values = pol_cols, guide = "none") +
    ggplot2::labs(title = sprintf("Polar heliograph CV%d (rho = %.3f)", k, object$correlations[k]),
                  x = NULL, y = NULL) +
    theme_cca_publication(base_size = 8.8) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(), legend.position = "bottom")
}

.cca_plot_bipolar_balance <- function(object, components, block, pal, grayscale) {
  comps <- .cca_resolve_components(object, components, max_default = 4L)
  dg <- cca_weight_diagnostics(object, components = comps)
  sm <- dg$summary
  if (block != "both") sm <- sm[sm$block == toupper(block), , drop = FALSE]
  d_en <- rbind(
    data.frame(component = sm$component, block = sm$block,
               polarity = "Positive (+)", energy = sm$energy_positive, stringsAsFactors = FALSE),
    data.frame(component = sm$component, block = sm$block,
               polarity = "Negative (-)", energy = -sm$energy_negative, stringsAsFactors = FALSE)
  )
  comp_lvls <- rev(paste0("CV", sort(unique(comps))))
  d_en$component <- factor(d_en$component, levels = comp_lvls)
  d_en$block <- factor(d_en$block, levels = unique(sm$block))
  d_en$polarity <- factor(d_en$polarity, levels = c("Positive (+)", "Negative (-)"))
  d_en$lbl <- sprintf("%+.2f", d_en$energy)
  d_en$hjust_lbl <- ifelse(d_en$energy >= 0, -0.15, 1.15)
  max_pos <- max(c(0.5, d_en$energy[d_en$energy >= 0]), na.rm = TRUE)
  min_neg <- min(c(-0.5, d_en$energy[d_en$energy < 0]), na.rm = TRUE)

  ggplot2::ggplot(d_en, ggplot2::aes(x = .data[["energy"]], y = .data[["component"]], fill = .data[["polarity"]])) +
    ggplot2::geom_vline(xintercept = 0, color = "#0f172a", linewidth = 0.6) +
    ggplot2::geom_col(width = 0.62, alpha = 0.92) +
    ggplot2::geom_text(
      data = d_en[abs(d_en$energy) > 0.02, , drop = FALSE],
      ggplot2::aes(label = .data[["lbl"]], hjust = .data[["hjust_lbl"]]),
      size = 2.5, fontface = "bold", color = "#1e293b", show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~ block, ncol = 1) +
    ggplot2::scale_x_continuous(limits = c(min_neg * 1.35, max_pos * 1.22)) +
    ggplot2::scale_fill_manual(values = .cca_bipolar_palette(pal, grayscale), name = "Polarity (r_s)") +
    ggplot2::labs(
      title = "Bipolar variance balance (sum of squared loadings)",
      x = "Signed sum of squared structure coefficients", y = NULL
    ) +
    theme_cca_publication(base_size = 8.8) +
    ggplot2::theme(legend.position = "bottom")
}

.cca_plot_heliograph_panel <- function(object, components, block, top, variable_labels,
                                       thresholds, pal, grayscale) {
  .cca_require("patchwork", "cca_plot(type = 'heliograph_panel')")
  m <- length(object$correlations)
  comps_h <- if (is.null(components)) {
    if (m >= 3L) c(2L, 3L) else if (m == 2L) c(1L, 2L) else 1L
  } else {
    .cca_resolve_components(object, components, max_default = 2L)
  }
  comps_h <- head(comps_h, 2L)
  helio_list <- lapply(comps_h, function(k) {
    .cca_plot_heliograph(object, component = k, block = block, top = top,
                         variable_labels = variable_labels, thresholds = thresholds,
                         pal = pal, grayscale = grayscale)
  })
  bal_comps <- seq_len(min(4L, m))
  p_bal <- .cca_plot_bipolar_balance(object, components = bal_comps, block = block,
                                     pal = pal, grayscale = grayscale)
  all_plots <- c(helio_list, list(p_bal))
  w_rel <- c(rep(1.22, length(helio_list)), 0.92)
  patchwork::wrap_plots(all_plots, widths = w_rel, guides = "collect") +
    patchwork::plot_annotation(
      title = "Bipolar polar heliographs and quadratic participation balance"
    ) & ggplot2::theme(legend.position = "bottom")
}

.cca_plot_fit <- function(object, type, component = 1, component_y = 2,
                          components = NULL, block = c("x", "y", "both"), top = 20,
                          threshold = 0.3, thresholds = c(0.30, 0.45),
                          groups = NULL, labels = NULL, variable_labels = NULL,
                          bootstrap = NULL, permutation = NULL, grayscale = FALSE,
                          group_ellipse = FALSE, ellipse_level = 0.95,
                          style = c("dot", "lollipop"), polarity_color = FALSE,
                          sparsity_mode = c("binary", "signed", "magnitude"),
                          vif_cap = 28, squish = 1.05, order_by_component = NULL,
                          metric = "loadings", ...) {
  .cca_require("ggplot2", "cca_plot")
  block_missing <- missing(block)
  block <- match.arg(block)
  style <- match.arg(style)
  sparsity_mode <- match.arg(sparsity_mode)
  pal <- cca_palette(8, grayscale)
  type <- match.arg(type, c("scree", "scores", "circle", "dual_circle", "loadings",
                            "cross_loadings", "coefficients", "biplot", "redundancy",
                            "heatmap", "network", "sparsity", "influence",
                            "communality", "suppressor", "bipolar_lollipop",
                            "weight_vs_loading", "tornado", "synoptic_heatmap",
                            "heliograph", "bipolar_balance", "heliograph_panel"))
  if (block_missing && type %in% c("suppressor", "bipolar_lollipop",
                                   "weight_vs_loading", "tornado",
                                   "synoptic_heatmap", "heliograph",
                                   "bipolar_balance", "heliograph_panel")) {
    block <- "both"
  }
  ## Companion objects feed the uncertainty layers of the scree figure only;
  ## accepting them for other types in silence hides a no-op argument.
  ignorados <- c(if (!is.null(bootstrap)) "bootstrap", if (!is.null(permutation)) "permutation")
  if (length(ignorados) && type != "scree")
    warning(sprintf(paste0("%s = is only used by type = 'scree'; it is ignored for ",
                    "type = '%s'. Call cca_plot(object, 'scree', %s = ...) for the ",
                    "uncertainty figure."),
                    paste(ignorados, collapse = " and "), type, ignorados[1]), call. = FALSE)
  if (type == "circle" && block == "both") type <- "dual_circle"
  if (type %in% c("weight_vs_loading", "tornado")) {
    return(.cca_plot_tornado(object, components, block, top, variable_labels,
                             order_by_component, pal, grayscale))
  }
  if (type == "suppressor") {
    return(.cca_plot_suppressor(object, components, block, top, variable_labels,
                                vif_cap, pal, grayscale))
  }
  if (type == "bipolar_lollipop") {
    return(.cca_plot_bipolar_lollipop(object, components, block, top, variable_labels,
                                      thresholds, order_by_component, metric, pal, grayscale))
  }
  if (type == "synoptic_heatmap") {
    return(.cca_plot_synoptic_heatmap(object, components, block, top, variable_labels,
                                      squish, pal, grayscale))
  }
  if (type == "heliograph") {
    k_h <- if (!is.null(components)) components[1L] else if (!missing(component)) component else min(2L, length(object$correlations))
    return(.cca_plot_heliograph(object, k_h, block, top, variable_labels, thresholds, pal, grayscale))
  }
  if (type == "bipolar_balance") {
    return(.cca_plot_bipolar_balance(object, components, block, pal, grayscale))
  }
  if (type == "heliograph_panel") {
    return(.cca_plot_heliograph_panel(object, components, block, top, variable_labels, thresholds, pal, grayscale))
  }
  if (type == "scree") {
    d <- data.frame(component = seq_along(object$correlations), correlation = object$correlations)
    p <- ggplot2::ggplot(d, ggplot2::aes(component, correlation)) +
      ggplot2::geom_line(linewidth = 0.5) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_x_continuous(breaks = d$component) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = "Canonical dimension", y = "Canonical correlation", title = "Canonical correlation profile")
    if (!is.null(bootstrap) && inherits(bootstrap, "cca_bootstrap")) {
      ci <- bootstrap$correlation_ci
      ## The CI table carries 'estimate'; the base plot maps 'correlation'.
      ## Declare the aesthetics without inheriting, so the layer draws.
      if (!"correlation" %in% names(ci) && "estimate" %in% names(ci)) ci$correlation <- ci$estimate
      p <- p + ggplot2::geom_errorbar(data = ci, inherit.aes = FALSE,
                                      ggplot2::aes(x = component, ymin = lower, ymax = upper),
                                      width = 0.12)
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
    if (block == "both" || !is.null(components) || style != "dot" ||
        isTRUE(polarity_color) || !is.null(variable_labels)) {
      return(.cca_plot_coeff_multi(object, type, components %||% component, block, top,
                                   variable_labels, style, polarity_color, pal, grayscale))
    }
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
  if (type == "communality") {
    ## Same quantities cca_table(type = "communality") reports, drawn per variable.
    d <- rbind(data.frame(block = "X", variable = names(object$communality_x),
                         communality = as.numeric(object$communality_x), row.names = NULL),
               data.frame(block = "Y", variable = names(object$communality_y),
                         communality = as.numeric(object$communality_y), row.names = NULL))
    d <- d[is.finite(d$communality), , drop = FALSE]
    if (!nrow(d)) .cca_stop("Communalities are unavailable for this fit.")
    d$variable <- factor(d$variable, levels = d$variable[order(d$communality)])
    return(ggplot2::ggplot(d, ggplot2::aes(communality, variable, fill = block)) +
      ggplot2::geom_col(width = 0.7) + ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_fill_manual(values = c("X" = pal[1], "Y" = pal[2]), guide = "none") +
      ggplot2::labs(x = "Communality (sum of squared loadings)", y = NULL,
                    title = "Variable communalities") + theme_cca_publication())
  }
  if (type == "sparsity") {
    .cca_assert_orig_weights(object, "both", "type = 'sparsity'")
    d <- rbind(.cca_long_matrix(object$weights_x, "X", "weight"), .cca_long_matrix(object$weights_y, "Y", "weight"))
    d$component <- factor(d$component, levels = colnames(object$weights_x))
    if (sparsity_mode == "binary") {
      d$selected <- abs(d$value) > 1e-8
      return(ggplot2::ggplot(d, ggplot2::aes(component, variable, fill = selected)) + ggplot2::geom_tile() +
        ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_fill_manual(values = c("FALSE" = "white", "TRUE" = pal[1])) +
        ggplot2::labs(x = "Component", y = NULL, fill = "Selected", title = "Canonical sparsity map") + theme_cca_publication())
    }
    if (sparsity_mode == "signed") {
      lim <- max(abs(d$value), 1e-6, na.rm = TRUE)
      return(ggplot2::ggplot(d, ggplot2::aes(component, variable, fill = .data[["value"]])) +
        ggplot2::geom_tile(color = "white", linewidth = 0.35) +
        ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_fill_gradient2(low = pal[2], mid = "white", high = pal[1], midpoint = 0, limits = c(-lim, lim)) +
        ggplot2::labs(x = "Component", y = NULL, fill = "Weight", title = "Canonical signed weight sparsity") + theme_cca_publication())
    }
    d$abs_val <- abs(d$value)
    return(ggplot2::ggplot(d, ggplot2::aes(component, variable, fill = .data[["abs_val"]])) +
      ggplot2::geom_tile(color = "white", linewidth = 0.35) +
      ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_fill_gradient(low = "white", high = pal[1]) +
      ggplot2::labs(x = "Component", y = NULL, fill = "|Weight|", title = "Canonical weight magnitude map") + theme_cca_publication())
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

.cca_plot_other <- function(object, type, component = 1, block = c("x", "y", "both"),
                            top = 20, quantity = NULL, highlight_zero = FALSE,
                            grayscale = FALSE, ...) {
  .cca_require("ggplot2", "cca_plot")
  block <- match.arg(block)
  pal <- cca_palette(8, grayscale)
  if (inherits(object, "cca_bootstrap")) {
    ## Same rule as plot.cca_bootstrap(): an explicit quantity wins over type.
    if (!is.null(quantity)) type <- quantity
    if (type %in% c("bootstrap", "scree")) {
      d <- object$correlation_ci
      return(ggplot2::ggplot(d, ggplot2::aes(component, estimate)) + ggplot2::geom_point(size=2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=lower,ymax=upper), width=.12) +
        ggplot2::scale_x_continuous(breaks=d$component) + ggplot2::coord_cartesian(ylim=c(0,1)) +
        ggplot2::labs(x="Canonical dimension",y="Canonical correlation",title="Bootstrap uncertainty") + theme_cca_publication())
    }
    quantity <- if (type %in% c("loading","cross_loading","weight")) type else "loading"
    d <- object$variable_draws
    blk_sel <- if (block == "both") c("x", "y") else tolower(block)
    d <- d[d$quantity == quantity & tolower(d$block) %in% blk_sel & d$component == paste0("CV",component),,drop=FALSE]
    if (!nrow(d)) .cca_stop("Requested bootstrap quantity is unavailable.")
    summarize_block <- function(db, blk_lbl) {
      sp <- split(db$value, db$variable)
      sm <- do.call(rbind, lapply(sp, function(z) {
        q <- unname(stats::quantile(z, c(.025, .975), na.rm = TRUE))
        c(median = stats::median(z, na.rm = TRUE), lower = q[1], upper = q[2])
      }))
      dd <- data.frame(block = blk_lbl, variable = rownames(sm), median = sm[, "median"],
                       lower = sm[, "lower"], upper = sm[, "upper"], row.names = NULL,
                       stringsAsFactors = FALSE)
      idx <- .cca_top(dd$median, top)
      dd[idx, , drop = FALSE]
    }
    if (block == "both") {
      dd <- do.call(rbind, lapply(c("x", "y"), function(b) {
        db <- d[tolower(d$block) == b, , drop = FALSE]
        if (nrow(db)) summarize_block(db, toupper(b)) else NULL
      }))
      if (!nrow(dd)) .cca_stop("Requested bootstrap quantity is unavailable.")
      dd$block <- factor(dd$block, levels = unique(dd$block))
    } else {
      dd <- summarize_block(d, toupper(block))
    }
    dd$variable <- factor(dd$variable, levels = unique(dd$variable[order(dd$median)]))
    dd$crosses_zero <- factor(ifelse(dd$lower <= 0 & dd$upper >= 0, "Crosses 0", "Consistent sign"),
                              levels = c("Consistent sign", "Crosses 0"))
    p <- if (isTRUE(highlight_zero)) {
      ggplot2::ggplot(dd, ggplot2::aes(median, variable, color = .data[["crosses_zero"]])) +
        ggplot2::scale_color_manual(values = c("Consistent sign" = pal[1], "Crosses 0" = pal[2]), name = "CI vs zero")
    } else {
      ggplot2::ggplot(dd, ggplot2::aes(median, variable))
    }
    p <- p + ggplot2::geom_vline(xintercept=0,linewidth=.3) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin=lower,xmax=upper), orientation="y", width=.15) + ggplot2::geom_point(size=2) +
      ggplot2::labs(x=paste("Bootstrap",gsub("_"," ",quantity)),y=NULL,
                    title=paste("Bootstrap", if (block == "both") "X and Y" else toupper(block), paste0("CV",component)))
    if (block == "both") p <- p + ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y")
    return(p + theme_cca_publication())
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
  ## Inference and validation objects are served by the same figure helpers as
  ## the plot() methods; the mapping is done on type only, never by calling
  ## plot() from here (that would loop back into cca_plot()).
  rotas <- c(cca_bootstrap = "bootstrap", cca_validation = "validation",
             cca_permutation = "permutation", cca_comparison = "comparison",
             cca_robust_comparison = "comparison", cca_tuning = "tuning",
             cca_stability = "stability")
  chaves_inferencia <- c("bootstrap", "validation", "permutation", "comparison",
                         "tuning", "stability")
  if (inherits(object, "cca_weight_diagnostics")) {
    if (identical(type, "scree")) type <- "suppressor"
    object <- object$fit
  }
  if (!inherits(object, "cca_fit")) {
    if (identical(type, "scree")) {
      for (cl in names(rotas)) if (inherits(object, cl)) { type <- unname(rotas[[cl]]); break }
    }
  } else if (type %in% chaves_inferencia) {
    .cca_stop(sprintf(paste0("type '%s' is an inference figure, not a figure of a fitted model. ",
                            "Build the object first and plot it: ",
                            "plot(cca_%s(fit, ...)) for bootstrap/validation/permutation, ",
                            "or plot(cca_compare(...)), plot(cca_tune(...)), plot(cca_stability(...))."),
             type, if (type %in% c("validation", "permutation")) type else "bootstrap"))
  }
  p <- if (inherits(object, "cca_fit")) .cca_plot_fit(object, type = type, ...) else .cca_plot_other(object, type = type, ...)
  if (isTRUE(interactive)) {
    if (type %in% c("heliograph", "heliograph_panel")) {
      .cca_stop(sprintf(
        "interactive = TRUE is not supported for type = '%s' because polar/composite layouts cannot be converted reliably to plotly. Use interactive = FALSE or inspect type = 'bipolar_lollipop' / 'synoptic_heatmap' interactively.",
        type
      ))
    }
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
  if (inherits(plot, "plotly"))
    .cca_stop(paste0("cca_export() writes static ggplot figures. A plotly widget cannot be ",
                     "saved as a static file; build the figure with cca_plot(interactive = FALSE) ",
                     "or keep the widget for HTML output."))
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
  type <- match.arg(type, c("correlations","coefficients","loadings","cross_loadings","redundancy","communality","diagnostics","participation"))
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
  if (type == "participation") return(cca_weight_diagnostics(object)$table)
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
