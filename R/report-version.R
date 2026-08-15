# Reporting and engine audit ------------------------------------------------

#' Report availability and versions of optional CCA engines
#' @export
cca_version_info <- function() {
  pkgs <- c("CCA","CCP","PMA","RGCCA","mixedCCA","multiCCA","corpcor","robustbase","MASS",
            "ggplot2","plotly","tourr","spinifex","gganimate","rstan","posterior","reticulate","caret","flextable","officer","rmarkdown")
  ver <- vapply(pkgs,function(p) if(requireNamespace(p,quietly=TRUE)) as.character(utils::packageVersion(p)) else NA_character_,character(1))
  out <- data.frame(package=pkgs,available=!is.na(ver),version=ver,stringsAsFactors=FALSE)
  py_available <- FALSE; py_version <- NA_character_
  if (requireNamespace("reticulate", quietly = TRUE)) {
    has_py <- tryCatch(isTRUE(reticulate::py_available(initialize = FALSE)), error = function(e) FALSE)
    if (has_py) {
      py_available <- tryCatch(reticulate::py_module_available("cca_zoo"), error = function(e) FALSE)
      if (isTRUE(py_available)) {
        py_version <- tryCatch({
          md <- reticulate::import("importlib.metadata", convert = TRUE)
          as.character(md$version("cca-zoo"))
        }, error = function(e) NA_character_)
      }
    }
  }
  rbind(out, data.frame(package = "python:cca-zoo", available = py_available,
                        version = py_version, stringsAsFactors = FALSE))
}

.cca_report_markdown <- function(object, bootstrap=NULL, validation=NULL, permutation=NULL, title="Canonical Correlation Analysis Report") {
  cor_tab <- .cca_markdown_table(cca_table(object,"correlations"))
  red_tab <- .cca_markdown_table(cca_table(object,"redundancy"))
  load_tab <- cca_table(object,"loadings"); load_tab <- load_tab[order(abs(load_tab$value),decreasing=TRUE),]; load_tab <- head(load_tab,30)
  load_md <- .cca_markdown_table(load_tab)
  scope <- object$inference_scope
  txt <- c(
    paste0("# ",title),"",
    "## Analysis identity","",
    paste0("- Method: **",object$method,"**"),
    paste0("- Engine: **",object$engine,"**"),
    paste0("- Paired observations: ",object$n),
    paste0("- Variables: X = ",object$p,"; Y = ",object$q),
    paste0("- Inferential scope: ",scope),"",
    "## Canonical association","",cor_tab,"",
    "The canonical correlations quantify association between fitted canonical variates. They should be interpreted together with redundancy, stability and out-of-sample performance rather than as stand-alone evidence of scientific relevance.","",
    "## Variance extracted and redundancy","",red_tab,"",
    "## Leading structure coefficients","",load_md,"",
    "## Interpretation safeguards","",
    "- Canonical weights are not interpreted in isolation when predictors are collinear.",
    "- Regularized, sparse, resistant, kernel and deep estimators require resampling or method-specific inference rather than automatic reuse of classical Wilks/Bartlett tests.",
    "- Outliers are flagged for investigation; they are not removed automatically.",
    "- Training correlation is separated from generalization performance."
  )
  if(!is.null(bootstrap)) txt <- c(txt,"","## Bootstrap uncertainty","",.cca_markdown_table(bootstrap$correlation_ci))
  if(!is.null(validation) && nrow(validation$results)) {
    sm_signed <- aggregate(correlation ~ component, validation$results,
                           function(z) c(mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE), n = sum(is.finite(z))))
    sm_abs <- aggregate(absolute_correlation ~ component, validation$results,
                        function(z) c(mean_abs = mean(z, na.rm = TRUE)))
    txt <- c(txt,"","## Out-of-sample validation","",
             "Signed held-out correlations preserve the orientation learned in each training fold; their absolute magnitudes are reported separately.","",
             "### Signed correlations","",paste(capture.output(print(sm_signed,row.names=FALSE)),collapse="\n"),"",
             "### Absolute magnitudes","",paste(capture.output(print(sm_abs,row.names=FALSE)),collapse="\n"))
  }
  if(!is.null(permutation)) {
    pt <- data.frame(component=seq_along(permutation$observed),observed=permutation$observed,p_value=permutation$p_value,p_maxT=permutation$p_maxT %||% NA_real_)
    txt <- c(txt,"","## Permutation assessment","",.cca_markdown_table(pt),"",permutation$note %||% "")
  }
  if(inherits(object,"cca_bayes") && !is.null(object$posterior_correlation_ci)) {
    d <- data.frame(component=seq_len(nrow(object$posterior_correlation_ci)),object$posterior_correlation_ci)
    txt <- c(txt,"","## Bayesian posterior canonical correlations","",.cca_markdown_table(d))
  }
  paste(txt,collapse="\n")
}

#' Generate a reproducible CCA analysis report
#' @export
cca_report <- function(object, file = "cca_report.md",
                       format = c("markdown","html","word"),
                       bootstrap = NULL, validation = NULL, permutation = NULL,
                       title = "Canonical Correlation Analysis Report",
                       include_figures = TRUE, quiet = TRUE) {
  if(!inherits(object,"cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  format <- match.arg(format)
  md <- .cca_report_markdown(object,bootstrap,validation,permutation,title)
  outdir <- dirname(file); if(!dir.exists(outdir)) dir.create(outdir,recursive=TRUE)
  if(include_figures) {
    .cca_require("ggplot2","cca_report(include_figures=TRUE)")
    figdir <- file.path(outdir,paste0(tools::file_path_sans_ext(basename(file)),"_figures")); dir.create(figdir,showWarnings=FALSE,recursive=TRUE)
    figs <- list(scree=cca_plot(object,"scree",bootstrap=bootstrap,permutation=permutation),scores=cca_plot(object,"scores"),redundancy=cca_plot(object,"redundancy"))
    links <- character()
    for(nm in names(figs)) {
      fn <- file.path(figdir,paste0(nm,".png")); cca_export(figs[[nm]],fn,preset="double",dpi=600)
      links <- c(links,paste0("### ",tools::toTitleCase(nm),"\n\n![](",basename(figdir),"/",basename(fn),")"))
    }
    md <- paste(md,"\n\n## Figures\n\n",paste(links,collapse="\n\n"))
  }
  if(format=="markdown") { writeLines(md,file,useBytes=TRUE); return(invisible(normalizePath(file,winslash="/",mustWork=FALSE))) }
  .cca_require("rmarkdown","cca_report rendering")
  ext <- if(format=="html") ".html" else ".docx"
  target <- if(tolower(tools::file_ext(file)) %in% c("html","docx")) file else paste0(tools::file_path_sans_ext(file),ext)
  rmd <- tempfile(fileext=".Rmd")
  yaml <- c("---",paste0("title: \"",gsub('"','\\\\"',title),"\""),
            paste0("output: ",if(format=="html") "html_document" else "word_document"),"---","")
  writeLines(c(yaml,md),rmd,useBytes=TRUE)
  rmarkdown::render(rmd,output_file=basename(target),output_dir=dirname(target),quiet=quiet,envir=new.env(parent=globalenv()))
  invisible(normalizePath(target,winslash="/",mustWork=FALSE))
}
