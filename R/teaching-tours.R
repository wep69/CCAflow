# Teaching, simulation and projection tours ---------------------------------

#' Simulate controlled canonical-correlation scenarios
#' @export
cca_simulate <- function(n = 150, p = 6, q = 5,
                         scenario = c("linear","weak","null","collinear","highdim","sparse",
                                      "outliers","heavy_tail","monotone","nonmonotone","binary",
                                      "zero_inflated","missing","functional","three_blocks"),
                         rho = 0.8, seed = 1) {
  scenario <- match.arg(scenario)
  .cca_with_seed(seed, {
    if (scenario == "highdim") { p <- max(p, ceiling(n * .8)); q <- max(q, ceiling(n * .7)) }
    z <- matrix(stats::rnorm(n * 3), n, 3)
    Ax <- matrix(stats::rnorm(3 * p), 3, p); Ay <- matrix(stats::rnorm(3 * q), 3, q)
    X <- z %*% Ax + matrix(stats::rnorm(n*p),n,p)
    Y <- rho * z %*% Ay + sqrt(max(1-rho^2,.01)) * matrix(stats::rnorm(n*q),n,q)
    if (scenario == "weak") Y <- 0.2*z%*%Ay + matrix(stats::rnorm(n*q),n,q)
    if (scenario == "null") Y <- matrix(stats::rnorm(n*q),n,q)
    if (scenario == "collinear") {
      if (p >= 3) X[,3] <- X[,1] + .01*stats::rnorm(n)
      if (q >= 3) Y[,3] <- Y[,1] - .01*stats::rnorm(n)
    }
    if (scenario == "sparse") {
      X <- matrix(stats::rnorm(n*p),n,p); Y <- matrix(stats::rnorm(n*q),n,q)
      latent <- stats::rnorm(n); X[,seq_len(min(2,p))] <- X[,seq_len(min(2,p)),drop=FALSE] + latent
      Y[,seq_len(min(2,q))] <- Y[,seq_len(min(2,q)),drop=FALSE] + rho*latent
    }
    if (scenario == "outliers") {
      idx <- seq_len(max(3, floor(.05*n))); X[idx,] <- X[idx,] + 8; Y[idx,] <- Y[idx,] - 6
    }
    if (scenario == "heavy_tail") {
      X <- z %*% Ax + matrix(stats::rt(n*p,df=3),n,p); Y <- rho*z%*%Ay + matrix(stats::rt(n*q,df=3),n,q)
    }
    if (scenario == "monotone") {
      X <- matrix(stats::rnorm(n*p),n,p); latent <- X[,1]
      Y <- matrix(stats::rnorm(n*q,sd=.5),n,q); Y[,1] <- exp(latent/2)+stats::rnorm(n,sd=.2)
    }
    if (scenario == "nonmonotone") {
      X <- matrix(stats::runif(n*p,-2,2),n,p); latent <- X[,1]
      Y <- matrix(stats::rnorm(n*q,sd=.5),n,q); Y[,1] <- latent^2 + stats::rnorm(n,sd=.15)
    }
    if (scenario == "binary") { X <- 1*(X>0); Y <- 1*(Y>0) }
    if (scenario == "zero_inflated") {
      X[stats::runif(length(X))<.35] <- 0; Y[stats::runif(length(Y))<.35] <- 0
    }
    if (scenario == "missing") {
      X[sample(length(X),floor(.05*length(X)))] <- NA; Y[sample(length(Y),floor(.05*length(Y)))] <- NA
    }
    if (scenario == "functional") {
      tt <- seq(0,1,length.out=40); a <- stats::rnorm(n); b <- rho*a+sqrt(1-rho^2)*stats::rnorm(n)
      XF <- outer(a,sin(2*pi*tt)) + matrix(stats::rnorm(n*length(tt),sd=.2),n)
      YF <- outer(b,cos(2*pi*tt)) + matrix(stats::rnorm(n*length(tt),sd=.2),n)
      return(list(X=XF,Y=YF,time_x=tt,time_y=tt,truth=list(rho=rho),scenario=scenario))
    }
    colnames(X) <- paste0("X",seq_len(ncol(X))); colnames(Y) <- paste0("Y",seq_len(ncol(Y)))
    if (scenario == "three_blocks") {
      Z3 <- z %*% matrix(stats::rnorm(3*4),3,4)+matrix(stats::rnorm(n*4),n,4)
      colnames(Z3) <- paste0("Z",1:4)
      return(list(blocks=list(X=X,Y=Y,Z=Z3),truth=list(rho=rho),scenario=scenario))
    }
    list(X=X,Y=Y,truth=list(rho=rho),scenario=scenario)
  })
}

#' Synthetic agronomy soil-by-plant teaching data
#' @export
cca_agri_data <- function(n = 120, seed = 1) {
  .cca_with_seed(seed, {
    fertility <- stats::rnorm(n); water <- stats::rnorm(n); stress <- stats::rnorm(n)
    X <- cbind(
      soil_N = 25 + 5*fertility + stats::rnorm(n,sd=2),
      soil_P = 12 + 3*fertility + stats::rnorm(n,sd=1.5),
      soil_K = 80 + 10*fertility + stats::rnorm(n,sd=5),
      soil_moisture = 24 + 4*water + stats::rnorm(n,sd=2),
      soil_EC = 1.2 + .25*stress - .10*water + stats::rnorm(n,sd=.12),
      soil_pH = 6.1 + .25*fertility + stats::rnorm(n,sd=.15)
    )
    Y <- cbind(
      photosynthesis = 18 + 2.0*fertility + 1.4*water - 1.8*stress + stats::rnorm(n,sd=1.3),
      stomatal_conductance = .30 + .04*water - .035*stress + stats::rnorm(n,sd=.025),
      chlorophyll = 42 + 3.5*fertility - 1.2*stress + stats::rnorm(n,sd=2.2),
      biomass = 55 + 6*fertility + 4*water - 4*stress + stats::rnorm(n,sd=4),
      yield = 4.5 + .7*fertility + .45*water - .5*stress + stats::rnorm(n,sd=.35)
    )
    treatment <- factor(rep(LETTERS[1:4], length.out=n))
    list(X=as.data.frame(X),Y=as.data.frame(Y),treatment=treatment,
         description="Synthetic agronomy data: soil properties versus plant physiology and yield.")
  })
}

.cca_basis <- function(object, block="x", n_axes=2) {
  block <- match.arg(block,c("x","y"))
  W <- if(block=="x") object$weights_x else object$weights_y
  p <- if(block=="x") object$p else object$q
  if (is.null(W) || nrow(W)!=p) {
    L <- if(block=="x") object$loadings_x else object$loadings_y
    W <- L
  }
  n_axes <- min(n_axes,ncol(W),nrow(W))
  Q <- qr.Q(qr(W[,seq_len(n_axes),drop=FALSE]))[,seq_len(n_axes),drop=FALSE]
  rownames(Q) <- rownames(W); colnames(Q) <- paste0("Axis",seq_len(n_axes)); Q
}

.cca_start_basis <- function(data,n_axes=2) {
  Z <- scale(data); sv <- svd(Z,nu=0,nv=min(n_axes,ncol(Z)))
  sv$v[,seq_len(min(n_axes,ncol(sv$v))),drop=FALSE]
}

#' Interpolate between two projection bases
#' @export
cca_basis_tour <- function(data, start, target, frames = 40) {
  data <- .cca_as_matrix(data,"data")
  start <- as.matrix(start); target <- as.matrix(target)
  if (!all(dim(start)==dim(target)) || nrow(start)!=ncol(data)) .cca_stop("start/target basis dimensions must match the data variables.")
  ts <- seq(0,1,length.out=frames); bases <- vector("list",frames); projections <- vector("list",frames)
  for(i in seq_along(ts)) {
    M <- (1-ts[i])*start + ts[i]*target
    Q <- qr.Q(qr(M))[,seq_len(ncol(M)),drop=FALSE]
    bases[[i]] <- Q; projections[[i]] <- data %*% Q
  }
  out <- list(data=data,start=start,target=target,t=ts,bases=bases,projections=projections)
  class(out) <- "cca_tour"
  out
}

#' Canonical projection tour
#' @export
cca_tour <- function(object, block = c("x","y"), frames = 40,
                     start = c("pca","identity"), renderer = c("none","spinifex"),
                     animate = FALSE) {
  if(!inherits(object,"cca_fit")) .cca_stop("object must inherit from 'cca_fit'.")
  block <- match.arg(block); start <- match.arg(start); renderer <- match.arg(renderer)
  data <- if(block=="x") object$Xs else object$Ys
  target <- .cca_basis(object,block,min(2,ncol(data)))
  st <- if(start=="pca") .cca_start_basis(data,ncol(target)) else diag(ncol(data))[,seq_len(ncol(target)),drop=FALSE]
  out <- cca_basis_tour(data,st,target,frames); out$fit <- object; out$block <- block; out$type <- "canonical"
  if(renderer=="spinifex") {
    .cca_require("spinifex","cca_tour(renderer='spinifex')")
    out$rendered <- cca_tour_render(out, engine = "spinifex", animate = if (animate) "plotly" else "none")
  }
  out$animate <- animate
  out
}

#' Canonical-association guided tour
#' @export
cca_guided_tour <- function(object, block = c("x","y"), frames = 40, start = "pca") {
  block <- match.arg(block); tr <- cca_tour(object,block,frames,start,"none",FALSE)
  opp <- if(block=="x") object$scores_y[,1] else object$scores_x[,1]
  tr$objective <- vapply(tr$projections,function(P) abs(stats::cor(P[,1],opp,use="complete.obs")),numeric(1))
  tr$type <- "guided-canonical-association"; class(tr) <- c("cca_guided_tour","cca_tour"); tr
}

#' Grand-tour wrapper through tourr
#' @export
cca_grand_tour <- function(object, block = c("x","y"), ...) {
  .cca_require("tourr","cca_grand_tour")
  block <- match.arg(block); data <- if(block=="x") object$Xs else object$Ys
  structure(list(data=data,tour_path=tourr::grand_tour(),block=block,fit=object,args=list(...)),class="cca_grand_tour")
}

#' Compare two fitted canonical bases through a projection path
#' @export
cca_compare_tour <- function(object1, object2, block = c("x","y"), frames = 40) {
  block <- match.arg(block)
  data <- if(block=="x") object1$Xs else object1$Ys
  b1 <- .cca_basis(object1,block,min(2,ncol(data))); b2 <- .cca_basis(object2,block,ncol(b1))
  if(nrow(b1)!=nrow(b2)) .cca_stop("Models do not share a compatible variable space.")
  tr <- cca_basis_tour(data,b1,b2,frames); tr$fit1<-object1;tr$fit2<-object2;tr$block<-block;tr$type<-"model-comparison"
  tr
}

#' Classical-versus-resistant outlier tour
#' @export
cca_outlier_tour <- function(X, Y, block = c("x","y"), frames = 40, robust_engine = "winsor") {
  cmp <- cca_robust_compare(X,Y,robust_engine=robust_engine)
  tr <- cca_compare_tour(cmp$classical,cmp$resistant,match.arg(block),frames)
  tr$type <- "classical-to-resistant"; tr$comparison <- cmp; tr
}

#' Select static keyframes from a tour
#' @export
cca_tour_keyframes <- function(tour, n = 6, strategy = c("even","objective"),
                               plots = FALSE, group = NULL) {
  if(!inherits(tour,"cca_tour")) .cca_stop("tour must inherit from 'cca_tour'.")
  strategy <- match.arg(strategy); N <- length(tour$bases); n <- min(n,N)
  if(strategy=="objective" && !is.null(tour$objective)) {
    topn <- max(0L, n - 2L)
    extra <- if (topn) order(tour$objective,decreasing=TRUE)[seq_len(min(topn, length(tour$objective)))] else integer()
    idx <- sort(unique(c(1L,N,extra)))
    if (length(idx) > n) idx <- idx[seq_len(n)]
  } else idx <- unique(round(seq(1,N,length.out=n)))
  out <- list(index=idx,t=tour$t[idx],bases=tour$bases[idx],projections=tour$projections[idx],
              objective=if(!is.null(tour$objective)) tour$objective[idx] else NULL)
  if (isTRUE(plots)) out$plots <- lapply(idx, function(i) cca_tour_render(tour, "static", frame = i, group = group))
  class(out) <- "cca_tour_keyframes"
  out
}

#' Four-panel canonical association explorer
#' @export
cca_explorer <- function(object, component = 1, top = 15, interactive = FALSE) {
  net <- tryCatch(cca_plot(object,"network",threshold=.3,interactive=interactive),
                  error=function(e) cca_plot(object,"network",threshold=.15,interactive=interactive))
  out <- list(
    scores=cca_plot(object,"scores",component=component,interactive=interactive),
    loadings_x=cca_plot(object,"loadings",component=component,block="x",top=top,interactive=interactive),
    loadings_y=cca_plot(object,"loadings",component=component,block="y",top=top,interactive=interactive),
    network=net
  )
  class(out)<-"cca_explorer"; out
}

#' Canonical Correlation Learning Lab
#' @export
cca_teach <- function(X = NULL, Y = NULL, level = c("intro","intermediate","advanced"), seed = 1) {
  level <- match.arg(level)
  if(is.null(X)||is.null(Y)) { d<-cca_agri_data(seed=seed); X<-d$X;Y<-d$Y }
  fit <- cca_classic(X,Y,ncomp=min(3,ncol(.cca_as_matrix(X)),ncol(.cca_as_matrix(Y))))
  content <- list(
    intro=c(
      "CCA searches for paired linear combinations U = Xa and V = Yb with maximum correlation.",
      "The first pair maximizes cross-block association; later classical pairs are constrained relative to earlier canonical variates.",
      "Canonical weights define the combinations, while structure coefficients describe correlations between original variables and canonical scores."
    ),
    intermediate=c(
      "Interpret canonical correlations together with structure coefficients, cross-loadings, variance extracted and redundancy.",
      "Large canonical correlation does not imply that either block captures much variance from the other block.",
      "Bootstrap stability and held-out association help distinguish stable signal from sample-specific optimization."
    ),
    advanced=c(
      "Regularization is essential when dimensions are large relative to sample size or covariance matrices are ill-conditioned.",
      "Sparse methods trade association strength for variable selection and require stability assessment.",
      "Resistant, copula, kernel, probabilistic, Bayesian, functional and multiblock variants address different departures from classical assumptions; model comparison should use appropriate resampling."
    )
  )
  chosen <- switch(level,intro=content$intro,intermediate=c(content$intro,content$intermediate),advanced=unlist(content,use.names=FALSE))
  equations <- c(
    objective="rho = (a' Sxy b) / sqrt[(a' Sxx a)(b' Syy b)]",
    redundancy="Redundancy_X,k = mean(loadings_X,k^2) * rho_k^2",
    generalization="Compare rho_train with rho_validation under a resampling design compatible with the experimental unit."
  )
  out <- list(level=level,narrative=chosen,equations=equations,fit=fit,
              tables=list(correlations=cca_table(fit,"correlations"),loadings=cca_table(fit,"loadings"),redundancy=cca_table(fit,"redundancy")),
              plots=list(scree=cca_plot(fit,"scree"),scores=cca_plot(fit,"scores"),circle=cca_plot(fit,"circle")),
              tour=cca_guided_tour(fit,"x",frames=30))
  class(out)<-"cca_teaching_lab"; out
}

#' @export
print.cca_tour <- function(x,...) {cat("CCAflow tour:",x$type %||% "basis path","\n  Frames:",length(x$bases),"\n");invisible(x)}
#' @export
print.cca_teaching_lab <- function(x,...) {cat("CCAflow teaching lab -",x$level,"\n\n",paste0("- ",x$narrative,collapse="\n"),"\n");invisible(x)}

.cca_tour_basis_array <- function(tour) {
  a <- simplify2array(tour$bases)
  if (length(dim(a)) == 2L) a <- array(a, dim = c(nrow(a), ncol(a), 1L))
  a
}

#' Render a projection tour or a static publication keyframe
#' @export
cca_tour_render <- function(tour, engine = c("spinifex", "static"), frame = NULL,
                            animate = c("none", "plotly", "gganimate"),
                            group = NULL, point_size = 1.7, alpha = 0.75, ...) {
  if (!inherits(tour, "cca_tour")) .cca_stop("tour must inherit from 'cca_tour'.")
  engine <- match.arg(engine); animate <- match.arg(animate)
  if (!is.null(group) && length(group) != nrow(tour$data)) .cca_stop("group must have one value per observation.")
  if (engine == "static") {
    .cca_require("ggplot2", "cca_tour_render(engine='static')")
    if (is.null(frame)) frame <- length(tour$projections)
    frame <- max(1L, min(as.integer(frame), length(tour$projections)))
    P <- as.data.frame(tour$projections[[frame]])
    if (ncol(P) == 1L) P$V2 <- 0
    names(P)[1:2] <- c("Axis1", "Axis2")
    P$group <- if (is.null(group)) factor("Observations") else factor(group)
    p <- ggplot2::ggplot(P, ggplot2::aes(x = Axis1, y = Axis2, shape = group)) +
      ggplot2::geom_point(size = point_size, alpha = alpha) +
      ggplot2::coord_equal() +
      ggplot2::labs(x = "Projection axis 1", y = "Projection axis 2",
                    title = sprintf("%s tour, frame %d/%d", tour$type %||% "Projection", frame, length(tour$projections)),
                    subtitle = sprintf("Path position = %.3f", tour$t[frame])) +
      theme_cca_publication()
    if (is.null(group)) p <- p + ggplot2::guides(shape = "none")
    return(p)
  }
  .cca_require("spinifex", "cca_tour_render(engine='spinifex')")
  arr <- .cca_tour_basis_array(tour)
  # The tour path is already interpolated by cca_basis_tour(). Passing angle = 0
  # makes spinifex use the supplied frames directly, which also avoids a known
  # tourr (>= 1.2.8) interpolate() failure with R >= 4.5 on pre-computed arrays.
  dots <- list(...)
  if (!"angle" %in% names(dots)) dots$angle <- 0
  g <- do.call(spinifex::ggtour, c(list(basis_array = arr, data = tour$data), dots)) +
    spinifex::proto_default(aes_args = if (is.null(group)) list(alpha = alpha) else list(shape = factor(group), alpha = alpha)) +
    spinifex::proto_basis()
  if (animate == "plotly") return(spinifex::animate_plotly(g))
  if (animate == "gganimate") {
    .cca_require("gganimate", "cca_tour_render(animate='gganimate')")
    return(spinifex::animate_gganimate(g))
  }
  g
}

#' @export
plot.cca_tour <- function(x, frame = NULL, group = NULL, ...) {
  cca_tour_render(x, engine = "static", frame = frame, group = group, ...)
}
