#' Phylogenetically weighted regression: a method for modeling
#' non-stationarity on evolutionary trees
#'
#' \code{pwr.cv} performs model cross validation on the holdout set
#' \code{holdout}, K-fold cross validation can be specified using
#' \code{fold}. \code{pwr.treePlot} plots the phylogeny and allows for
#' plotting of PWR model coefficients at the tips of the tree. PWR
#' model coefficients to be pllotted must be appended to the phylo4d
#' object. \code{pwr.phylobubbles} plots either circles or squares
#' corresponding to the magnitude of each cell of a phylo4d object,
#' and allows the plotting of PWR moel coefficients.  Model
#' coefficients to be plotted must be appended to the phylo4d object
#' @param formula Model formula to be fit
#' @param phy4d \code{phylobase} format object to be used for fitting
#' @param bw Bandwidth for model fit
#' @param wfun Weighting function to be used; myust be one of
#'     \code{brownian}, \code{martins}, \code{gaussian}, or
#'     \code{Gaus}
#' @param verbose Whether to be verbose during model fitting (default
#'     \code{FALSE})
#' @param formula Model formula to be fit
#' @param phy4d \code{phylobase} format object to be used for fitting
#' @param holdout vector of integers indicating position of tips on
#'     the tree to be held out during model fitting
#' @param coef coefficients to extract from the model; must be one of
#'     \code{formula}, \code{slope}, or \code{all}
#' @param model Model type to be fit; must be either \code{pwr}, or
#'     \code{pgls}
#' @param which.tips vector of integers indicating position of tips on
#'     the tree to be used, defaults to all species on the tree
#' @param phy A phylo4 or phylo4d object
#' @note Some information about this
#' @return Phylogenetically weighted model object; use \code{pwr.coef}
#'     to extract coefficients, and \code{pwr.tp} and
#'     \code{pwr.phylobubbles} for plotting.
#' @author Jim Regetz, T. Jonathan Davies, Elizabeth M. Wolkovich,
#'     Brian J. McGill, and a tiny bit of cosmetic tidying from Will
#'     Pearse. Plotting code (\code{pwr.tp} and
#'     \code{pwr.phylobubbles}) only slightly adapted from code in
#'     \code{phylobase}; thank you to the authors of \code{phylobase}!
#' 
#' @seealso caper:pgls phylobase:tp phylobase:phylobubbles
#' @references Davies, T. J., Regetz, J., Wolkovich, E. M., & McGill,
#'     B. J. (2019). Phylogenetically weighted regression: A method
#'     for modelling non-stationarity on evolutionary trees. Global
#'     Ecology and Biogeography, 28(2), 275-285.
#' @examples
#' # Below are some examples in order of increasing complexity: (1) a
#' # simple demonstration, and (2) some cross-validation examples.
#' 
#' ##########################################
#' # (1) A simple example ###################
#' ##########################################
#' # Make up some data
#' # - Note the deep bifurcation in the slope term
#' tree <- read.tree(text=
#'     "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);")
#' dat <- data.frame(
#'    x = c(0, 1, 2, 3, 0, 1, 2, 3),
#'    y = rnorm(8, c(1, 2, 3, 4, 2, 1.5, 1, 0.5), sd=0.25)
#' )
#'
#' # Combine in phylobase format, and then run a scatterplot
#' c.data <- phylo4d(tree, dat)
#' plot(tipData(c.data)$x, tipData(c.data)$y, pch=rep(c(1,16), each=4),
#'     xlab="Predictor", ylab="Response", bty="l")
#' text(tipData(c.data)$x, tipData(c.data)$y, labels=rownames(tipData(c.data)),
#'     pos=3, cex=0.8)
#'
#' # Phylogenetically weighted-regression with various weighting functions
#' pwr.b <- pwr.coef(pwr.fit(y ~ x, combined.data, wfun="brownian"))
#' pwr.m <- pwr.coef(pwr.fit(y ~ x, combined.data, wfun="martins"))
#' pwr.g <- pwr.coef(pwr.fit(y ~ x, combined.data, wfun="gaussian"))
#' pwr.G <- pwr.coef(pwr.fit(y ~ x, combined.data, wfun="Gauss"))
#' #...note use of 'pwr.coef' to get estimates from the models
#' 
#' # Run a PGLS
#' pgls <- gls(y ~ x, data=dat, correlation=corBrownian(phy=tree))
#'
#' # Contrast the PWR and PGLS estimates
#' pwr.treePlot(addData(c.data, data.frame(gest=coef(pgls)[2], glb=confint(pgls)[2,1],
#'    gub=confint(pgls)[2,2], pwr.m))[,c("est", "lb", "ub", "gest", "glb",
#'    "gub")], show.tip.label=TRUE, lower=-1., upper=1.5, pex=2, aex=3.5)
#' @importFrom ape read.tree
#' @importFrom nlme gls corMatrix
#' @importFrom phylobase phylo4d
#' @importMethodsFrom phylobase phylo4d
#' @importFrom subplex subplex
#' @importFrom grid grid.newpage
#' @rdname pwr
#' @name pwr
#' @aliases pwr pwr.fit pwr.coef pwr.treePlot pwr.phylobubbles
#######################################
# Fitting functions ###################
#######################################
#' @export
pwr.fit <- function(formula, phy4d, which.tips, bw, wfun, verbose=FALSE) {
    if (missing(which.tips)) {
        which.tips <- seq.int(nTips(phy4d))
    }
    if (missing(bw)) {
        if (verbose) cat("Computing optimal bandwidth...")
        bw <- .get.opt.bw(formula, phy4d, wfun)
        if (verbose) cat("done.\n")
    }
    if (verbose) cat("Using bandwidth of", bw, "\n",
        "Generating weights matrix\n")
    wts <- .getWts(phy4d, bw, wfun)
    if (verbose) cat("Running PWR...\n")
    if (verbose) pb <- txtProgressBar(min=0, max=nTips(phy4d), style=3)
    ans <- lapply(which.tips, function(i) {
        if (verbose) setTxtProgressBar(pb, i)
        .pwr.wts(formula, phy4d, wts[,i])
    })
    if (verbose) close(pb)
    attr(ans, "weights") <- wts
    class(ans) <- "phy.wt.reg"
    return(ans)
}
#' @export
pwr.coef <- function(pwrres) {
    if(!inherits(pwrres, "phy.wt.reg"))
        stop("Error, ", pwrres, " is not a phylogenetically weighted regression")
    est <- sapply(seq_along(pwrres), function(i) {
        if (nrow(pwrres[[i]])==2) {
            pwrres[[i]][2,1]
        } else NA
    })
    lb <- sapply(seq_along(pwrres), function(i) {
        if (nrow(pwrres[[i]])==2) {
            pwrres[[i]][2,1] - 1.96*pwrres[[i]][2,2]
        } else NA
    })
    ub <- sapply(seq_along(pwrres), function(i) {
        if (nrow(pwrres[[i]])==2) {
            pwrres[[i]][2,1] + 1.96*pwrres[[i]][2,2]
        } else NA
    })
    data.frame(est, lb, ub)
}
#' @export
# Cross-validation
pwr.cv <-  function(formula, phy4d, holdout, coef=c("formula","slope","all"), model=c("pwr","pgls")) {
    coef <- match.arg(coef)
    model <- match.arg(model)

    if(model=="pwr"){
        # extract training set from phy4d
        phy.train <- phy4d[setdiff(seq(nTips(phy4d)), holdout)]
        if (missing(bwidth)) {
            bwidth <- get.opt.bw(formula, phy.train, wfun=wfun, method=method)
        }
        # get weights vectors for *all* species, but only keep rows
        # corresponding to training set
        wts.train <- getWts(phy4d, bwidth, wfun)[-holdout,]
        # extract training data from phy4d
        dat.train <- tipData(phy.train)
        # loop over each point and do a weighted least squares regression
        # with weights based on distance
        yhat <- sapply(holdout, function(p) {
            w <- wts.train[,p]
            b <- lm(formula, data=data.frame(dat.train, w), weights=w)
            predict(b, tipData(phy4d)[p,])
        })
        if(coef=="formula")
            return(data.frame(
                y=tipData(phy4d)[holdout, as.character(formula)[[2]]],
                yhat.pwr=yhat)
                )
        if(coef=="slope")
            return(data.frame(
                slope=tipData(phy4d)[holdout, "true_slope"],
                slope.pwr=yhat)
                )
        if(coef=="all")
            return(data.frame(
                slope=tipData(phy4d)[holdout,],
                slope.pwr=yhat)
                )
        } else  {
            # extract training set from phy4d
            p4.train <- phy4d[setdiff(seq(nTips(phy4d)), holdout)]
            phy.train <- suppressWarnings(as(p4.train, "phylo"))
            # extract training data from phy4d
            dat.train <- tipData(p4.train)
            # loop over each point and do a weighted least squares regression
            # with weights based on distance
            pgls <- do.call(gls, list(model=formula, data=dat.train,
                                      correlation=corBrownian(phy=phy.train)))
            if(coef=="all")
                return(data.frame(
                    y=tipData(phy4d)[holdout, ],
                    yhat.pgls=predict(pgls, tipData(phy4d)[holdout,]))
                    )
            if(coef=="formula")
                return(data.frame(
                    y=tipData(phy4d)[holdout, as.character(formula)[[2]]],
                    yhat.pgls=predict(pgls, tipData(phy4d)[holdout,]))
                    )
            if(coef=="slope")
                return(data.frame(
                    slope=tipData(phy4d)[holdout, "true_slope"],
                    slope.pgls=pgls$coefficients[[2]])
                    )
        }
}

#######################################
# Plotting functions ##################
#######################################
#' @export
#' @importFrom grid grid.newpage
# modified version of phylobase:::treePlot
pwr.treePlot <- function (phy, show.tip.label=TRUE, show.node.label=FALSE,
    tip.order=NULL, tip.plot.fun="bubbles", plot.at.tip=TRUE,
    edge.color="black", node.color="black", tip.color="black",
    edge.width=1, newpage=TRUE, margins=c(1.1, 1.1, 1.1, 1.1), ...) {
    if (!inherits(phy, "phylo4"))
        stop("treePlot requires a phylo4 or phylo4d object")
    if (!isRooted(phy))
        stop("treePlot function requires a rooted tree.")
    if (newpage)
        grid.newpage()
    type <- "phylogram"
    Nedges <- nEdges(phy)
    Ntips <- nTips(phy)
    if (!is.null(tip.order) && length(tip.order) > 1) {
        if (length(tip.order) != Ntips) {
            stop("tip.order must be the same length as nTips(phy)")
        }
        if (is.numeric(tip.order)) {
            tip.order <- tip.order
        } else {
            if (is.character(tip.order)) {
                tip.order <- as.numeric(names(tipLabels(phy))[
                    match(tip.order, tipLabels(phy))])
            }
        }
        tip.order <- rev(tip.order)
    }
    if (!hasEdgeLength(phy) || type == "cladogram") {
        edgeLength(phy) <- rep(1, Nedges)
    }
    xxyy <- phyloXXYY(phy, tip.order)
    pushViewport(plotViewport(margins=margins))
    pwr.phylobubbles(type=type, show.node.label=show.node.label,
        rot=0, edge.color=edge.color, node.color=node.color,
        tip.color=tip.color, edge.width=edge.width,
        show.tip.label=show.tip.label, newpage=TRUE, ..., XXYY=xxyy)
    upViewport()
}
#' @export
pwr.phylobubbles <- function (type=type, place.tip.label="right",
    show.node.label=show.node.label, show.tip.label=show.tip.label,
    edge.color=edge.color, node.color=node.color, tip.color=tip.color,
    edge.width=edge.width, newpage=TRUE, cex=1, pex=1, aex=1, ..., XXYY,
    square=FALSE, show.estimates=TRUE, lower=NULL, upper=NULL) {

    nVars <- 1
    lab.right <- ifelse(place.tip.label %in% c("right", "both"),
        TRUE, FALSE) && show.tip.label
    lab.left <- ifelse(place.tip.label %in% c("left", "both"),
        TRUE, FALSE) && show.tip.label
    phy <- XXYY$phy
    tmin <- min(tdata(phy, type="tip"), na.rm=TRUE)
    tmax <- max(tdata(phy, type="tip"), na.rm=TRUE)
    pedges <- edges(phy)
    tip.order <- XXYY$torder
    tipdata <- tdata(phy, type="tip")[tip.order, , drop=FALSE]
    dlabwdth <- max(stringWidth(colnames(tipdata))) * 1.2
    if (convertWidth(dlabwdth, "cm", valueOnly=TRUE) < 2) {
        dlabwdth <- unit(2, "cm")
    }
    phyplotlayout <- grid.layout(nrow=2, ncol=2, heights=unit.c(unit(1,
        "null"), dlabwdth), widths=unit(c(1, 1), c("null",
        "null"), list(NULL, NULL)))
    pushViewport(viewport(layout=phyplotlayout, name="phyplotlayout"))
    pushViewport(viewport(layout.pos.row=1:2, layout.pos.col=2,
        height=unit(1, "npc") + convertUnit(dlabwdth, "npc"),
        name="bubbleplots", default.units="native"))
    tys <- XXYY$yy[pedges[, 2] <= nTips(phy)]
    tys <- tys[match(names(tipLabels(phy))[tip.order], XXYY$torder)]
    maxr <- ifelse(ncol(tipdata) > nTips(phy), 1/ncol(tipdata),
        1/nTips(phy))
    tipdataS <- apply(tipdata, 2, function(x) (maxr * x)/max(abs(x),
        na.rm=TRUE))
    if (is.null(lower)) {
        lower <- min(tipdata)
    }
    if (is.null(upper)) {
        upper <- max(tipdata)
    }
    tipdataS2 <- .rescale(tipdata, lower, upper)
    if (nVars == 1) {
        xpos <- 0.5
    } else {
        xpos <- seq(0 + maxr + 0.12, 1 - maxr - 0.12, length.out=nVars)
    }
    xrep <- rep(xpos, each=length(tys))
    yrep <- rep(tys, nVars)
    ccol <- ifelse(tipdata < 0, "black", "white")
    naxs <- matrix(xrep, ncol=nVars)
    nays <- matrix(yrep, ncol=nVars)
    dnas <- is.na(tipdataS)
    naxs <- naxs[dnas]
    nays <- nays[dnas]
    tipdataS[is.na(tipdataS)] <- 0 + 0.001
    if (lab.right) {
        tiplabwidth <- max(stringWidth(tipLabels(phy)))
    } else {
        tiplabwidth <- unit(0, "null", NULL)
    }
    bublayout <- grid.layout(nrow=2, ncol=2, widths=unit.c(unit(1,
        "null", NULL), tiplabwidth), heights=unit.c(unit(1,
        "null", NULL), dlabwdth))
    pushViewport(viewport(x=0.5, y=0.5, width=0.95, height=1,
        layout=bublayout, name="bublayout"))
    pushViewport(viewport(name="bubble_plots", layout=bublayout,
        layout.pos.col=1, layout.pos.row=1))

    # plot x-axis
    labs <- pretty(c(lower, upper))
    vals <- .rescale(labs, lower, upper)
    labs <- format(labs[0 <= vals & vals <= 1], nsmall=1)
    vals <- vals[0 <= vals & vals <= 1]
    ex <- -0.02 * aex
    grid.segments(x0=min(vals), x1=max(vals), y0=ex, y1=ex)
    grid.segments(x0=vals, x1=vals, y0=ex+0.01, y1=ex,
        gp=gpar(col="black"))
    grid.text(labs, x=vals, y=ex-0.01*aex, gp=gpar(cex=0.5*cex))
    grid.text("Coefficient", x=0.5, y=ex-0.02*aex, gp =
        gpar(cex=0.5*cex))

    # plot interesting results
    grid.segments(x0=tipdataS2$gest, x1=tipdataS2$gest, y0=0,
        y1=1, gp=gpar(col="grey"))
    grid.segments(x0=tipdataS2$glb, x1=tipdataS2$glb, y0=0, y1=1,
        gp=gpar(col="grey", lty="dashed"))
    grid.segments(x0=tipdataS2$gub, x1=tipdataS2$gub, y0=0, y1=1,
        gp=gpar(col="grey", lty="dashed"))
    grid.segments(x0=tipdataS2$lb, x1=tipdataS2$ub, y0=tys, y1=tys,
        gp=gpar(col="red"))
    if (!is.null(tipdataS2$lb.1) & !is.null(tipdataS2$ub.1)) {
        grid.segments(x0=tipdataS2$lb.1, x1=tipdataS2$ub.1,
            y0=tys+0.3*diff(tys[1:2]), y1=tys+0.3*diff(tys[1:2]),
            gp=gpar(col="grey"))
    }
    if (!is.null(tipdataS2$lb.2) & !is.null(tipdataS2$ub.2)) {
        grid.segments(x0=tipdataS2$lb.2, x1=tipdataS2$ub.2,
            y0=tys+0.6*diff(tys[1:2]), y1=tys+0.6*diff(tys[1:2]),
            gp=gpar(col="green"))
    }
    if (!is.null(tipdataS2$simcoef)) {
        grid.points(tipdataS2$simcoef, yrep, pch=1, size =
            unit(0.03*pex, "npc"))

    }
    if (show.estimates) {
        grid.points(tipdataS2$est, yrep, pch=16, size=unit(0.02*pex,
            "npc"), gp=gpar(col="red"))

    }

    if (length(naxs) > 0) {
        grid.points(naxs, nays, pch=4)
    }
    upViewport()
    if (lab.right) {
        pushViewport(viewport(name="bubble_tip_labels", layout=bublayout,
            layout.pos.col=2, layout.pos.row=1))
        tt <- sub("_", " ", tipLabels(phy)[tip.order])
        grid.text(tt, 0.1, tys, just="left", gp=gpar(cex=0.5*cex))
        upViewport()
    }
    pushViewport(viewport(name="bubble_data_labels", layout=bublayout,
        layout.pos.col=1, layout.pos.row=2))
    datalaboffset <- convertUnit(unit(15, "mm"), "npc", valueOnly=TRUE)
    upViewport(3)
    pushViewport(viewport(layout.pos.row=2, layout.pos.col=1,
        name="bubblelegend"))
    yyy <- phylobase:::.bubLegendGrob(tipdata, tipdataS)
    grid.draw(yyy)
    upViewport()
    pushViewport(viewport(layout.pos.row=1, layout.pos.col=1,
        name="tree"))
    plotOneTree(XXYY, "phylogram", show.tip.label=FALSE, show.node.label,
        edge.color, node.color, tip.color, edge.width, rot=0)
    upViewport(2)

}

#######################################
# Internal functions ##################
#######################################
# Run PWR given weight function and bandwidth
.pwr.wfn <- function(formula, phy4d, wfun, bwidth, holdout=FALSE) {
    wts <- .getWts(phy4d, bwidth, wfun)
    bs <- list()
    res <- c()
    yhat <- c()
    # loop over each point and do a weighted least squares regression
    # with weights based on distance
    for (p in 1:nTips(phy4d)) {
        w <- wts[,p]
        if (holdout) {
            w[p] <- 0
        }
        b <- lm(formula, data=data.frame(tipData(phy4d), w), weights=w)
        bs[[p]] <- coef(b)
        res[p] <- resid(b)[p]
        yhat[p] <- fitted(b)[p]
    }
    coef <- do.call("rbind", bs)
    return(list(coef=coef, resid=res, fitted=yhat))
}
# Run PWR given explicit weight values
.pwr.wts <- function(formula, phy4d, wts, verbose=FALSE, plot=FALSE) {
    tree <- suppressWarnings(as(phy4d, "phylo"))
    dat <- data.frame(tipData(phy4d), wts=wts)
    # do straight pwr
    pwr.fit <- lm(formula, data=dat, weights=wts)
    if (verbose) {
        cat("\nPWR confidence intervals:\n")
        print(confint(pwr.fit))
    }
    return(coef(summary(pwr.fit)))
}
# optimal bandwidth finder
.get.opt.bw <- function(formula, phy4d, wfun, interval, method="subplex",
    verbose=FALSE) {

    if (wfun == "brownian") {
        if (verbose) {
            message("no bandwidth for brownian distance; returning NULL")
        }
        return(NULL)
    }
    # set bounds on possible bandwidth values - these work fine for MS
    # purposes but may need to be revisited for generality
    if (missing(interval)) {
        dist <- cophenetic(suppressWarnings(as(phy4d, "phylo")))
        if (wfun %in% c("gaussian", "Gauss")) {
            lo <- 0.01
            hi <- max(dist)/1.5
        } else if (wfun == "exponential") {
            lo <- 0.01
            hi <- 20
        } else if (wfun == "martins") {
            lo <- 0
            hi <- -log(1e-6)/min(dist[lower.tri(dist)])
        }
        interval <- c(lo, hi)
    }
    if (verbose) message("-- Running ", method, " --")
    if (method=="subplex") {
        optfn <- function(logbwidth) {
            sum(.pwr.wfn(formula, phy4d, wfun, exp(logbwidth), TRUE)$resid^2)
        }
        runtime <- system.time(res <- subplex(-1, optfn))
        if (res$convergence!=0) {
           warning(paste("bandwidth optimization problem (code ",
               res$convergence, ")", sep=""))
        }
        opt.bw <- exp(res$par)
    } else if (method=="optimize") {
        optfn <- function(bwidth) {
            sum(.pwr.wfn(formula, phy4d, wfun, bwidth, TRUE)$resid^2)
        }
        runtime <- system.time(res <- optimize(optfn, interval=interval))
        opt.bw <- res$minimum
    } else if (method=="L-BFGS-B") {
        optfn <- function(bwidth) {
            sum(.pwr.wfn(formula, phy4d, wfun, bwidth, TRUE)$resid^2)
            #vals <- .pwr.wfn(formula, phy4d, wfun, bwidth, TRUE)
            #if (any(is.na(vals$coef[, "x"]))) Inf else sum(vals$resid^2)
        }
        runtime <- system.time(res <- optim(1, optfn, lower=0,
            method="L-BFGS-B"))
        opt.bw <- res$par
    } else if (method=="nlm") {
        optfn <- function(bwidth) {
            sum(.pwr.wfn(formula, phy4d, wfun, bwidth, TRUE)$resid^2)
        }
        runtime <- system.time(res <- nlm(optfn, 0))
        opt.bw <- res$minimum
    } else {
        stop("invalid optimization algorithm")
    }
    if (verbose) {
        print(runtime)
        message("bandwidth = ", opt.bw)
    }
    opt.bw
}
# create full matrix of weights
.getWts <- function(phy4d, bw, wfun) {
    data <- tipData(phy4d)
    tree <- suppressWarnings(as(phy4d, "phylo"))
    dist <- cophenetic(tree)
    switch(wfun,
        gaussian = {
            sqd <- sqrt(dist)
            dnorm(t(t(sqd) / apply(sqd, 2, sd)) / bw)
        },
        exponential = {
            exp(-dist/bw)
        },
        Gauss = {
            exp((-0.5) * ((dist^2)/(bw^2)))
        },
        brownian = {
            corMatrix(Initialize(corBrownian(phy=tree), data))
        },
        martins = {
            corMatrix(Initialize(corMartins(bw, phy=tree), data))
        },
        stop("invalid distance weighting scheme")
    )
}
.rescale <- function(x, lower=NULL, upper=NULL, na.rm=TRUE) {
    if (is.null(lower)) {
        lower <- min(x, na.rm=na.rm)
    }
    if (is.null(upper)) {
        upper <- max(x, na.rm=na.rm)
    }
    (x - lower) / (upper - lower)
}
