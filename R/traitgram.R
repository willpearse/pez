##' Traitgram for comparative community object
##'
##' A wrapper for the \code{\link{traitgram}} function in the
##' \code{picante} package.
##'
##' @param object A \code{\link{comparative.comm}} object.
##' @param trait Which trait to plot.  If \code{\link{missing}}, use
##' the first trait.  If a positive \code{\link{numeric}} vector of
##' \code{\link{length}} one, use the \code{as.integer(trait)}th trait.  If a
##' \code{\link{numeric}} vector, use it instead of the trait data in
##' \code{object}.  If a \code{\link{character}} vector of
##' \code{\link{length}} one, use the trait with that name.  If a
##' \code{\link{function}} pass the trait data frame through that
##' function and use the result (see \code{\link{princompOne}} for a
##' useful function).  If an \code{\link{expression}}, evaluate that
##' expression in the environment of the trait data and use the result.  If a
##' \code{\link{character}} vector, then convert to an expression and
##' evaluate in the environment of the trait data and use the result.
##' @param moreArgs List of more arguments to pass on to \code{trait}
##' (if its a \code{\link{function}}).
##' @param ... Additional arguments to be passed on to
##' \code{\link{traitgram}}.
##' @return See \code{\link{traitgram}}
##' @importFrom picante traitgram
##' @export
traitgram.cc <- function(object, trait, moreArgs = NULL, ...) {
    if(is.null(object$data)) stop("must supply trait information")
    if(missing(trait)) {
        if(is.null(dim(object$data))) {
            tt <- object$data
        } else {
            tt <- object$data[, 1]
        }
    } else if(is.numeric(trait)) {
        if(length(trait) == 1) {
            if(trait < 1) stop("trait can't be a negative number")
            tt <- object$data[, as.integer(trait)]
        } else {
            tt <- trait
        }
    } else if(is.function(trait)) {
        tt <- do.call(trait, c(list(object$data), moreArgs))
    } else if(is.language(trait)) {
        tt <- with(object$data, eval(trait))
    } else if(is.character(trait)) {
        trait <- parse(text = paste(trait, collapse = "; "))
        tt <- with(object$data, eval(trait))
    }

                                        # FIXME: not sure this this naming is necessary
    if(is.null(rownames(object$data))) {
        names(tt) <- names(object$data)
    } else {
        names(tt) <- rownames(object$data)
    }

                                        # resolve possible polytomies
    pp <- multi2di(object$phy)

                                        # plot
    traitgram(tt, pp, ...)
}

##' First axis of a principal components analysis
##'
##' A wrapper for \code{\link{princomp}}
##' 
##' @param x A matrix-like object
##' @param ... Arguments to pass on to \code{\link{princomp}}
##' @return The first axis of a PCA
##' @export
princompOne <- function(x, ...) princomp(x, ...)$scores[,1]
