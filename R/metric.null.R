#' Compare any metric(s) with null distributions
#'
#' Compare any calculated metric within \code{pez} to values expected
#' under some null distribution. This is a very light wrapper around
#' the utility functions \code{\link{.generic.null}} and
#' \code{\link{generic.metrics}} (which is, itself, a very simple
#' function!).
#'
#' @note \code{comp.fun} can be \emph{anything}; much ink has been
#' written about the use of standard effect sizes in eco-phylogenetic
#' analyses (\emph{e.g.}, Kembel 2009). That this function makes it
#' easy for you to compute Standard Effect Sizes does not necessarily
#' mean that you should (see Pearse et al. 2013).
#'
#' Calculating null permutations on a dispersion metric makes little
#' sense, since (by definition; see Pearse et al. 2014) a dispersion
#' metric \emph{require} the use of a null distribution to be
#' calculated. There is nothing to stop you doing so, however! The
#' code makes no attempt to stop you calculating null dissimilarity
#' metrics, but I am not certain that doing so is a good idea using
#' this code as I don't know what to do with such null models!
#' @param data data \code{\link{comparative.comm}} object
#' @param permute number of null permutations to perform (default
#' 1000)
#' @param metrics vector of functions to be calculated on \code{data};
#' see \link{pez.metrics} for a list of them.
#' @param null.model one of "taxa.labels", "richness", "frequency",
#' "sample.pool", "phylogeny.pool", "independentswap", or
#' "independentswap". These correspond to the null models available in
#' \code{\link{picante}}
#' @param permute number of null permutations
#' @param comp.fun comparison function to compare observed values with
#' null values. Default is \code{\link{.ses}}; this is a Standard
#' Effect Size (obs - mean)/SEmean. You may supply your own function;
#' it should take the observed site-metric matrix as its first
#' argument, and a site-metric-permutation array as its second. See
#' the internals of \code{\link{generic.null}} for an example of its
#' use.
#' @return Output from \code{comp.fun}, by default a vector of
#' standard effect sizes
#' @param ... additional arguments (e.g, \code{dist},
#' \code{abundance.weighted}) to be passed to any metric functions
#' (see \code{\link{generic.metrics}} for possible arguments)
#' @author Will Pearse
#' @references Kembel S.W. (2009) Disentangling niche and neutral
#' influences on community assembly: assessing the performance of
#' community phylogenetic structure tests. Ecology letters, 12(9),
#' 949-960.
#' @references Pearse W.D., Jones F.A. & Purvis A. (2013) Barro
#' Colorado Island's phylogenetic assemblage structure across fine
#' spatial scales and among clades of different ages. Ecology, 94(12),
#' 2861-2872.
#' @export
#' @rdname generic.metric
#' @name generic.metric
generic.null <- function(data, metrics, null.model=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), permute=1000, comp.fun=.ses, ...){
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative.comm object")
    if(permute < 0) stop("Can't have negative null permutations!")
    null.model <- match.arg(null.model)
    if(!is.list(metrics) | !all(sapply(metrics, is.function)))
        stop("'metrics' must be a list of functions")
    
    #Calculate real values
    observed <- generic.metrics(data, metrics, ...)
    
    #Calculate null distributions and structure
    null <- .generic.null(data, metrics, null.model, permute, ...)

    #Perform comparison and return
    return(comp.fun(observed, null))
}

#' Calculate Standard Effect Sizes of metrics
#' @param observed observed metric values in site-metric matrix
#' (\emph{e.g.}, from \link{\code{generic.metrics}})
#' @param null null distributions (\emph{e.g.}, from
#' \link{\code{.metric.null}}) in a site-metric-permutation array
#' @return Vector of standard effect sizes
#' @export
#' @rdname generic.metric
#' @name generic.metric
.ses <- function(observed, null){
    means <- apply(null, 1:2, mean)
    ses <- apply(null, 1:2, function(x) sd(x)/sqrt(length(x)))
    return((observed - means)/sds)
}

#' Produce null randomisations and compute metrics across them
#' @return site-metric-permutation array
#' @export
#' @rdname generic.metric
#' @name generic.metric
.metric.null <- function(data, metrics, null.model=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), permute=1000, ...){
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative.comm object")
    if(permute < 0) stop("Can't have negative null permutations!")
    null.model <- match.arg(null.model)
    if(!is.list(metrics) | !all(sapply(metrics, is.function)))
        stop("'metrics' must be a list of functions")

    #Calculate null distributions and structure
    null <- apply(replicate(permute, randomizeMatrix(x$comm, null.model=null.model)), 3, comparative.comm, phy=x$phy)
    null <- lapply(null, function(x) do.call(cbind, lapply(metrics, function(y) y(x))))
    null <- array(unlist(null), c(dim(null[[1]]), length(null)))

    return(null)
}

#' Calculate arbitrary metrics across data
#' @return site-metric matrix
#' @export
#' @rdname generic.metric
#' @name generic.metric
generic.metric <- function(data, metrics, ...){
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative.comm object")
    if(!is.list(metrics) | !all(sapply(metrics, is.function)))
        stop("'metrics' must be a list of functions")
    return(do.call(cbind, sapply(metrics, function(x) x(data, ...))))
}
