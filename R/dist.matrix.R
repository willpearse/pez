#' Make ecological co-existence matrices
#' @param x an object
#' @export
#' @details Returns the 1 - co-existence of species. Look at how this
#' is calcualted; it incorporates abundances, and if you don't want it
#' to do so simply call it on a presence/absensence (1/0) matrix.
#' @rdname comm.dist
#' @family distances
comm.dist <- function(x) UseMethod("comm.dist", x)
# - Dispatch on comparative.comm calls the matrix method, which does the work
#' @export
#' @rdname comm.dist
comm.dist.matrix <- function(x){
	output <- matrix(ncol=ncol(x), nrow=ncol(x))
	for(i in seq(ncol(x))){
		for(j in seq(ncol(x))){
			output[i,j] <- 1 - sum(x[,i]>0 & x[,j]>0) / sum(x[,i]>0 | x[,j]>0)
		}
	}
        if(any(!is.finite(output))){
            output[!is.finite(output)] <- 1
            warning("Co-existence matric contains non-overlapping species; treating as maximally dissimilar")
        }
        return(as.dist(output))
}
#' @export
#' @rdname comm.dist
comm.dist.comparative.comm <- function(x) return(comm.dist(x$comm))

#' Make trait distance matrices
#'
#' @param x an object
#' @param dist.func a function for computing distances.  The default,
#' \code{dist.func.default}, returns a Euclidean distance of the
#' scaled and centred data.
#' @param alltogether should one multivariate distance matrix be
#' computed for all traits at once (\code{alltogether = TRUE}) or for
#' each trait at a time (\code{alltogether = FALSE})?
#' @param ... not used
#' @rdname traits.dist
#' @family distances
#' @export
traits.dist <- function(x, dist.func = dist.func.default, ...) UseMethod("traits.dist", x)
#' @rdname traits.dist
#' @export
traits.dist.comparative.comm <- function(x, dist.func = dist.func.default, alltogether = TRUE, ...){
    if(is.null(x$data)) stop("No trait data for which to compute a trait distance matrix")
    if(alltogether){
        return(traits.dist(x$data))
    } else {
        return(sapply(as.data.frame(x$data), dist.func)) 
    }
}
#' @export
#' @rdname traits.dist
traits.dist.default <- function(x, dist.func = dist.func.default, ...) dist.func(x)
#' @export
#' @rdname traits.dist
traits.dist.data.frame <- function(x, dist.func = dist.func.default, ...) dist.func(as.matrix(x))
#' @export
#' @rdname traits.dist
dist.func.default <- function(x) dist(scale(x, center=TRUE, scale=TRUE))

#' Make phylogenetic distance matrices
#'
#' @param x an object
#' @param ... not used
#' @export
#' @rdname phylo.dist
#' @family distances
phylo.dist <- function(x, ...) UseMethod("phylo.dist", x)
#' @export
#' @rdname phylo.dist
phylo.dist.phylo <- function(x, ...) as.dist(cophenetic(x))
#' @export
#' @rdname phylo.dist
phylo.dist.comparative.comm <- function(x, ...) phylo.dist(x$phy)


#' Make functional phylogenetic distance matrix
#'
#' @param x a \code{\link{comparative.comm}} object
#' @param phyloWeight phylogenetic weighting parameter (referred to as
#' \code{a} in Cadotte et al. (2013)
#' @param p exponent giving the exponent to use for combining
#' functional and phylogenetic distances (the default, \code{p = 2},
#' gives a Euclidean combination).
#' @param ... not currently used
#' @rdname funct.phylo.dist
#' @family distances
#' @export
funct.phylo.dist <- function(x, phyloWeight, p = 2, ...) UseMethod("funct.phylo.dist", x)
#' @export
funct.phylo.dist.comparative.comm <- function(x, phyloWeight, p, ...) {
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(is.null(data$data)) stop("'data' must contain trait data")
    if(phyloWeight < 0 | phyloWeight > 1) stop("'a' must be between 0 and 1")
    if(!is.numeric(p)) stop("'p' must be a numeric")
    
    FDist <- traits.dist(x)
    PDist <- phylo.dist(x)
    FDist <- FDist/max(FDist)
    PDist <- PDist/max(PDist)
    (phyloWeight * PDist^p + (1 - phyloWeight) * FDist^p)^(1/p)
}


#' Make environmental tolerance distance matrices
#'
#' @param x an object
#' @param ... (ignored)
#' @export
#' @rdname pianka.dist
#' @details Return's the Pianka's distance; i.e., niche overlap based
#' on environemntal tolerances contained within your dataset.
#' @family distances
pianka.dist <- function(x, ...) UseMethod("pianka.dist", x)

#' @rdname pianka.dist
#' @param env environmental variable to be used to calculate the
#' distance matrix
#' @export
pianka.dist.matrix <- function(x, env = NULL, ...){
    #Checks and assertions
    if(!is.numeric(x)) stop("Need numeric community matrix for Pianaka calculations")
    if(!is.factor(env)) stop("Pianaka calculations require a factor as the second argument")
    #Find the proportional occupancy
    propOcc <- matrix(nrow=ncol(x), ncol=length(levels(env)))
    colnames(propOcc) <- levels(env)
    #Matrices with one row become vectors; this confuses colSums
    for(j in seq(ncol(propOcc)))
        if(sum(env==colnames(propOcc)[j])>1) propOcc[,j] <- colSums(x[env==colnames(propOcc)[j],]) else propOcc[,j] <- x[env==colnames(propOcc)[j],]
    propOcc <- apply(propOcc, 2, function(x) x/rowSums(propOcc))
    #Get the Pianka overlap
    pianka <- matrix(ncol=ncol(x), nrow=ncol(x))
    #Matrix is symmetrical, so use this to speed things up
    for(j in seq(from=1, to=ncol(x)-1)){
        for(k in seq(from=j+1, to=ncol(x))){
            pianka[j,k] <- sum(propOcc[j,] * propOcc[k,]) / sqrt(sum(propOcc[j,]^2) + sum(propOcc[k,]^2))
            pianka[k,j] <- pianka[j,k]
        }
    }
    return(as.dist(pianka))
}

#' @export
#' @param alltogether whether to calculate distance matrices for all
#' the environmental variables together (default: TRUE) or calculate
#' and return each separately
#' @rdname pianka.dist
pianka.dist.comparative.comm <- function(x, alltogether = TRUE, ...){
	#Checks and assertions
	if(is.null(x$env)) stop("Cannot calculate Pianka distances without environmental data")
	if(any(!sapply(x$env, is.factor))) stop("Cannot calculate Pianka distances of environmental data non-factor-level environmental data")
	#Pianaka matrix for each environmental variable
	pianka <- array(dim=c(ncol(x$comm), ncol(x$comm), ncol(x$env)))
	for(i in seq(ncol(x$env)))
		pianka[,,i] <- as.matrix(pianka.dist.matrix(x$comm, x$env[,i]))
  if(alltogether)
    return(as.dist(apply(pianka, 1:2, mean))) else return(pianka)
}
