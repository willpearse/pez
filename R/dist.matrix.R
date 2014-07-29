#' Make ecological co-existence matrices
#' @export
comm.dist <- function(x) UseMethod("comm.dist", x)
# - Dispatch on comparative.comm calls the matrix method, which does the work
#' @export
comm.dist.matrix <- function(x){
	output <- matrix(ncol=ncol(x), nrow=ncol(x))
	for(i in seq(ncol(x))){
		for(j in seq(ncol(x))){
			output[i,j] <- 1 - sum(x[,i]>0 & x[,j]>0) / sum(x[,i]>0 | x[,j]>0)
		}
	}
	return(as.dist(output))
}
#' @export
comm.dist.comparative.comm <- function(x) return(comm.dist(x$comm))

#' @export
dist.func.default <- function(x) dist(scale(x, center=TRUE, scale=TRUE))

#' Make trait distance matrices
#' @export
traits.dist <- function(x, dist.func = dist.func.default, ...) UseMethod("traits.dist", x)
#' @export
traits.dist.comparative.comm <- function(x, dist.func = dist.func.default, altogether = TRUE){
    if(is.null(x$traits)) stop("No trait data for which to compute a trait distance matrix")
    if(altogether){
        return(traits.dist(x$traits))
    } else {
                                        # FIXME: is this change ok?  i
                                        # radically changed this
                                        # because i think we should
                                        # pass around dist objects (or
                                        # in this case a list of dist
                                        # objects) --
                                        # https://github.com/willpearse/pez/issues/2#issuecomment-50240752
        return(sapply(as.data.frame(x$traits), dist.func)) 
    }
}
#' @export
traits.dist.numeric <- function(x, dist.func = dist.func.default, ...) dist.func(x)
#' @export
traits.dist.matrix <- function(x, dist.func = dist.func.default, ...) dist.func(x)
#' @export
traits.dist.data.frame <- function(x, dist.func = dist.func.default, ...) dist.func(as.matrix(x))

#' Make phylogenetic distance matrices
#' TODO: not sure if this is necessary/useful yet
#' @export
phylo.dist <- function(x, ...) UseMethod("phylo.dist", x)
#' @export
phylo.dist.phylo <- function(x, ...) as.dist(cophenetic(x))
#' @export
phylo.dist.comparative.comm <- function(x, ...) phylo.dist(x$phy)


#' Make functional phylogenetic distance matrix
#'
#' @param phyloWeight phylogenetic weighting parameter (referred to as
#' \code{a} in Cadotte et al. (2013)
#' @param p exponent giving the norm to use for combining functional
#' and phylogenetic distances (\code{p = 2} gives a Euclidean
#' combination).
#' @export
funct.phylo.dist <- function(x, phyloWeight, p, ...) UseMethod("funct.phylo.dist", x)
#' export
funct.phylo.dist.comparative.comm <- function(x, phyloWeight, p, ...) {
    out <- traits.dist(x)
    out[] <- (phyloWeight * phylo.dist(x)^p +
              (1 - phyloWeight) * out[]^p)^(1/p)
    return(out)
}


#' Make environmental tolerance distance matrices
#' @export
pianka.dist <- function(x, ...) UseMethod("pianka.dist", x)
#' @export
pianka.dist.matrix <- function(comm, env=NULL){
	#Checks and assertions
	if(!is.numeric(comm)) stop("Need numeric community matrix for Pianaka calculations")
	if(!is.factor(env)) stop("Pianaka calculations require a factor as the second argument")
	#Find the proportional occupancy
	propOcc <- matrix(nrow=ncol(comm), ncol=length(levels(env)))
	colnames(propOcc) <- levels(env)
	#Matrices with one row become vectors; this confuses colSums
	for(j in seq(ncol(propOcc)))
		if(sum(env==colnames(propOcc)[j])>1) propOcc[,j] <- colSums(comm[env==colnames(propOcc)[j],]) else propOcc[,j] <- comm[env==colnames(propOcc)[j],]
	propOcc <- apply(propOcc, 2, function(x) x/rowSums(propOcc))
	#Get the Pianka overlap
	pianka <- matrix(ncol=ncol(comm), nrow=ncol(comm))
	#Matrix is symmetrical, so use this to speed things up
	for(j in seq(from=1, to=ncol(comm)-1)){
		for(k in seq(from=j+1, to=ncol(comm))){
			pianka[j,k] <- sum(propOcc[j,] * propOcc[k,]) / sqrt(sum(propOcc[j,]^2) + sum(propOcc[k,]^2))
			pianka[k,j] <- pianka[j,k]
		}
	}
	return(as.dist(pianka))
}
#' @export
pianka.dist.comparative.comm <- function(x, altogether=TRUE){
	#Checks and assertions
	if(is.null(x$env)) stop("Cannot calculate Pianka distances without environmental data")
	if(any(!sapply(x$env, is.factor))) stop("Cannot calculate Pianka distances of environmental data non-factor-level environmental data")
	#Pianaka matrix for each environmental variable
	pianka <- array(dim=c(ncol(x$comm), ncol(x$comm), ncol(x$env)))
	for(i in seq(ncol(x$env)))
		pianka[,,i] <- as.matrix(pianka.dist.matrix(x$comm, x$env[,i]))
  if(altogether)
    return(as.dist(apply(pianka, 1:2, mean))) else return(pianka)
}
