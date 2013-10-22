#' Make ecological co-existence matrices
#' @export
comm.dist <- function(x) UseMethod("comm.dist", x)
# - Dispatch on comparative.comm calls the matrix method, which does the work
comm.dist.matrix <- function(x){
	output <- matrix(ncol=ncol(x), nrow=ncol(x))
	for(i in seq(ncol(x))){
		for(j in seq(ncol(x))){
			output[i,j] <- 1 - sum(x[,i]>0 & x[,j]>0) / sum(x[,i]>0 | x[,j]>0)
		}
	}
	return(as.dist(output))
}
comm.dist.comparative.comm <- function(x) return(comm.dist(x$comm))


#' Make trait distance matrices
#' @export
traits.dist <- function(x, ...) UseMethod("traits.dist", x)
traits.dist.comparative.comm <- function(x, altogether=TRUE){
	if(is.null(x$traits)) stop("No trait data for which to compute a trait distance matrix")
  if(altogether){
    mat <- as.matrix(x$traits)
    return(traits.dist(mat))
  } else {
    traits <- array(dim=c(nrow(x$traits), nrow(x$traits), ncol(x$traits)))
    for(i in seq(ncol(x$traits)))
      traits[,,i] <- as.matrix(traits.dist(x$traits[,i]))
	  return(traits)
  }
}
traits.dist.numeric <- function(x){
  trait <- scale(x, center=TRUE, scale=TRUE)
  return(dist(trait))
}
traits.dist.matrix <- function(x){
  if(!is.numeric(x)) stop("Can only compute trait distance matrix of continuous traits")
  mat <- scale(x, center=TRUE, scale=TRUE)
  return(dist(mat))
}

#' Make environmental tolerance distance matrices
#' @export
pianka.dist <- function(x, ...) UseMethod("pianka.dist", x)
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