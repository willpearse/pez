#' Functional-phylogenetic community simulations
#'
#' Simulations for assessing the relative importance of phylogenetic 
#' versus functional information for understanding variation in
#' community composition.
#'
#' Simulations are based on the following procedure:
#' \describe{
#'
#'	\item{Phylogeny}{A tree is simulated using \code{rcoal} with
#'	default settings.}
#'
#'	\item{Trait evolution}{Two traits are simulated under Brownian motion. 
#'	One trait is called 'observed' and the other 'unknown'. Both traits
#'  influence species' probabilies of occurrence, but only the 'observed'
#'	trait is included in the output. The idea behind this distinction is that
#'	information about the 'unknown' trait will possibly be present in
#'	the phylogeny, which is assumed known.}
#'
#'	\item{Probabilities of occurrence}{For each species at each site, these
#'	probabilities depend logistically on two gradients; one is observed
#'	and the other is unknown, each corresponding to the observed and unknown
#'	traits (see the \code{sim.envObserved} and \code{sim.envUnknown} arguments).
#'	The corresponding traits are the logit-scale slopes of these logistic
#'	curves -- in other words each species depends on the gradient differently
#'	depending on their traits. When \code{p} equals one (zero), only the observed
#'	(unknown) gradient and trait determine probability of occurrence. When 
#'	\code{p} is between zero and one, both the observed and unknown
#'	trait-gradient combinations have some influence.}
#'
#'	\item{Occurrence}{Each species is present at each site with probability
#'	given by the corresponding probability of occurrence.}
#'
#' }
#'
#' @param n Number of sites to simulate.
#' @param m Number of species to simulate.
#' @param p Number between 0 and 1 giving the relative importance of
#'	observed versus unknown traits in determining community structure.
#' @param diverg.obs A vector with m elements giving the divergence
#'	effects on the observed trait for each species.
#' @param diverg.unk A vector with m elements giving the divergence
#'	effects on the unknown trait for each species.
#' @param sim.envObserved A numeric vector of simulated values for the 
#'	observed environmental variables.
#' @param sim.envUnknown A numeric vector of simulated values for the 
#'	unknown environmental variables.
#' @param site.names Character vector of site names.
#' @param spp.names Character vector of species names.
#' @return An object of class \code{fpcomSims} with components:
#'	\item{comm}{An n-by-m matrix of presences and absences}
#'	\item{probs}{An n-by-m matrix of probabilities of occurrence}
#'	\item{traits}{A length-m vector of values for the observed trait}
#'	\item{env}{A length-n vector of values for the observed gradient}
#'	\item{tree}{A phylo object with the phylogenetic tree relating
#'	the m species}
#' @export
fpcomSims <- function(n, m, p = 0.5,
	diverg.obs = rep(0, m),
	diverg.unk = rep(0, m),
	sim.envObserved = as.vector(scale(rnorm(n))),
	sim.envUnknown = as.vector(scale(rnorm(n))),
	site.names = numnames(n,"site"), 
	spp.names = numnames(m,"sp")
){
	
	names(sim.envObserved) <- names(sim.envUnknown) <- site.names
	
	# simulate evolutionary history
	sim.tree <- rcoal(m, tip.label = spp.names)
	sim.tree$tip.label <- spp.names
	
	# store the resulting traits for the species
	sim.traits.obs <- as.vector(t(chol(vcv(sim.tree))) %*% rnorm(m, sd = 1))
	sim.traits.unk <- as.vector(t(chol(vcv(sim.tree))) %*% rnorm(m, sd = 1))
	sim.traits.obs <- sim.traits.obs + diverg.obs
	sim.traits.unk <- sim.traits.unk + diverg.unk
	names(sim.traits.obs) <- spp.names

	# simulate the contemporary communities
	sim.eta <- (p*outer(sim.envObserved, sim.traits.obs)) + 
			   ((1-p)*outer(sim.envUnknown, sim.traits.unk))
	sim.eta <- scale(as.vector(sim.eta))
	dim(sim.eta) <- c(n, m)
	sim.p <- exp(sim.eta)/(1+exp(sim.eta))
	sim.comm <- matrix(rbinom(n*m, 1, sim.p), n, m)
	dimnames(sim.comm) <- dimnames(sim.p) <- list(site.names, spp.names)

	ecosysfunc <- (p*sim.envObserved) + ((1-p)*sim.envUnknown)

	out <- list(comm = sim.comm, 
		probs = sim.p,
		traits = structure(sim.traits.obs, unknown = sim.traits.unk),
		env = sim.envObserved,
		ecosysfunc = ecosysfunc,
		tree = sim.tree,
		tuning = p)
	class(out) <- "fpcomSims"
	return(out)
}

update.fpcomSims <- function(object, nnew = 1, 
	sim.envObserved = rnorm(nnew),
	sim.envUnknown = rnorm(nnew),
	site.names = numnames(length(object$env) + nnew, "site"),
	spp.names = numnames(length(object$traits),"sp"), ...){

	m <- ncol(object$probs)
	
	sim.traits.obs <- structure(object$traits, unknown = NULL)
	sim.traits.unk <- attr(object$traits, 'unknown')
	p <- object$tuning

	# simulate the new contemporary communities
	sim.eta <- (p*outer(sim.envObserved, sim.traits.obs)) + 
			   ((1-p)*outer(sim.envUnknown, sim.traits.unk))
	sim.eta <- scale(as.vector(sim.eta))
	dim(sim.eta) <- c(nnew, m)
	sim.p <- exp(sim.eta)/(1+exp(sim.eta))
	sim.comm <- matrix(rbinom(nnew*m, 1, sim.p), nnew, m)
	
	sim.p <- rbind(object$probs, sim.p)
	sim.comm <- rbind(object$comm, sim.comm)
	sim.envObserved <- c(object$env, sim.envObserved)
		
	dimnames(sim.comm) <- dimnames(sim.p) <- list(site.names, spp.names)
	names(sim.envObserved) <- site.names
	
	ecosysfunc <- sim.envObserved + sim.envUnknown
	
	out <- list(comm = sim.comm, 
		probs = sim.p,
		traits = object$traits,
		env = sim.envObserved,
		ecosysfunc = ecosysfunc,
		tree = object$tree,
		tuning = p)
	class(out) <- "fpcomSims"
	return(out)
}

#' Plot fpcomSims objects
#'
#' This plot method can produce four different graphical summaries of the
#' output of the \code{\link{fpcomSims}} function.
#'
#' \describe{
#'
#'	\item{\code{plottype} equals \code{"distance"}}{A scatterplot
#'	of the phylogenetic versus functional distances between the species 
#'	pairs are produced. Species numbers are used to label the points. The
#'	\code{\link{cophenetic}} function is used to compute the phylogenetic
#'	distances and the \code{\link{dist}} function is used to compute the
#'	functional Euclidean distances. Both distances are standardised such
#'	that the maximum distance is one.}
#'	
#'	\item{\code{plottype} equals \code{"traitgram"}}{A \code{\link{traitgram}}
#'	is produced, which combines both the phylogeny and the observed trait.}
#'
#'	\item{\code{plottype} equals \code{"gradient"}}{A scatterplot of the 
#'	probabilities of occurrence versus the observed gradient is produced.
#'	The species are identified by their numbers. Note that these probabilities
#'	of occurrence depend only on the two gradients and the two traits, and
#'	therefore are not subject to variation among sites with identical gradient
#'	values. Such variation is expressed in the \code{comm} element of
#'	\code{fpcomSims} objects, which gives not probability of occurrence but
#'	occurrence itself. Note that if the \code{p} argument to \code{fpcomSims}
#'	is one, then only the observed trait and gradient determine probability of
#'	occurrence and the plot consists of perfect sigmoid curves.  But if \code{p}
#'	is zero then only the unknown gradient and trait determine probability of
#'	occurrence and the plot is just noise. In this latter case we would expect
#'	phylogenetic distance to provide more important information about communities.}
#'
#'	\item{\code{plottype} equals \code{"ordination"}}{A special type of 
#'	ordination of the species is produced. This ordination is based on the
#'	probabilities of occurrence and not on occurrence itself, and so it is able
#'	to fully explain all variation on two axes (if \code{p} is not either zero
#'	of one). Therefore, such an ordination is not possible in practice with real
#'	data but it is useful in this case for fully representing distances between
#'	species in species distribution space. Technically it is a singular value
#'	decomposition of the logit transformed probabilities of occurrence. Note
#'	well the percentages of variation explained by the two axes, and that this
#'	ordination is gauranteed to have 100 percent of the variation explained by
#'	the first two axes.}
#' }
#'
#' @param x An \code{fpcomSims} object.
#' @param y Not used at the moment.
#' @param cex.add The cex graphics parameter for additional graphical. 
#'	elements not produced by the \code{\link{plot}} command.
#' @param plottype The type of plot to produce -- see details.
#' @param ... Additional parameters to be passed to \code{\link{plot}}.
#' @method plot fpcomSims
#' @return No return value, called for its side-effect of producing a plot.
#' @export
plot.fpcomSims <- function(x,y,cex.add=0.8,
	plottype=c("distance","traitgram","gradient","ordination"),...){
	
	if(plottype[1]=="distance"){
		PDist <- cophenetic(x$tree)/max(cophenetic(x$tree))
		PDist <- PDist[order(row.names(PDist)),order(row.names(PDist))]
		FDist <- as.matrix(dist(x$traits))/max(as.matrix(dist(x$traits)))
		plot(PDist,FDist,type="n",
			xlab="phylogenetic distance",
			ylab="functional distance",
			...)
		text(as.dist(PDist),as.dist(FDist),
			paste(as.dist(col(PDist)),as.dist(row(FDist)),sep=","),
			cex=cex.add)
	}
	
	if(plottype[1]=="traitgram"){
		traitgram(x$traits,x$tree,...)
	}
	
	if(plottype[1]=="gradient"){
		plot(x$env,x$probs[,1],type="n",ylim=c(0,1),
			xlab="gradient",
			ylab="probability of occurrence",
			...)
		for(i in 1:ncol(x$comm)){
			text(x$env,x$probs[,i],i,cex=cex.add)
		}
	}
	
	if(plottype[1]=="ordination"){
		x.ord <- svd(log(x$probs/(1-x$probs)))
		plot(x.ord$v[,1:2]%*%diag(x.ord$d[1:2]),type="n",asp=1,
			xlab=paste("Ordination axis I (",
				100*round(x.ord$d[1]/sum(x.ord$d),2),"%)",sep=""),
			ylab=paste("Ordination axis II (",
				100*round(x.ord$d[2]/sum(x.ord$d),2),"%)",sep=""),
			...)
		text(x.ord$v[,1:2]%*%diag(x.ord$d[1:2]),labels=names(x$traits),cex=cex.add)
	}
}



#' Generates numbered names
#'
#' Generates numbered names based on a character prefix that will
#' be alphanumerically sorted as expected.
#'
#' @param n Number of names.
#' @param prefix Character prefix to go in front of the numbers.
#' @return A character vector of numbered names.
numnames <- function(n, prefix = "name"){
	n <- as.integer(n)
	if(n < 1) stop("number of names must be one or more")
	zeropad.code <- paste("%0", nchar(n), ".0f", sep = "")
	numpart <- sprintf(zeropad.code, 1:n)
	paste(prefix, numpart, sep = "")
}

#' Functional phylogenetic distance
#'
#' Computates a functional phylogenetic distance matrix from 
#' phylogenetic and functional distance matrices and two
#' weighting parameters.
#'
#' @param data comparative.comm object with traits and a phylogeny
#' @param a A number between 0 and 1 giving the amount of weight.
#'	to put on \code{PDist} relative to \code{FDist}.
#' @param p A number giving the \code{p}-norm.
#' @return A distance matrix.
#' @note This function is not very user friendly yet.  There are
#'	no doubt many use cases that I've ignored.
#' @export
FPDist <- function(data, a=0.5, p=1){
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(is.null(data$traits)) stop("'data' must contain trait data")
    if(a < 0 | a > 1) stop("'a' must be between 0 and 1")
    if(!is.numeric(p)) stop("'p' must be a numeric")
    
    PDist <- cophenetic(data$phy)
    FDist <- as.matrix(dist(data$traits), rownames.force = TRUE)
	
    FDist <- FDist/max(FDist)
    PDist <- PDist/max(PDist)
    
    return(((a*(PDist^p)) + ((1-a)*(FDist^p)))^(1/p))
}


#' Null models for functional-phylogenetic diversity
#'
#' Simulate expectations (under a null model) of mean pairwise distance for 
#' a set of communities with different species richness.
#'
#' If \code{plot == TRUE}, then a surface is drawn giving the
#' null distribution.  Lighter shades of gray give larger intervals 
#' with categories: 0.005-0.995 = 99\%, 0.025-0.975 = 95\%, 0.05-0.95 = 90\%,
#' 0.25-0.75 = 50\%.
#'
#' @param tab The community by species presence table.
#' @param type The distance matrix on which the mean pairwise distance
#' is based. Either "trait" or "phy" to a phylogenetic or trait-based
#' distance matrix, or an actual matrix to use (e.g., one from
#' \code{\link{FPDist}})
#' @param n.sim The number of permutations of the presence vector used to 
#'  make the estimations.
#' @param plot TRUE or FALSE to make the plot of the expected average 
#'  mean pairwise distance, and the 5-95\% confidence interval.
#' @param disp99 Display the 99\% interval?
#' @note No serious checking of user-provided matrices is performed;
#' this is both useful and dangerous!
#' @return TODO
#' @export
ConDivSim<-function(data, type="traits", n.sim=100, plot = TRUE, disp99 = FALSE){
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    Dist <- NULL
    if(type == "traits")
        if(is.null(data$traits)) stop("'data' must contain trait data if using a trait matrix") else dist <- as.matrix(dist(data$traits))
    if(type == "phy")
        dist <- copheneitc(data$phy)
    if(is.null(dist))
        if(nrow(type)==length(data$phy$tip.label) & nrow(type)==ncol(type))
            dist <- type else stop("'Provided distance matrix is of incorrect dimension'")
    
    
    ## Calculate the observed species richness 
    SpeRich.obs <- rowSums(data$comm)
    Occur.obs <- colSums(data$comm)
	
    ## Check that numbers are high enough for randomisation
    ## (e.g. no communities with one species only)
    if(any(SpeRich.obs < 2)) 
        stop('all communities must have more than one species present')
    if(any(Occur.obs < 1)) 
        stop('all species must be present at more than zero communities')
    
    ## Calculate the range of species richness in the communities
    if(min(SpeRich.obs)>2){
        SpRich <- (min(SpeRich.obs)-1):(max(SpeRich.obs)+1)
    }
    else {
        SpRich <- min(SpeRich.obs):(max(SpeRich.obs)+1)
    }

    ## Calculate the regional species richness(gamma)
    NbSp <- ncol(data$comm)
    
    ## Create the matrix in which to store the results
    Randomiz <- matrix(0, length(SpRich), 11)
    colnames(Randomiz)<-c("SpRich","ExpMeanMPD",
                          "ExpQ0.005MPD","ExpQ0.025MPD",
                          "ExpQ0.05MPD","ExpQ0.25MPD",
                          "ExpQ0.75MPD","ExpQ0.95MPD",
                          "ExpQ0.975MPD","ExpQ0.995MPD",
                          "ExpSdMPD")
    row.names(Randomiz) <- 1:length(SpRich)
    
    ## calculate the observed mean pairwise distance for the community
    ## (with the unweighted method from mpd in picante, it means only 
    ## the lower triangle!)
    MPD.obs <- mpd(data$comm, dist, abundance.weighted = FALSE)
    
    ## Loop on the species richness range to calculate the random 
    ## expectations of mean pairwise distance
    for (k in 1:length(SpRich)){
        
    	# Vector of local richness within the possible range SpRich
    	locrich <- SpRich[k]
        
    	# Create the basic vector of presence
        Vec.Pres <- c(rep(1, locrich), rep(0, NbSp-locrich))
        
        # Estimate the random expectations of mean pairwise distances 
        # for the random community with locrich as the species richness
        Sim.Div <- NULL
        
        for (i in 1:n.sim) {
            Vec.Pres.Sim <- sample(Vec.Pres)
            sample.dis <- dist[as.logical(Vec.Pres.Sim), as.logical(Vec.Pres.Sim)]
            
            # line in mpd from picante with abundance.weighted = FALSE
            Sim.Div[i] <- mean(sample.dis[lower.tri(sample.dis)])
      	}
        
        # Store the results in Randomiz
        Randomiz[k,1]<-locrich
        Randomiz[k,2]<-mean(Sim.Div)
        Randomiz[k, 3:10] <- quantile(Sim.Div, c(0.005, 0.025, 0.05, 0.25, 0.75, 0.95, 0.975, 0.995))
    	Randomiz[k,11]<-sd(Sim.Div)
    }
	
    ## Make the graph if required
    if(plot == TRUE){
        
        plot(Randomiz[,1:2],
             ylim = c(min(Randomiz[, 3], MPD.obs), max(Randomiz[, 10], MPD.obs)),
             type="n", xlab="Species richness",ylab="Mean pairwise distance")
        
        ## Surface for each of the quantile pair:
        if(disp99){
            polygon( ## 0.005-0.995 = 99%
                    c(Randomiz[,1], Randomiz[length(SpRich):1, 1]), 
                    c(Randomiz[,3], Randomiz[length(SpRich):1, 10]),
                    col=grey(0.9), border=NA)
    	}
    	
        polygon( ## 0.025-0.975 = 95%
                c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
                c(Randomiz[,4], Randomiz[length(SpRich):1,9]),
                col=grey(0.8),border=NA)
    	
        polygon( ## 0.05-0.95 = 90%
                c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
                c(Randomiz[,5], Randomiz[length(SpRich):1,8]),
                col=grey(0.7),border=NA)
        
        polygon( ## 0.25-0.75 = 50%
                c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
                c(Randomiz[,6], Randomiz[length(SpRich):1,7]),
                col=grey(0.6),border=NA)
    	
        ## Average of the expected mean pairwise distance
        lines(Randomiz[,1:2], col = "black", lwd = 2)
        
        ## Labelling the communities
        if(is.null(row.names(data$comm)))
            points(SpeRich.obs, MPD.obs) else
        text(SpeRich.obs, MPD.obs, labels = row.names(data$comm))
        
        
    }
    
    return(Randomiz)
}
