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
#' @param object a \code{\link{comparative.comm}} object
#' @param type character string giving the type of distance matrix on
#' which the mean pairwise distance is based. Either "trait" or "phy"
#' to a phylogenetic or trait-based distance matrix, or an actual
#' matrix to use (e.g., one from \code{\link{funct.phylo.dist}})
#' @param n.sim The number of permutations of the presence vector used to 
#'  make the estimations.
#' @param plot TRUE or FALSE to make the plot of the expected average 
#'  mean pairwise distance, and the 5-95\% confidence interval.
#' @param disp99 Display the 99\% interval?
#' @note No serious checking of user-provided matrices is performed;
#' this is both useful and dangerous!
#' @return TODO
#' @export
ConDivSim<-function(object, type="traits", n.sim=100, plot = TRUE, disp99 = FALSE){
    #Assertions and argument handling
    if(!inherits(object, "comparative.comm"))
        stop("'data' must be a comparative community ecology object")
    Dist <- NULL
    if(type == "traits")
        if(is.null(object$data)) {
            stop("'object' must contain trait data if using a trait matrix")
        } else {
            dist <- as.matrix(dist(object$data))
        }
    if(type == "phy")
        dist <- cophenetic(object$phy)
    if(is.null(dist))
        if(nrow(type)==length(object$phy$tip.label) & nrow(type)==ncol(type))
            dist <- type else stop("'Provided distance matrix is of incorrect dimension'")
    
    
    ## Calculate the observed species richness 
    SpeRich.obs <- rowSums(object$comm)
    Occur.obs <- colSums(object$comm)
	
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
    NbSp <- ncol(object$comm)
    
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
    MPD.obs <- mpd(object$comm, dist, abundance.weighted = FALSE)
    
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
        if(is.null(row.names(object$comm)))
            points(SpeRich.obs, MPD.obs) else
        text(SpeRich.obs, MPD.obs, labels = row.names(object$comm))
        
        
    }
    
    return(Randomiz)
}
