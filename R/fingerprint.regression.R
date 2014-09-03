#' Plot a fingerprint regression (sensu Cavender-Bares et al. 2004 Figure 6)
#' 
#' \code{fingerprint.regression} calculates traits' phylogenetic
#' inertia, and then plots this as a function of trait similarity
#' among communities
#' 
#' @param data a comparative community ecology object on which to run the analyses
#' @param eco.rnd null distribution with which to compare your community data - one of:
#' taxa.labels, richness, frequency, sample.pool, phylogeny.pool, independentswap, trialswap
#' (as implemented in 'picante')
#' @param eco.method regression to perform - one of: lm, quantile, mantel
#' @param eco.permute the number of null permutations to perform
#' @param evo.method how to measure phylogenetic inertia - one of:
#' lambda, delta, kappa.
#' @param eco.swap number of independent swap iterations to perform
#' (if using that randomisation); default is 1000
#' @param ... additional parameters to pass on to model fitting functions
#' @details This is extremely unchcked, so beware!
#' @note Like the eco.trait and eco.env methods, this is a data-hungry
#' method. Warnings will be generated if any of the methods cannot be
#' fitted properly (the examples below give toy examples of this). In
#' such cases the summary and plot methods of these functions may
#' generate errors; use 'traceback()' to examine where these are
#' coming from, and consider whether you want to be working with the
#' data generating these errors. I am loathe to hide these errors or
#' gloss over them, because they represent the reality of your data!
#' @author Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
#' fingerprint.regression(data, eco.permute=10)
#' plot(fingerprint.regression(data, permute=10, method="lm))
#' plot(fingerprint.regression(data, permute=10, method="lm", altogether=FALSE))
#' }
#' @export
fingerprint.regression <- function(data, eco.rnd=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  eco.method=c("quantile", "lm", "mantel"), eco.permute=1000, evo.method=c("lambda", "delta", "kappa"), eco.swap=1000, ...){
  #Checks
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  eco.rnd <- match.arg(eco.rnd)
  evo.method <- match.arg(evo.method)
  eco.method <- match.arg(eco.method)
  if(eco.permute < 2) stop("This method relies on random perumtations; you must have at least 2")
  
  #Evolution of traits
  evolution <- phy.signal(data, traits=TRUE, method=evo.method)
  
  #Ecology of traits
  ecology <- eco.trait.regression(data, eco.rnd, eco.permute, eco.method, altogether=FALSE, ...)
  
  #Summarise ecology regressions
  obs.slopes <- numeric(ncol(data$data))
  rnd.slopes <- matrix(ncol=ncol(data$data), nrow=ecology$permute)
  for(i in seq(ncol(data$data))){
    obs.slopes[i] <- ecology[[i]]$obs.slope
    rnd.slopes[,i] <- ecology[[i]]$rnd.slopes
  }
  median.difference <- obs.slopes - apply(rnd.slopes, 2, median)
    
  #Prepare output and return 
  output <- list(evo=evolution, eco=list(raw=ecology, obs.slopes=obs.slopes, rnd.slopes=rnd.slopes, median.diff=median.difference), evo.method=evo.method, eco.method=eco.method)
  class(output) <- "fingerprint.regression"
  return(output)
}

#' #' Print a fingerprint.regression (...by summarising it...)
#' @method print fingerprint.regression
#' @param x \code{fingerprint.regression} object
#' @param ... arguments passed to
#' \code{summary.fingerprint.regression} (currently ignored)
#' @export
print.fingerprint.regression <- function(x, ...){
  summary(x, ...)
}

#' Summarise a fingerprint.regression
#' @method summary fingerprint.regression
#' @param object \code{fingerprint.regression} object
#' @param ... ignored
#' @export
summary.fingerprint.regression <- function(object, ...){
  cat("Phylogenetic inertia calculated using", object$evo.method, "(examine with model$evo):\n")
  print(summary(object$evo))
  cat("Ecological coexistence calculated calculated using", object$eco.method, "(examine with model$eco):\n")
  print(summary(object$eco$obs.slope))
}

#' Plot a fingerprint.regression
#' @method plot fingerprint.regression
#' @param x \code{fingerprint.regression} object
#' @param eco plot the observed slopes ("slope", the default), or the
#' median difference between the simulations and the observed values
#' ('corrected')
#' @param xlab label for x-axis (default "Ecological Trait Coexistence")
#' @param ylab label for y-axis (default "Phylogenetic inertia")
#' @param ... additional plotting arguments
#' @export
plot.fingerprint.regression <- function(x, eco=c("slope", "corrected"), xlab="Ecological Trait Coexistence", ylab="Phylogenetic inertia", ...){
  eco <- match.arg(eco)
  if(eco == "slope")
    eco <- x$eco$obs.slopes else eco <- x$eco$median.diff
  plot(x$evo ~ eco, xlab=xlab, ylab=ylab, ...)
}
