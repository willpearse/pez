#' Plot a jcb.plot (sensu Cavender-Bares et al. 2004 Figure 6)
#' 
#' \code{jcb.regression} calculates traits' phylogenetic inertia, and then plots this as a function of trait similarity among communities
#' 
#' @param data a comparative community ecology object on which to run the analyses
#' @param eco.rnd null distribution with which to compare your community data - one of:
#' taxa.labels, richness, frequency, sample.pool, phylogeny.pool, independentswap, trialswap
#' (as implemented in 'picante')
#' @param eco.method regression to perform - one of: lm, quantile, mantel
#' @param eco.permute the number of null permutations to perform
#' @param eco.swap the number of independent swaps to perform (if doing so)
#' @param evo.method how to measure phylogenetic inertia - one of: lambda, delta, kappa, D.c
#' @param evo.permute number of permutations of phylogenetic inertia (D.c only)
#' @details This is extremely unchcked, so beware!
#' @notes Like the eco.trait and eco.env methods, this is a data-hungry method. Warnings will be generated if any of the methods cannot be fitted properly (the examples below give toy examples of this). In such cases the summary and plot methods of these functions may generate errors; use 'traceback()' to examine where these are coming from, and consider whether you want to be working with the data generating these errors. I am loathe to hide these errors or gloss over them, because they represent the reality of your data!
#' @author Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
#' jcb.regression(data, eco.permute=10)
#' plot(jcb.regression(data, permute=10, method="lm))
#' plot(jcb.regression(data, permute=10, method="lm", altogether=FALSE))
#' }
#' @import caper
#' @import picante
#' @import quantreg
#' @import vegan
#' @export
jcb.regression <- function(data, eco.rnd=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  eco.method=c("quantile", "lm", "mantel"), eco.permute=1000, eco.swap=1000, evo.method=c("lambda", "delta", "kappa"), evo.permute=1000, ...){
  #Checks
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  eco.rnd <- match.arg(eco.rnd)
  evo.method <- match.arg(evo.method)
  eco.method <- match.arg(eco.method)
  if(eco.permute < 2) stop("This method relies on random perumtations; you must have at least 2")
  
  #Evolution of traits
  evolution <- phy.signal(data, traits=TRUE, method=evo.method)
  
  #Ecology of traits
  ecology <- eco.trait.regression(data, eco.rnd, eco.permute, eco.method, eco.swap, altogether=FALSE, ...)
  
  #Summarise ecology regressions
  obs.slopes <- numeric(ncol(data$traits))
  rnd.slopes <- matrix(ncol=ncol(data$traits), nrow=ecology$permute)
  for(i in seq(ncol(data$traits))){
    obs.slopes[i] <- ecology[[i]]$obs.slope
    rnd.slopes[,i] <- ecology[[i]]$rnd.slopes
  }
  median.difference <- obs.slopes - apply(rnd.slopes, 2, median)
    
  #Prepare output and return 
  output <- list(evo=evolution, eco=list(raw=ecology, obs.slopes=obs.slopes, rnd.slopes=rnd.slopes, median.diff=median.difference), evo.method=evo.method, eco.method=eco.method)
  class(output) <- "jcb.regression"
  return(output)
}

#' #' Print a jcb.regression (...by summarising it...)
#' @method print jcb.regression
#' @S3method print jcb.regression
#' @export
print.jcb.regression <- function(x, ...){
  summary(x, ...)
}

#' Summarise a jcb.regression
#' @method summary jcb.regression
#' @S3method summary jcb.regression
#' @export
summary.jcb.regression <- function(x, ...){
  cat("Phylogenetic inertia calculated using", x$evo.method, "(examine with model$evo):\n")
  print(summary(x$evo))
  cat("Ecological coexistence calculated calculated using", x$eco.method, "(examine with model$eco):\n")
  print(summary(x$eco$obs.slope))
}

#' Plot a jcb.regression
#' @method plot jcb.regression
#' @S3method plot jcb.regression
#' @export
plot.jcb.regression <- function(x, eco=c("slope", "corrected"), xlab="Ecological Trait Coexistence", ylab="Phylogenetic inertia", ...){
  eco <- match.arg(eco)
  if(eco == "slope")
    eco <- x$eco$obs.slopes else eco <- x$eco$median.diff
  plot(x$evo ~ eco, xlab=xlab, ylab=ylab, ...)
}
