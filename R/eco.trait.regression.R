#' Run an eco.trait.regression regression (Cavender-Bares et al. 2004)
#' 
#' \code{eco.trait.regression} runs a regression of ecological
#' community coexistence versus species trait-distances
#' 
#' @param data a comparative community ecology object on which to run the regression
#' @param randomisation what kind of null distributions to compare your data with - one of:
#' taxa.labels, richness, frequency, sample.pool, phylogeny.pool, independentswap, trialswap
#' (as implemented in 'picante')
#' @param permute the number of null permutations to perform
#' @param method how to compare distance matrices - one of: quantile
#' (quantile regression), lim (linear regression), mantel (Mantel
#' test)
#' @param altogether use distance matrix based on all traits (default
#' TRUE), or perform separate regressions for each trait
#' @param ... additional parameters to pass on to model fitting functions
#' @details This is extremely unchcked, so beware!
#' @author Will Pearse, Jeannine Cavender-Bares
#' @note Like the eco.trait and eco.env methods, this is a data-hungry
#' method. Warnings will be generated if any of the methods cannot be
#' fitted properly (the examples below give toy examples of this). In
#' such cases the summary and plot methods of these functions may
#' generate errors; use 'traceback()' to examine where these are
#' coming from, and consider whether you want to be working with the
#' data generating these errors. I am loathe to hide these errors or
#' gloss over them, because they represent the reality of your data!
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
#' eco.trait.regression(data, permute=10)
#' plot(eco.trait.regression(data, permute=10, method="lm))
#' plot(eco.trait.regression(data, permute=10, method="lm", altogether=FALSE))
#' }
#' @importFrom quantreg rq
#' @importFrom vegan mantel
#' @export
eco.trait.regression <- function(data,
  randomisation=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  permute=0, method=c("quantile", "lm", "mantel"), altogether=TRUE, ...){
	#Assertions and argument handling
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
	randomisation <- match.arg(randomisation)
  method <- match.arg(method)
	if(permute < 0) stop("Can't have negative null permutations!")
	
	#Setup matrices
	eco.matrix <- comm.dist(data$comm)
  if(!is.null(data$data)) traits.matrix <- traits.dist(data, altogether) else stop("'data' must contain trait data for a trait regression!")
	
  #Observed eco.trait.regression
  if(altogether)
      observed <- .eco.trait.regression(eco.matrix, traits.matrix, NULL, method, ...) else {
      #Do separately for all traits
      observed <- vector("list", ncol(data$data))
      for(i in seq(ncol(data$data)))
          observed[[i]] <- .eco.trait.regression(eco.matrix, traits.matrix, i, method, ...)
  }
  
  #Randomisations
  if(altogether){
    #Using mean of traits
  	randomisations <- vector(mode="list", length=permute)
  	#This won't execute if permute is 0...
  	for(i in seq(from=1, length.out=permute)){
            curr.rnd <- .eco.null(data$comm, randomisation)
            rnd.mat <- comm.dist(curr.rnd)
            randomisations[[i]] <- .eco.trait.regression(rnd.mat, traits.matrix, NULL, method, ...)
  	}
  } else {
    #Separately for each trait
    # - preallocate
    randomisations <- vector(mode="list", length=ncol(data$data))
    for(i in seq_along(randomisations)) randomisations[[i]] <- vector("list", permute)
    for(i in seq(from=1, length.out=permute)){
        curr.rnd <- .eco.null(data$comm, randomisation)
        rnd.mat <- comm.dist(curr.rnd)
        for(j in seq(ncol(data$data)))
            randomisations[[j]][[i]] <- .eco.trait.regression(rnd.mat, traits.matrix, j, method, ...)
    }
  }
  
  #Prepare output (...and return)
  if(altogether)
    output <- .prepare.regression.output(observed, randomisations, method, permute, "eco.trait.regression") else {
      output <- vector("list", ncol(data$data))
      for(i in seq_along(output)){
        output[[i]] <- .prepare.regression.output(observed[[i]], randomisations[[i]], method, permute, "eco.trait.regression")
        output[[i]]$altogether <- altogether
      }
      output$type <- "eco.trait.regression"
      class(output) <- "ecophyl.regression.list"
    }
  output$data <- data
  output$altogether <- altogether
  output$permute <- permute;output$method<-method
  return(output)
}


#Perform one set of EcoPhy regressions
.eco.trait.regression <- function(eco.mat, trait.mat, which.trait=NULL, method=c("quantile", "lm", "mantel"), ...){
  method <- match.arg(method)
  #Check to see if we're doing this across all traits
  # - now we're passing around distance matrices by default you have to be extra careful...
  if(!is.null(which.trait) & !inherits(trait.mat, "dist"))
      trait.mat <- as.dist(trait.mat[,,which.trait])
  
  if(method == 'lm')
    model <- lm(as.numeric(eco.mat) ~ as.numeric(trait.mat), ...)
  
  if(method == "quantile")
    model <- rq(as.numeric(eco.mat) ~ as.numeric(trait.mat), ...)
  
  if(method == "mantel")
    model <- mantel(eco.mat, trait.mat, ...)
  
  return(model)
}
