#' Run an eco.phy.regression regression (Cavender-Bares et al. 2004)
#' 
#' \code{eco.phy.regression} runs an eco.phy.regression regression, in the sense of Cavender-Bares et al. (2004), 
#' on a comparative.comm object.
#' 
#' @param data a comparative community ecology object on which to run the regression
#' @param randomisation what kind of null distributions to compare your data with - one of:
#' taxa.labels, richness, frequency, sample.pool, phylogeny.pool, independentswap, trialswap
#' (as implemented in 'picante')
#' @param permute the number of null permutations to perform
#' @details This is extremely unchcked, so beware!
#' @author Will Pearse, Jeannine Cavender-Bares
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample)
#' eco.phy.regression(data, permute=10)
#' }
#' @import picante quantreg vegan
#' @export
eco.phy.regression <- function(data,
  randomisation=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  permute=0, method=c("quantile", "lm", "mantel"), ...){
	#Assertions and argument handling
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
	randomisation <- match.arg(randomisation)
  method <- match.arg(method)
	if(permute < 0) stop("Can't have negative null permutations!")
	
	#Setup matrices
	eco.matrix <- comm.dist(data$comm)
	phy.matrix <- dist(cophenetic(data$phy))
	
	#Observed eco.phy.regression
	observed <- .eco.phy.regression(eco.matrix, phy.matrix, method, ...)
	
	#Randomisations
	randomisations <- vector(mode="list", length=permute)
	#This won't execute if permute is 0...
	for(i in seq(from=1, length.out=permute)){
	  curr.rnd <- .eco.null(data$comm, randomisation)
		rnd.mat <- comm.dist(curr.rnd)
	  if(any(is.na(rnd.mat))){
	    warning("NAs in permuted community matrix; skipping this iteration")
	    next()
	  }
		randomisations[[i]] <- .eco.phy.regression(rnd.mat, phy.matrix, method, ...)
	}
  
  #Prepare output (...and return)
  output <- .prepare.regression.output(observed, randomisations, method, permute, "eco.phy.regression")
  output$data <- data
  return(output)
}

#Perform one set of EcoPhy regressions
.eco.phy.regression <- function(eco.mat, phy.mat, method=c("quantile", "lm", "mantel"), ...){
  method <- match.arg(method)
  if(method == 'lm')
	  model <- lm(as.numeric(eco.mat) ~ as.numeric(phy.mat), ...)
  
  if(method == "quantile")
    model <- rq(as.numeric(eco.mat) ~ as.numeric(phy.mat), ...)
  
  if(method == "mantel")
    model <- mantel(eco.mat, phy.mat, ...)
	
  return(model)
}

#' Prints an eco.phy.regression object (...by summarising it...)
#' @method print eco.phy.regression
#' @S3method print eco.phy.regression
#' @export
print.eco.phy.regression <- function(x, ...){
  summary(x, ...)
}


#' @method summary eco.phy.regression
#' @S3method summary eco.phy.regression
#' @export
summary.eco.phy.regression <- function(x, ...){
  .summary.regression(x, "eco.phy.regression regression")
}

#' @method plot eco.phy.regression
#' @S3method plot eco.phy.regression
#' @export
plot.eco.phy.regression <- function(x, ...){
  eco.matrix <- as.numeric(comm.dist(x$data$comm))
  phy.matrix <- as.numeric(as.dist(cophenetic(x$data$phy)))
  .plot.regression(phy.matrix, eco.matrix, x$observed, x$randomisations, x$method, x$permute,
        xlab="Phylogenetic Distance", ylab="Ecological Co-existnce", ...)
}