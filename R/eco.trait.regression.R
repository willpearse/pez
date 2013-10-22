#' Run an eco.trait.regression regression (Cavender-Bares et al. 2004)
#' 
#' \code{eco.trait.regression} runs a regression of ecological community coexistence versus species trait-distances
#' 
#' @param data a comparative community ecology object on which to run the regression
#' @param randomisation what kind of null distributions to compare your data with - one of:
#' taxa.labels, richness, frequency, sample.pool, phylogeny.pool, independentswap, trialswap
#' (as implemented in 'picante')
#' @param permute the number of null permutations to perform
#' @param altogether use distance matrix based on all traits (default TRUE), or perform separate regressions for each trait
#' @details This is extremely unchcked, so beware!
#' @author Will Pearse, Jeannine Cavender-Bares
#' @note Like the eco.trait and eco.env methods, this is a data-hungry method. Warnings will be generated if any of the methods cannot be fitted properly (the examples below give toy examples of this). In such cases the summary and plot methods of these functions may generate errors; use 'traceback()' to examine where these are coming from, and consider whether you want to be working with the data generating these errors. I am loathe to hide these errors or gloss over them, because they represent the reality of your data!
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
#' eco.trait.regression(data, permute=10)
#' plot(eco.trait.regression(data, permute=10, method="lm))
#' plot(eco.trait.regression(data, permute=10, method="lm", altogether=FALSE))
#' }
#' @import picante quantreg vegan
#' @export
eco.trait.regression <- function(data,
  randomisation=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  permute=0, method=c("quantile", "lm", "mantel"), swap.iter=1000, altogether=TRUE, ...){
	#Assertions and argument handling
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
	randomisation <- match.arg(randomisation)
  method <- match.arg(method)
	if(permute < 0) stop("Can't have negative null permutations!")
	
	#Setup matrices
	eco.matrix <- comm.dist(data$comm)
  if(!is.null(data$traits)) traits.matrix <- traits.dist(data, altogether) else stop("'data' must contain trait data for a trait regression!")
	
	#Observed eco.trait.regression
  if(altogether)
	  observed <- .eco.trait.regression(eco.matrix, traits.matrix, NULL, method, ...) else {
      #Do separately for all traits
      observed <- vector("list", ncol(data$traits))
      for(i in seq(ncol(data$traits)))
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
  		randomisations[[i]] <- .eco.trait.regression(eco.matrix, traits.matrix, NULL, method, ...)
  	}
  } else {
    #Separately for each trait
    # - preallocate
    randomisations <- vector(mode="list", length=ncol(data$traits))
    for(i in seq_along(randomisations)) randomisations[[i]] <- vector("list", permute)
    for(i in seq(from=1, length.out=permute)){
      curr.rnd <- .eco.null(data$comm, randomisation)
      rnd.mat <- comm.dist(curr.rnd)
      for(j in seq(ncol(data$traits)))
        randomisations[[j]][[i]] <- .eco.trait.regression(eco.matrix, traits.matrix, j, method, ...)
    }
  }
  
  #Prepare output (...and return)
  if(altogether)
    output <- .prepare.regression.output(observed, randomisations, method, permute, "eco.trait.regression") else {
      output <- vector("list", ncol(data$traits))
      for(i in seq_along(output)){
        output[[i]] <- .prepare.regression.output(observed[[i]], randomisations[[i]], method, permute, "eco.trait.regression")
        output[[i]]$altogether <- altogether
      }
      class(output) <- "eco.trait.regression.list"
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
  if(!is.null(which.trait)) trait.mat <- as.dist(trait.mat[,,which.trait])
  
  if(method == 'lm')
    model <- lm(as.numeric(eco.mat) ~ as.numeric(trait.mat), ...)
  
  if(method == "quantile")
    model <- rq(as.numeric(eco.mat) ~ as.numeric(trait.mat), ...)
  
  if(method == "mantel")
    model <- mantel(eco.mat, trait.mat, ...)
  
  return(model)
}

#\code{print.eco.trait.regression} prints an eco.trait.regression regression (...by summarising it...))
#' @method print eco.trait.regression
#' @S3method print eco.trait.regression
#' @export
print.eco.trait.regression <- function(x, ...){
  summary(x)
}

#' \code{summary.eco.trait.regression} summarises an eco.trait.regression regression
#' @method summary eco.trait.regression
#' @S3method summary eco.trait.regression
#' @export
summary.eco.trait.regression <- function(x, ...){
  .summary.regression(x, "eco.trait.regression regression")
}

#' \code{plot.eco.trait.regression} plots an eco.trait.regression regression
#' @method plot eco.trait.regression
#' @S3method plot eco.trait.regression
#' @export
plot.eco.trait.regression <- function(x, ...){
  #Beware calling without knowing which environmental distance matrix we should be using
  if(!x$altogether) stop("Cannot call 'plot' on an element of an eco.trait.regression.list - uses plot(list, no.trait) instead")
  eco.matrix <- as.numeric(comm.dist(x$data$comm))
  trait.matrix <- as.numeric(traits.dist(x$data, TRUE))
  .plot.regression(trait.matrix, eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                   xlab="Trait Distance", ylab="Ecological Co-existnce", ...)
}

#' Print basic information about a list of eco.trait.regressions
#' @method print eco.trait.regression.list
#' @S3method print eco.trait.regression.list
#' @export
print.eco.trait.regression.list <- function(x, ...){
  cat("\neco.trait.regression list:\n")
  cat("Traits:\n")
  cat(colnames(x$data$traits), sep=", ")
  cat("\nTo examine each regression of each trait, use something like 'x[[1]]', or 'print(x[[2]])'\n")
  cat("To display all at once, call something like 'summary(regression.list)'\n")
}

#' Summarise an eco.trait.regression regression
#' @method summary eco.trait.regression.list
#' @export
summary.eco.trait.regression.list <- function(x, ...){
  cat("\neco.trait.regression list:\n")
  for(i in seq(ncol(x$data$traits))){
    cat("\n\n**", names(x$data$traits)[i], "**\n")
    summary(x[[i]], ...)
  }
}

#' Plot an eco.env.regression regression
#' @method plot eco.trait.regression.list
#' @export
plot.eco.trait.regression.list <- function(x, which=NULL, ...){
  eco.matrix <- as.numeric(comm.dist(x$data$comm))
  trait.matrix <- traits.dist(x$data, FALSE)
  if(is.null(which)){
    for(i in seq(ncol(x$data$traits)))
      .plot.regression(as.numeric(as.dist(trait.matrix[,,i])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                       xlab="Trait Distance", ylab="Ecological Co-existnce", main=names(x$data$traits)[i], ...)
  } else .plot.regression(as.numeric(as.dist(trait.matrix[,,which])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                          xlab="Trait Distance", ylab="Ecological Co-existnce", ...)
}