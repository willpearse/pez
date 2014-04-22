#' Run an eco.env.regression regression (Cavender-Bares et al. 2004)
#' 
#' \code{eco.env.regression} runs a regression of ecological community coexistence versus species' trait-distances' 
#' overlap based on environmental tolerances
#' 
#' @param data a comparative community ecology object on which to run the regression
#' @param randomisation what kind of null distributions to compare your data with - one of:
#' taxa.labels, richness, frequency, sample.pool, phylogeny.pool, independentswap, trialswap
#' (as implemented in 'picante')
#' @param permute the number of null permutations to perform
#' @param altogether use distance matrix based on all traits (default TRUE), or perform separate regressions for each trait
#' @details This is extremely unchcked, so beware!
#' @note Like the eco.trait and eco.env methods, this is a data-hungry method. Warnings will be generated if any of the methods cannot be fitted properly (the examples below give toy examples of this). In such cases the summary and plot methods of these functions may generate errors; use 'traceback()' to examine where these are coming from, and consider whether you want to be working with the data generating these errors. I am loathe to hide these errors or gloss over them, because they represent the reality of your data!
#' @author Will Pearse, Jeannine Cavender-Bares
#' @examples \dontrun{
#' data(phylocom)
#' dummy.env <- data.frame(letters=letters[1:6], repeated=rep("A", 6), row.names=rownames(phylocom$sample))
#' data <- comparative.comm(phylocom$phy, phylocom$sample, env=dummy.env)
#' eco.env.regression(data, permute=10)
#' plot(eco.env.regression(data, permute=10, method="lm))
#' plot(eco.env.regression(data, permute=10, method="lm", altogether=FALSE))
#' }
#' @import picante quantreg vegan
#' @export
eco.env.regression <- function(data,
                                 randomisation=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
                                 permute=0, n.clades=1, method=c("quantile", "lm", "mantel"), swap.iter=1000, altogether=TRUE, ...){
  #Assertions and argument handling
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  randomisation <- match.arg(randomisation)
  method <- match.arg(method)
  if(permute < 0) stop("Can't have negative null permutations!")
  
  #Setup matrices
  eco.matrix <- comm.dist(data$comm)
  if(!is.null(data$env)) traits.matrix <- pianka.dist(data, altogether) else stop("'data' must contain environmental data for an environmental regression!")
  
  #Observed eco.env.regression
  if(altogether){
    observed <- .eco.env.regression(eco.matrix, traits.matrix, NULL, method, ...)} else {
      #Do separately for all traits
      observed <- vector("list", ncol(data$env))
      for(i in seq(ncol(data$env)))
        observed[[i]] <- .eco.env.regression(eco.matrix, traits.matrix, i, method, ...)
    }
  
  #Randomisations
  if(altogether){
    #Using mean of traits
    randomisations <- vector(mode="list", length=permute)
    #This won't execute if permute is 0...
    for(i in seq(from=1, length.out=permute)){
      curr.rnd <- .eco.null(data$comm, randomisation)
      rnd.mat <- comm.dist(curr.rnd)
      if(any(is.na(rnd.mat))){
        warning("NAs in permuted community matrix; skipping this iteration")
        next()
      }
      randomisations[[i]] <- .eco.env.regression(rnd.mat, traits.matrix, NULL, method, ...)
    }
  } else {
    #Separately for each trait
    # - preallocate
    randomisations <- vector(mode="list", length=ncol(data$env))
    for(i in seq_along(randomisations)) randomisations[[i]] <- vector("list", permute)
    for(i in seq(from=1, length.out=permute)){
      curr.rnd <- .eco.null(data$comm, randomisation)
      rnd.mat <- comm.dist(curr.rnd)
      if(any(is.na(rnd.mat))){
        warning("NAs in permuted community matrix; skipping this iteration")
        next()
      }
      for(j in seq(ncol(data$env)))
        randomisations[[j]][[i]] <- .eco.env.regression(rnd.mat, traits.matrix, j, method, ...)
    }
  }
  
  #Prepare output (...and return)
  if(altogether){
    output <- .prepare.regression.output(observed, randomisations, method, permute, "eco.env.regression")
    output$altogether <- TRUE
  } else{
      output <- vector("list", ncol(data$env))
      for(i in seq_along(output)){
        output[[i]] <- .prepare.regression.output(observed[[i]], randomisations[[i]], method, permute, "eco.env.regression")
        output[[i]]$altogether <- FALSE
      }
      output$altogether <- FALSE
      output$type <- "eco.env.regression.list"
      class(output) <- "ecophyl.regression.list"
    }
  output$data <- data
  output$permute<-permute;output$method<-method
  return(output)
}

#Perform one set of EcoPhy regressions
.eco.env.regression <- function(eco.mat, trait.mat, which.trait=NULL, method=c("quantile", "lm", "mantel"), ...){
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
