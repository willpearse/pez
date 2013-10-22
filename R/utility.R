#Null models for community data
#Internal use only; mostly based on picante's code
.eco.null <- function(comm, method=c("taxa.labels", "richness", "frequency", "independentswap", "trialswap"), swap.iter=1000){
  #Checks
  method <- match.arg(method)
  
  if(method == "taxa.labels")
    curr.rnd <- comm[,sample(seq(ncol(comm)))]
  if(method == "richness"){
    curr.rnd <- .C("frequency", m = as.numeric(comm), as.integer(nrow(comm)), as.integer(ncol(comm)), PACKAGE = "picante")
    curr.rnd <- matrix(curr.rnd$m, nrow = nrow(comm), dimnames = list(rownames(comm), colnames(comm)))
  }
  if(method == "frequency"){
    curr.rnd <- .C("richness", m = as.numeric(comm), as.integer(nrow(comm)), as.integer(ncol(comm)), PACKAGE = "picante")
    curr.rnd <- matrix(curr.rnd$m, nrow = nrow(comm), dimnames = list(rownames(comm), colnames(comm)))
  }
  if(method == "independentswap"){
    curr.rnd <- .C("independentswap", m = as.numeric(comm), as.integer(swap.iter), as.integer(nrow(comm)), as.integer(ncol(comm)), PACKAGE = "picante")
    curr.rnd <- matrix(curr.rnd$m, nrow = nrow(comm), dimnames = list(rownames(comm), colnames(comm)))
  }
  if(method == "trialswap"){
    curr.rnd <- .C("trialswap", m = as.numeric(comm), as.integer(swap.iter), as.integer(nrow(comm)), as.integer(ncol(comm)), PACKAGE = "picante")
    curr.rnd <- matrix(curr.rnd$m, nrow = nrow(comm), dimnames = list(rownames(comm), colnames(comm)))
  }
  return(curr.rnd)
}

#Prearing regressions for output
.prepare.regression.output <- function(observed, randomisations=NULL, method=c("lm", "quantile", "mantel"), permute=0, class=NULL){
  #Checks
  method <- match.arg(method)
  
  if(method == "lm"){
    if(permute > 0)
      rnd.slopes <- sapply(randomisations, function(x) coef(x)[2]) else rnd.slopes <- NULL
    obs.slope <- coef(observed)[2]
  }
  
  if(method == "quantile"){
    #Beware quantile regressions with only one tau!
    if(is.null(dim(coef(observed))))
      obs.slope <- coef(observed)[2] else obs.slope <- coef(observed)[2,]
    if(permute > 0){
      if(is.null(dim(coef(observed))))
        rnd.slopes <- sapply(randomisations, function(x) coef(x)[2]) else rnd.slopes <- sapply(randomisations, function(x) coef(x)[2,])
    } else rnd.slopes <- NULL
  }
  
  if(method == "mantel"){
    obs.slope <- observed$statistic
    if(permute > 0)
      rnd.slopes <- sapply(randomisations, function(x) x$statistic) else rnd.slopes <- NULL
  }
  
  output <- list(observed=observed, randomisations=randomisations, obs.slope=obs.slope, rnd.slopes=rnd.slopes,
                 method=method, permute=permute, randomisation=randomisations)
  if(!is.null(class)) class(output) <- class
  return(output)
}

#Printing summaries of regressions
.summary.regression <- function(x, header){
  cat("\n", header, "\n", sep="")
  cat("Method: ", x$method, "\n")
  if(x$permute > 0)
    cat("Randomisation: ", x$method, "; Permutations: ", x$permute, "\n") else cat("Randomisation: NONE\n")
  if(x$method == "quantile" && length(x$obs.slope)>1)
    cat("Observed slopes (at specified taus): ", paste(round(x$obs.slope,2), collapse=","), "\n") else cat("Observed slope: ", round(x$obs.slope,2), "\n")
  if(x$permute > 0 ){
    if(x$method == "quantile" && length(x$obs.slope)>1){
      cat("Random slope means (at specified taus) +/- SD:\n")
      for(i in seq(nrow(x$rnd.slopes)))
        cat(round(mean(x$rnd.slopes[i,]),2), " +/- ", round(sd(x$rnd.slopes[i,]),4), "\n")
    } else cat("Random slope mean +/-SD: ", round(mean(x$rnd.slopes),2), " +/- ", round(sd(x$rnd.slopes),4), "\n")
  }
  cat("Observed model summary:\n")
  if(x$method == "mantel")
    print(x$observed) else print(summary(x$observed))
  cat("\n")
}

#Plotting summaries of regressions
.plot.regression <- function(x, y, observed, randomisations, method=c("quantile", "lm", "mantel"), permute=0, ...){
  method <- match.arg(method)
  plot(y ~ x, ...)
  if(method == "lm"){
    abline(observed, lwd=3)
    #Easiest way to silence lapply...
    if(permute>0 && method=="lm")
      silent<-lapply(randomisations, abline, col="red")
  }
  if(method == "quantile"){
    #Check to see if we've got more than one tau value and plot accordingly
    if(is.null(dim(coef(observed)))){
      abline(coef(observed), lwd=3)
      for(j in seq(from=1,length.out=permute))
        abline(coef(randomisations[[j]]), col="red")
    } else {
      for(j in seq(ncol(coef(observed)))){
        abline(coef(observed)[,j], lwd=3)
        for(k in seq_along(randomisations))
          abline(coef(randomisations[[k]])[,j], col="red")
      }
    }
  }
}

#Trim a phylogeny (ape work-around)
#' @export
drop_tip <- function(tree, spp)
  if(length(setdiff(tree$tip.label, spp)) >0) return(drop.tip(tree, spp)) else return(tree)