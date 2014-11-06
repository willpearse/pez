#Null models for community data
#Internal use only; mostly based on picante's code
#' @importFrom picante randomizeMatrix
.eco.null <- function(comm, method=c("taxa.labels", "richness", "frequency", "independentswap", "trialswap"), swap.iter=1000){
  #Checks
  method <- match.arg(method)
  
  if(method == "taxa.labels")
    curr.rnd <- comm[,sample(seq(ncol(comm)))] else{
        curr.rnd <- randomizeMatrix(comm, method)
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
  if(!is.null(class)) output$type <- class
  class(output) <- "eco.xxx.regression"
  return(output)
}

#Printing summaries of regressions
#ADD RETURN STATEMENTS FOR FALL-THROUGH LIKE OTHERS
#' @method summary eco.xxx.regression
#' @param object \code{eco.xxx.regression} object
#' @export
#' @rdname eco.xxx.regression
#' @importFrom quantreg summary.rq print.rq
summary.eco.xxx.regression <- function(object, ...){
    x <- object  ## FIXME:  lazy hack
    cat("\n", x$type, "\n", sep="")
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
        print(x$observed)
    if(x$method == "quantile")
        print.rq(x$observed)
    if(x$method == "lm")
        print(summary(x$observed))
    cat("\n")
}

#' @method print eco.xxx.regression
#' @param x \code{eco.xxx.regression} object
#' @export
#' @rdname eco.xxx.regression
print.eco.xxx.regression <- function(x, ...){
    summary(x, ...)
}

.plot.regression <- function(x, y, observed, randomisations,
                             method=c("quantile", "lm", "mantel"),
                             permute=0, ...){
    method <- match.arg(method)
    plot(y ~ x, ...)
    if(method == "lm"){
        abline(observed, lwd=3)
                                        # Easiest way to silence
                                        # lapply...
        if(permute>0 && method=="lm")
            silent<-lapply(randomisations, abline, col="red")
    }
    if(method == "quantile"){
                                        # Check to see if we've got
                                        # more than one tau value and
                                        # plot accordingly
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

#' @method plot eco.xxx.regression
#' @export
#' @rdname eco.xxx.regression
plot.eco.xxx.regression <- function(x, ...){
    if(x$type == "eco.env.regression"){
        #Beware calling without knowing which environmental distance matrix we should be using
        if(!x$altogether) stop("Cannot call 'plot' on an element of an eco.env.regression.list - uses plot(list, no.trait) instead")
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        env.matrix <- as.numeric(pianka.dist(x$data, TRUE))
        .plot.regression(env.matrix, eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                         xlab="Pianka's Distance", ylab="Ecological Co-existnce", ...)
    }
    
    if(x$type == "eco.trait.regression"){
        #Beware calling without knowing which environmental distance matrix we should be using
        if(!x$altogether) stop("Cannot call 'plot' on an element of an eco.trait.regression.list - uses plot(list, no.trait) instead")
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        trait.matrix <- as.numeric(traits.dist(x$data, TRUE))
        .plot.regression(trait.matrix, eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                         xlab="Trait Distance", ylab="Ecological Co-existnce", ...)
    }

    if(x$type == "eco.phy.regression"){
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        phy.matrix <- as.numeric(as.dist(cophenetic(x$data$phy)))
        .plot.regression(phy.matrix, eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                         xlab="Phylogenetic Distance", ylab="Ecological Co-existnce", ...)
    }
}

#' List of eco.xxx.regressions
#' @method summary eco.xxx.regression.list
#' @param object \code{eco.xxx.regression.list} object
#' @param ... additional arguments to plotting functions
#' @export
#' @rdname eco.xxx.regression.list
#' @name eco.xxx.regression.list
summary.eco.xxx.regression.list <- function(object, ...){
    x <- object ## FIXME:  lazy hack
    if(x$type == "eco.env.regression.list"){
        cat("\neco.env.regression list:\n")
        for(i in seq(ncol(object$data$env))){
            cat("\n\n**", names(object$data$env)[i], "**\n")
            summary(object[[i]], ...)
        }
        return()
    }

    if(x$type == "eco.trait.regression.list"){
        cat("\neco.trait.regression list:\n")
        for(i in seq(ncol(object$data$data))){
            cat("\n\n**", names(object$data$data)[i], "**\n")
            summary(x[[i]], ...)
        }
        return()
    }

}

#' @method print eco.xxx.regression.list
#' @param x \code{eco.xxx.regression.list} object
#' @export
#' @rdname eco.xxx.regression.list
print.eco.xxx.regression.list <- function(x, ...){
    if(x$type == "eco.env.regression.list"){
        cat("\neco.env.regression list:\n")
        cat("Traits:\n")
        cat(colnames(x$data$env), sep=", ")
        cat("\nTo examine each regression of each trait, use something like 'x[[1]]', or 'print(x[[2]])'\n")
        cat("To display all at once, call something like 'summary(regression.list)'\n")
        return()
    }
    if(x$type == "eco.trait.regression.list"){
        cat("\neco.trait.regression list:\n")
        cat("Traits:\n")
        cat(colnames(x$data$data), sep=", ")
        cat("\nTo examine each regression of each trait, use something like 'x[[1]]', or 'print(x[[2]])'\n")
        cat("To display all at once, call something like 'summary(regression.list)'\n")
        return()
    }
}

#' @method plot eco.xxx.regression.list
#' @export
#' @rdname eco.xxx.regression.list
plot.eco.xxx.regression.list <- function(x, ...){
    if(x$type == "eco.env.regression.list"){
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        env.matrix <- pianka.dist(x$data, FALSE)
        if(is.null(which)){
            for(i in seq(ncol(x$data$env)))
                .plot.regression(as.numeric(as.dist(env.matrix[,,i])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                                 xlab="Pianka's Distance", ylab="Ecological Co-existnce", main=names(x$data$env)[i], ...)
        } else .plot.regression(as.numeric(as.dist(env.matrix[,,which])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                                xlab="Pianka's Distance", ylab="Ecological Co-existnce", ...)
        return()
    }

    if(x$type == "eco.trait.regression.list"){
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        trait.matrix <- traits.dist(x$data, FALSE)
        if(is.null(which)){
            for(i in seq(ncol(x$data$data)))
                .plot.regression(as.numeric(as.dist(trait.matrix[,,i])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                                 xlab="Trait Distance", ylab="Ecological Co-existnce", main=names(x$data$data)[i], ...)
        } else .plot.regression(as.numeric(as.dist(trait.matrix[,,which])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                                xlab="Trait Distance", ylab="Ecological Co-existnce", ...)
        return()
    }
}

##' Trim a phylogeny
##'
##' This is a weak wrapper around \code{ape}'s
##' \code{\link{drop.tip}}. Importantly, if asked to drop no species
##' from a phylogeny, it will just return the phylogeny (not an empty
##' phylogeny, as \code{\link{drop.tip}}) will.
##'
##' @param tree An \code{\link[ape:phylo]{phylo}} object
##' @param spp A vector of species (one, many, or none) to be removed
##' from \code{tree}
##' @return \code{\link[ape:phylo]{phylo}} object
##' @seealso \code{\link[ape:drop.tip]{drop.tip}} \code{\link[ape:extract.clade]{extract.clade}}
##' @export
drop_tip <- function(tree, spp)
  if(length(setdiff(tree$tip.label, spp)) >0) return(drop.tip(tree, spp)) else return(tree)

#' @method summary phy.structure
#' @export
summary.phy.structure <- function(object, ...){
    #Argument checking
    if(!inherits(object, "phy.structure"))
        stop("'", substitute(deparse(object)), "' not of class 'phy.structure'")
    
    #Shape
    if(object$type == "shape"){
        cat("\nShape metrics in this object:\n")
        if(!is.null(object$psv))
            cat("\tPSV (psv)\n")
        if(!is.null(object$psr))
            cat("\tPSR (psr)\n")
        if(!is.null(object$mpd))
            cat("\tMean Phylogenetic Distance (mpd)\n")
        if(!is.null(object$mntd))
            cat("\tMean Mean Nearest Taxon Distance (mntd)\n")
        if(!is.null(object$pd))
            cat("\tPhylogenetic Distance (pd)\n")
        if(!is.null(object$pd.ivs))
            cat("\tPhylogenetic Distance stanardised according to number of species in sample (pd.ivs)\n")
        if(!is.null(object$colless))
            cat("\tColless' index (colless)\n")
        if(!is.null(object$gamma))
            cat("\tGamma (gamma)\n")
        if(!is.null(object$taxon))
            cat("\tTaxonomic diversity (taxon)\n")
        if(!is.null(object$eigen.sum))
            cat("\tSum of dominant eigenvector (eigen.sum)\n")
        if(!is.null(object$eed))
            cat("\tEED (eed)\n")
        if(!is.null(object$hed))
            cat("\tHED (hed)\n")
        if(!is.null(object$dist.fd))
            cat("\tFunctional trait distance (dist.fd)\n")
        cat("Use something like 'output$psv' to work with each measure\n")
        cat("...or coef(output) to get a table of metric values\n")
        return()
    }
    #Evenness
    if(object$type == "evenness"){
        cat("\nEvenness metrics in this object:\n")
        if(!is.null(object$pse))
            cat("\tPhylogenetic Species Evenness (pse)\n")
        if(!is.null(object$rao))
            cat("\tRao's quadratic entropy (rao)\n")
        if(!is.null(object$taxon))
            cat("\tTaxonomic diversity (taxon)\n")
        if(!is.null(object$entropy))
            cat("\tPhylogenetic entropy (entropy)\n")
        if(!is.null(object$pae))
            cat("\tPAC (pac)\n")
        if(!is.null(object$iac))
            cat("\tIAC (iac)\n")
        if(!is.null(object$Haed))
            cat("\tHaed (Haed)\n")
        if(!is.null(object$Eaed))
            cat("\tEaed (Eaed)\n")
        if(!is.null(object$pst))
            cat("\tPst - Simpson's diversity (pst)\n")
        if(!is.null(object$lambda))
            cat("\tLambda transformation (lambda)\n")
        if(!is.null(object$delta))
            cat("\tDelta transformation (delta)\n")
        if(!is.null(object$kappa))
            cat("\tKappa transformation (kappa)\n")
        if(!is.null(object$mpd))
            cat("\tMean Phylogenetic Distance (mpd)\n")
        if(!is.null(object$mntd))
            cat("\tMean Mean Nearest Taxon Distance (mntd)\n")
        if(!is.null(object$dist.fd))
            cat("\tFunctional trait distance (dist.fd)\n")
        cat("Use something like 'output$rao' to work with each measure\n")
        cat("...or coef(output) to get a table of metric values\n")
        return()
    }
    #Dispersion
    if(object$type == "dispersion"){
        cat("\nDispersion metrics:\n")
        if(!is.null(object$sesmpd)){
            cat("\t SESmpd (sesmpd)\n")
            print(object$sesmpd)
        }
        if(!is.null(object$sesmntd)){
            cat("\tSESmntd (sesmntd)\n")
            print(object$sesmntd)
        }
        if(!is.null(object$sespd)){
            cat("\tSESpd (sespd)\n")
            print(object$sespd)
        }
        if(!is.null(object$innd)){
            cat("\tINND (innd)\n")
            print(object$innd)
        }
        if(!is.null(object$d)){
            cat("\tD (d)\n")
            print(object$d)
        }
        cat("Use something like 'output$d' to work with each measure\n")
        cat("...or coef(output) to get a table of metric values\n")
        return()
    }
    #Dissimilarity
    if(object$type == "dissimilarity"){
        cat("\nDissimilarity metrics in this object:\n")
        if(!is.null(object$unifrac))
            cat("\tUniFrac (unifrac)\n")
        if(!is.null(object$pcd))
            cat("\tPCD (pcd)\n")
        if(!is.null(object$phylosor))
            cat("\tPhyloSor (phylosor)\n")
        if(!is.null(object$comdist))
            cat("\tComdist (comdist)\n")
        cat("Use something like 'output$unifrac' to work with each measure\n")
        return()
    }
    cat("No valid phylogenetic structure metrics found.\n")
}

#' @export
print.phy.structure <- function(x, ...){
    summary(x)
}

#' @method coef phy.structure
#' @export
coef.phy.structure <- function(object, ...){
    #Shape
    if(object$type == "dissimilarity")
      cat("\nCannot produce a simple summary of dissimilarity matrices. Sorry.\n")
    else
      return(object$coefs)
}

.removeErrors <- function(object) {
    if(inherits(object, "try-error")) return(NULL)
    object
}

#' Manipulate internals of \code{comparative.comm} object
#' @param expr expression to be evaluated within the scope of
#' \code{data}
#' @param ... ignored
#' @export
#' @rdname cc.manip
within.comparative.comm <- function(data, expr, ...){
	
	parent <- parent.frame()
	e <- evalq(environment(), data, parent)
	eval(substitute(expr), e)
	
	cc <- as.list(e)
	class(cc) <- class(data)
        return(cc)
}
