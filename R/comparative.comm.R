#' \code{comparative.comm} creates a community comparative ecology object
#' 
#' @param phy phylogeny (in ape::phylo format) of species
#' @param comm a community matrix (as used in vegan) with species as columns and rows as community samples
#' @param traits a data.frame of species traits (with row names matching those of 'comm')
#' @param env a data.frame of environmental data (all continuous
#' variables) with row names matching those of 'comm'
#' @param warn whether to warn if species/sites are dropped when
#' creating object (default: TRUE)
#' @param vcv whether to calculate variance-covariance matrix and
#' store in object to save time (but perhaps not memory) for later
#' calculations (default: FALSE). Storing a VCV now may cause problems
#' if you intend to use this object with the package \code{caper}.
#' @param force.root if phylogeny is unrooted, a root.edge of value
#' force.root will be added (default: -1, which means this will never
#' happen). Rarely needed, even more rarely advisable.
#' @details Basic checking of whether the input data match up is
#' performed; you need only supply 'comm' and 'phy', nothing else is
#' mandatory.  Subsetting according to communities (rows) and species
#' (columns) is possible - see 'examples' for details.
#' @return comparative.comm object
#' @examples \dontrun{
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
#' data[1:3,]
#' data[,1:3]
#' }
#' @importFrom ape is.rooted cophenetic.phylo
#' @importFrom ade4 scalewt
#' @export
comparative.comm <- function(phy, comm, traits=NULL, env=NULL, warn=TRUE, vcv=FALSE, force.root=-1){
  #Assertions and argument handling
  names <- list()
  #Phylogeny
  if(!inherits(phy, "phylo")) 
    stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
  if(! is.rooted(phy)){
    if(force.root != -1){
      phy$root.edge <- force.root
    } else {
      stop("'", deparse(substitute(phy)), "' is not rooted or has a basal polytomy.")
    }
  }
  names$phy <- deparse(substitute(phy))
  if(any(duplicated(phy$tip.label))) stop('Duplicate tip labels present in phylogeny')
  #Community matrix
  if(! is.matrix(comm)) stop("'", deparse(substitute(comm)), "' must be an object of class 'matrix'.")
  if(is.null(colnames(comm))) stop("'", deparse(substitute(comm)), "' must have column names (preferably species!)")
  if(is.null(rownames(comm))) stop("'", deparse(substitute(comm)), "' must have row names")
  if(any(is.na(comm))) stop("'", deparse(substitute(comm)), "' contains NAs")
  names$comm <- deparse(substitute(comm))
  #Traits
  if(!is.null(traits)){
    if(! is.data.frame(traits)) stop("'", deparse(substitute(traits)), "' must be an object of class 'data.frame'.")
    if(is.null(colnames(traits))) stop("'", deparse(substitute(traits)), "' must have row names (preferably species!)")
    if(is.null(rownames(traits))) stop("'", deparse(substitute(traits)), "' must have row names")
    names$data <- deparse(substitute(traits))
  } else names$data <- NULL
  #Environment
  if(!is.null(env)){
    if(! is.data.frame(env)) stop("'", deparse(substitute(env)), "' must be an object of class 'data.frame'.")
    if(is.null(colnames(env))) stop("'", deparse(substitute(env)), "' must have row names (preferably sites!)")
    if(is.null(rownames(env))) stop("'", deparse(substitute(env)), "' must have row names")
    names$env <- deparse(substitute(env))
  } else names$env <- NULL
  
  #Create intersection/drop lists and warn
  species.to.keep <- intersect(phy$tip.label, colnames(comm))
  if(!is.null(traits)) species.to.keep <- intersect(species.to.keep, rownames(traits))
  if(!is.null(env)) sites.to.keep <- intersect(rownames(comm), rownames(env)) else sites.to.keep <- rownames(comm)
  if(warn){
    if(length(setdiff(phy$tip.label, species.to.keep)) > 0)
      warning("Mismatch between phylogeny and other data, dropping ", length(setdiff(phy$tip.label, species.to.keep)), " tips")
    if(length(setdiff(colnames(comm), species.to.keep)) > 0)
      warning("Mismatch between community matrix and other data, dropping ", length(setdiff(colnames(comm), species.to.keep)), " columns")
    if(!is.null(traits) & length(setdiff(rownames(traits), species.to.keep)) > 0)
      warning("Mismatch between traits and other data, dropping ", length(setdiff(rownames(traits), species.to.keep)), " columns")
    if(length(setdiff(rownames(comm), sites.to.keep)) > 0)
      warning("Mismatch between community matrix and data, dropping ", length(setdiff(rownames(comm), sites.to.keep)), " rows")
    if(!is.null(env) & length(setdiff(rownames(env), sites.to.keep)) > 0)
      warning("Mismatch between env. data and other data, dropping ", length(setdiff(rownames(env), sites.to.keep)), " rows")
  }
  
  #Subset data and keep record
  # - make sure to re-order factor levels where necessary
  comm.sp.lost <- setdiff(colnames(comm), species.to.keep)
  comm <- comm[, colnames(comm) %in% species.to.keep]
  comm.sites.lost <- setdiff(rownames(comm), sites.to.keep)
  comm <- comm[rownames(comm) %in% sites.to.keep, ]
  if(any(dim(comm)==0))
    stop("ERROR: community data has no sites in common with rest of data")
  phy.sp.lost <- setdiff(phy$tip.label, species.to.keep)
  if(length(phy.sp.lost) == length(phy$tip.label))
    stop("ERROR: phylogeny has no species in common with rest of data")
  phy <- drop_tip(phy, phy.sp.lost)
  if(!is.null(traits)){
    traits.sp.lost <- setdiff(rownames(traits), species.to.keep)
    traits <- traits[rownames(traits) %in% species.to.keep, , drop = FALSE]
    if(any(dim(traits)==0))
      stop("ERROR: trait data has no species in common with rest of data")
  } else traits.sp.lost <- character(0)
  if(!is.null(env)){
    env.sites.lost <- setdiff(rownames(env), sites.to.keep)
    env <- env[rownames(env) %in% sites.to.keep, , drop = FALSE]
    if(any(dim(env)==0))
      stop("ERROR: environmental data has no sites in common with rest of data")
  } else env.sites.lost <- character(0)
  
  #Put species and sites in same order(s)
  # - leave phy alone (if caper needs it altered it will alter it)
  # - match everything to phylogeny's order (makes later work easier)
  comm <- comm[order(rownames(comm)), ]
  comm <- comm[, match(phy$tip.label, colnames(comm))]
  if(!is.null(traits)){
    traits <- traits[rownames(traits), , drop = FALSE]
    traits <- traits[, colnames(traits), drop = FALSE]
  }
  if(!is.null(env)){
    env <- env[rownames(env), , drop = FALSE]
    env <- env[, colnames(env), drop = FALSE]
  }
  
  #Handle VCV, make output, and return
  output <- list(phy=phy, comm=comm, data=traits, env=env, dropped=list(comm.sp.lost = comm.sp.lost,
                                                                 comm.sites.lost=comm.sites.lost,
                                                                 phy.sp.lost=phy.sp.lost,
                                                                 traits.sp.lost=traits.sp.lost,
                                                                 env.sites.lost=env.sites.lost), names=names)
  if(vcv) output$vcv <- cophenetic(phy) else output$vcv <- NULL
  class(output) <- c("comparative.comm", "comparative.data")
  return(output)
}

############################
# some useful generics######
############################
#' @param x \code{comparative.comm} object to be printed
#' @param ... not currently used
#' @rdname comparative.comm
#' @method print comparative.comm
#' @export
print.comparative.comm <- function(x, ...){
    #Argument checking
    if(!inherits(x, "comparative.comm"))
        stop("'", substitute(deparse(x)), "' not of class 'comparative.comm'")
    
    # basic summary data
    cat("Comparative community dataset of", ncol(x$comm), "taxa:\n")
    cat("Phylogeny:", x$names$phy, "\n")
    cat("   ", x$phy$Nnode, " internal nodes", sep='')
    if(!is.null(x$vcv)) cat(', VCV matrix present\n') else cat("\n")
    cat("Community data:", x$names$comm, "\n")
    cat("    ", nrow(x$comm), " sites, ", ncol(x$comm), " taxa\n")
    
    if(!is.null(x$data)){
	    cat("Trait data:", x$names$data, "\n")
    	cat("   ", ncol(x$data), " variables\n")
    } else cat("Trait data: None\n")
    
    if(!is.null(x$env)){
    	cat("Environmental data:", x$names$env, "\n")
    	cat("   ", nrow(x$env), " sites, ", ncol(x$env), " variables\n")
    } else cat("Environmental data: None\n")
  
	# report on mismatch on merge
    if(length(x$dropped$phy.sp.lost)>0 | length(x$dropped$comm.sp.lost)>0 | length(x$dropped$data.sp.lost)>0){
  	 cat('Dropped taxa:\n')
  	 cat('   ', x$names$phy , ' : ', length(x$dropped$phy.sp.lost), '\n')
  	 cat('   ', x$names$comm , ' : ', length(x$dropped$comm.sp.lost), '\n')
     if(!is.null(x$data))
        cat('   ', x$names$data , ' : ', length(x$dropped$data.sp.lost), '\n')
    }
  if(length(x$dropped$comm.sites.lost)>0 | length(x$dropped$env.sites.lost)>0){
  	 cat('Dropped sites:\n')
       cat('   ', x$names$comm , ' : ', length(x$dropped$comm.sites.lost), '\n')
  	 if(!is.null(x$env))
       cat('   ', x$names$env , ' : ', length(x$dropped$env.sites.lost), '\n')
  }
}

#' Slice out the species or site that you want from your comparative.comm object
#' @param x \code{comparative.comm} object to be subset
#' @param sites numbers of sites to be kept or dropped from \code{x};
#' cannot be given as names, but rather numbers. For example,
#' \code{x[1:5,]}, or \code{x[-1:-5,]}, but not \code{x[c("site a",
#' "site b"),]}.
#' @param spp numbers of species to be kept or dropped from
#' \code{x}; cannot be given as names, but rather numbers. For
#' example, \code{x[,1:5]}, or \code{x[,-1:-5]}, but not
#' \code{x[c("sp a", "sp b"),]}.
#' @param warn whether to warn if species/sites are dropped when
#' creating object (default: TRUE)
#' @export
"[.comparative.comm" <- function(x, sites, spp, warn=FALSE) {
  #Assertions and setup
  if(!inherits(x, "comparative.comm"))
    stop("'", substitute(deparse(x)), "' not of class 'comparative.comm'")
    
  
  #Handle species
  if(!missing(spp)){
    if(is.null(spp)) stop("Null indices not permitted on comparative community data objects")
    spp.to.keep <- colnames(x$comm)[spp]
    comm <- x$comm[, spp.to.keep]
    phy <- drop_tip(x$phy, setdiff(x$phy$tip.label, spp.to.keep))
    if(!is.null(x$data))
      traits <- x$data[spp.to.keep, ] else traits <- NULL
    new.x <- comparative.comm(phy, comm, traits, x$env, warn=warn)
  } else new.x <- x
  
  #Handle sites
  if(!missing(sites)){
    if(is.null(sites)) stop("Null indices not permitted on comparative community data objects")
    sites.to.keep <- rownames(new.x$comm)[sites]
    comm <- new.x$comm[sites.to.keep, ]
    if(!is.null(x$env))
      env <- new.x$env[sites.to.keep, ] else env <- new.x$env
    new.x <- comparative.comm(x$phy, comm, new.x$data, env, warn=FALSE)
  }
  
  #Warn of dropped species/sites (if asked)
  if(warn){
    orig.species <- colnames(x$comm)
    orig.sites <- rownames(x$comm)
    if(length(setdiff(orig.species, new.x$phy$tip.label)) > 0)
      warning("Mismatch between phylogeny and other data, dropping ", length(setdiff(orig.species, new.x$phy$tip.label)), " tips")
    if(length(setdiff(orig.sites, colnames(new.x$comm))) > 0)
      warning("Mismatch between community matrix and other data, dropping ", length(setdiff(orig.species, colnames(new.x$comm))), " columns")
    if(!is.null(new.x$data) & length(setdiff(orig.species, rownames(new.x$data))) > 0)
      warning("Mismatch between traits and other data, dropping ", length(setdiff(orig.species, rownames(new.x$data))), " columns")
    if(length(setdiff(rownames(new.x$comm), orig.sites)) > 0)
      warning("Mismatch between community matrix and data, dropping ", length(setdiff(orig.sites, rownames(new.x$comm))), " rows")
    if(!is.null(new.x$env) & length(setdiff(rownames(new.x$env), orig.sites)) > 0)
      warning("Mismatch between env. data and other data, dropping ", length(setdiff(orig.sites, rownames(new.x$env))), " rows")
  }
	
  #Return
	return(new.x)
}
