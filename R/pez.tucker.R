#' Calculate Tucker et al. metrics: (phylogenetic) richness,
#' regularity, and divergence
#'
#' As described in Tucker et al. (2016), richness metrics measure how
#' much diversity there is ("sum of"), divergence metrics how
#' different species are ("mean of"), and regularity how variable
#' species are ("variance of"). Note that these are only the *alpha*
#' diversity metrics of Tucker et al.; for the beta-diversity metrics
#' please see \code{\link{pez.dissimilarity}}.
#'
#' Most of these metrics do not involve comparison with some kind of
#' evolutionary-derived expectation for phylogenetic shape. Those that
#' do, however, such as PSV or Colless' index, make no sense unless
#' applied to a phylogenetic distance matrix - their null expectation
#' *requires* it. Using square-rooted distance matrices, or distance
#' matrices that incorporate trait information, can be an excellent
#' thing to do, but (for the above reasons), \code{pez} won't give you
#' an answer for metrics for which WDP thinks it makes no
#' sense. \code{pd} can (...up to you whether you should!...) be used
#' with a square-rooted distance matrix, but the results *will always
#' be wrong* if you do not have an ultrametric tree (branch lengths
#' proportional to time) and you will be warned about this. WDP
#' strongly feels you should only be using ultrametric phylogenies in
#' any case, but code to fix this bug is welcome.
#' @param data \code{\link{comparative.comm}} object
#' @param abundance.weighted Whether to use abundance-weighted forms of these
#'     metrics (default: FALSE).
#' @param which.eigen The eigen vector to calculate for the PhyloEigen
#'     metric (\code{eigen.sum})
#' @param sqrt.phy If TRUE (default is FALSE) your phylogenetic
#'     distance matrix will be square-rooted; specifying TRUE will
#'     force the square-root transformation on phylogenetic distance
#'     matrices (in the spirit of Leitten and Cornwell, 2014). See
#'     `details' for details about different metric calculations when
#'     a distance matrix is used.
#' @param traitgram If not NULL (default), a number to be passed to
#'     \code{funct.phylo.dist} (\code{phyloWeight}; the `a'
#'     parameter), causing analysis on a distance matrix reflecting
#'     both traits and phylogeny (0 --> only phylogeny, 1 --> only
#'     traits; see \code{funct.phylo.dist}). If a vector of numbers is
#'     given, \code{pez.shape} iterates across them and returns a
#'     \code{data.frame} with coefficients from each iteration. See
#'     `details' for details about different metric calculations when
#'     a distance matrix is used.
#' @param traitgram.p A value for `p' to be used in conjunction with
#'     \code{traitgram} when calling \code{funct.phylo.dist}.
#' @param ext.dist Supply an external species-level distance matrix
#'     for use in calculations. See `details' for comments on the use
#'     of distance matrices in different metric calculations.
#' @param q value for \emph{q} in \code{scheiner} (default 0.0001)
#' @note As mentioned above, \code{dist.fd} is calculated using a
#'     phylogenetic distance matrix if no trait data are available, or
#'     if you specify \code{sqrt.phy}. It is not calculated by default
#'     because it generates warning messsages (which WDP is loathe to
#'     suppress) which are related to the general tendency for a low
#'     rank of phylogenetic distance matrices. Much ink has been
#'     written about this, and in par this problem is why the
#'     \code{eigen.sum} measure came to be suggested.
#'
#' Many of these metrics will cause (inconsequential) warnings
#'     if given assemblages with only one species in them, and return
#'     NA/NaN values depending on the metric. I consider these
#'     `features', not bugs.
#' @return \code{data.frame}s with values for each metric.
#' @author Will Pearse
#' @seealso \code{\link{pez.shape}} \code{\link{pez.evenness}}
#'     \code{\link{pez.dispersion}} \code{\link{pez.dissimilarity}}
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus
#'     M.R. (2014). Metrics and Models of Community Phylogenetics. In:
#'     Modern Phylogenetic Comparative Methods and Their Application
#'     in Evolutionary Biology. Springer Berlin Heidelberg,
#'     pp. 451-464.
#' @references \code{PSV,PSR,PSE} Helmus M.R., Bland T.J., Williams
#'     C.K. & Ives A.R. (2007). Phylogenetic measures of
#'     biodiversity. American Naturalist, 169, E68-E83.
#' @references \code{PD} Faith D.P. (1992). Conservation evaluation
#'     and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references \code{colless} Colless D.H. (1982). Review of
#'     phylogenetics: the theory and practice of phylogenetic
#'     systematics. Systematic Zoology, 31, 100-104.
#' @references \code{gamma} Pybus O.G. & Harvey P.H. (2000) Testing
#'     macro-evolutionary models using incomplete molecular
#'     phylogenies. _Proceedings of the Royal Society of
#'     London. Series B. Biological Sciences 267: 2267--2272.
#' @references \code{taxon} Clarke K.R. & Warwick R.M. (1998). A
#'     taxonomic distinctness index and its statistical
#'     properties. J. Appl. Ecol., 35, 523-531.
#' @references \code{eed,hed} (i.e., \emph{Eed, Hed}) Cadotte M.W.,
#'     Davies T.J., Regetz J., Kembel S.W., Cleland E. & Oakley
#'     T.H. (2010). Phylogenetic diversity metrics for ecological
#'     communities: integrating species richness, abundance and
#'     evolutionary history. Ecology Letters, 13, 96-105.
#' @references \code{innd,mipd} Ness J.H., Rollinson E.J. & Whitney
#'     K.D. (2011). Phylogenetic distance can predict susceptibility
#'     to attack by natural enemies. Oikos, 120, 1327-1334.
#' @references \code{scheiner} Scheiner, S.M. (20120). A metric of
#'     biodiversity that integrates abundance, phylogeny, and
#'     function.  Oikos, 121, 1191-1202.
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' (output <- pez.richness(data))
#' (output <- pez.divergence(data))
#' (output <- pez.regularity(data))
#' @rdname pez.tucker
#' @name pez.tucker
#' @importFrom ape is.binary.tree
#' @export
pez.richness <- function(data, abundance.weighted=FALSE, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, which.eigen=1, q=0.0001)
{
    ## Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
        stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
    if(!is.null(traitgram)){
        if(length(traitgram) > 1){
            output <- vector("list", length(traitgram))
            for(i in seq_along(output))
                output[[i]] <- cbind(Recall(data, sqrt.phy, traitgram=traitgram[i], which.eigen=which.eigen, traitgram.p=traitgram.p, q=q), traitgram[i], sites(data))
            output <- do.call(rbind, output)
            names(output)[ncol(output)-1] <- "traitgram"
            names(output)[ncol(output)] <- "site"
            rownames(output) <- NULL
            return(output)
        } else {
            dist <- as.matrix(funct.phylo.dist(data, traitgram, traitgram.p))
            traitgram <- TRUE
        }      
    } else traitgram <- FALSE
    
    if(!is.null(ext.dist)){
        dist <- .check.ext.dist(ext.dist, species(data), ncol(data$comm))
        ext.dist <- TRUE
    } else ext.dist <- FALSE
    
    ## Setup
    if(!abundance.weighted)
        data$comm[data$comm > 0] <- 1
    coefs <- data.frame(row.names=rownames(data$comm))
    if(sqrt.phy)
        data <- .sqrt.phy(data)
    if(traitgram==FALSE & ext.dist==FALSE)
        dist <- cophenetic(data$phy)
    
    ## Filter metrics according to suitability and calculate
    functions <- setNames(c(.pd,.psr), c("pd","psr"))
    if(sqrt.phy == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "colless", "gamma", "eed", "hed")]
    if(traitgram == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "pd", "colless", "gamma", "eed", "hed", "scheiner")]
    if(ext.dist == TRUE)
        functions <- functions[!names(functions) %in% c("pd", "psv", "psr", "colless", "gamma", "eed", "hed", "scheiner")]
    if(!is.binary.tree(data$phy) & "colless" %in% names(functions))
        warning("Cannot compute Colless' index with non-binary tree")
    output <- lapply(functions, function(x) try(x(data, dist=dist, abundance.weighted=abundance.weighted, which.eigen=which.eigen), silent=TRUE))
    
    ## Clean up output and return
    output <- Filter(function(x) !inherits(x, "try-error"), output)
    output <- do.call(cbind, output)
    return(as.data.frame(output))
}

#' @rdname pez.tucker
#' @name pez.tucker
#' @export
pez.divergence <- function(data, abundance.weighted=FALSE, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, which.eigen=1, q=0.0001)
{
    ## Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
        stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
    if(!is.null(traitgram)){
        if(length(traitgram) > 1){
            output <- vector("list", length(traitgram))
            for(i in seq_along(output))
                output[[i]] <- cbind(Recall(data, sqrt.phy, traitgram=traitgram[i], which.eigen=which.eigen, traitgram.p=traitgram.p, q=q), traitgram[i], sites(data))
            output <- do.call(rbind, output)
            names(output)[ncol(output)-1] <- "traitgram"
            names(output)[ncol(output)] <- "site"
            rownames(output) <- NULL
            return(output)
        } else {
            dist <- as.matrix(funct.phylo.dist(data, traitgram, traitgram.p))
            traitgram <- TRUE
        }      
    } else traitgram <- FALSE
    
    if(!is.null(ext.dist)){
        dist <- .check.ext.dist(ext.dist, species(data), ncol(data$comm))
        ext.dist <- TRUE
    } else ext.dist <- FALSE
    
    ## Setup
    if(!abundance.weighted)
        data$comm[data$comm > 0] <- 1
    coefs <- data.frame(row.names=rownames(data$comm))
    if(sqrt.phy)
        data <- .sqrt.phy(data)
    if(traitgram==FALSE & ext.dist==FALSE)
        dist <- cophenetic(data$phy)
    
    ## Filter metrics according to suitability and calculate
    functions <- setNames(c(.scheiner,.mpd,.mntd,.psv), c(".scheiner","mpd","mntd","psv"))
    if(abundance.weighted)
        functions <- append(functions, setNames(.pse,"pse"))
    if(sqrt.phy == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "colless", "gamma", "eed", "hed")]
    if(traitgram == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "pd", "colless", "gamma", "eed", "hed", "scheiner")]
    if(ext.dist == TRUE)
        functions <- functions[!names(functions) %in% c("pd", "psv", "psr", "colless", "gamma", "eed", "hed", "scheiner")]
    if(!is.binary.tree(data$phy) & "colless" %in% names(functions))
        warning("Cannot compute Colless' index with non-binary tree")
    output <- lapply(functions, function(x) try(x(data, dist=dist, abundance.weighted=abundance.weighted, which.eigen=which.eigen), silent=TRUE))
    
    ## Clean up output and return
    output <- Filter(function(x) !inherits(x, "try-error"), output)
    output <- do.call(cbind, output)
    return(as.data.frame(output))
}

#' @rdname pez.tucker
#' @name pez.tucker
#' @export
pez.regularity <- function(data, abundance.weighted=FALSE, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, which.eigen=1, q=0.0001)
{
    ## Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
        stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
    if(!is.null(traitgram)){
        if(length(traitgram) > 1){
            output <- vector("list", length(traitgram))
            for(i in seq_along(output))
                output[[i]] <- cbind(Recall(data, sqrt.phy, traitgram=traitgram[i], which.eigen=which.eigen, traitgram.p=traitgram.p, q=q), traitgram[i], sites(data))
            output <- do.call(rbind, output)
            names(output)[ncol(output)-1] <- "traitgram"
            names(output)[ncol(output)] <- "site"
            rownames(output) <- NULL
            return(output)
        } else {
            dist <- as.matrix(funct.phylo.dist(data, traitgram, traitgram.p))
            traitgram <- TRUE
        }      
    } else traitgram <- FALSE
    
    if(!is.null(ext.dist)){
        dist <- .check.ext.dist(ext.dist, species(data), ncol(data$comm))
        ext.dist <- TRUE
    } else ext.dist <- FALSE
    
    ## Setup
    if(!abundance.weighted)
        data$comm[data$comm > 0] <- 1
    coefs <- data.frame(row.names=rownames(data$comm))
    if(sqrt.phy)
        data <- .sqrt.phy(data)
    if(traitgram==FALSE & ext.dist==FALSE)
        dist <- cophenetic(data$phy)
    
    ## Filter metrics according to suitability and calculate
    functions <- setNames(c(.colless,.gamma,.iac,.vpd,.taxon,.vntd,.eed,.hed,.haed), c(".colless","gamma","iac","vpd","taxon","vntd","eed","hed","haed"))
    if(sqrt.phy == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "colless", "gamma", "eed", "hed")]
    if(traitgram == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "pd", "colless", "gamma", "eed", "hed", "scheiner")]
    if(ext.dist == TRUE)
        functions <- functions[!names(functions) %in% c("pd", "psv", "psr", "colless", "gamma", "eed", "hed", "scheiner")]
    if(!is.binary.tree(data$phy) & "colless" %in% names(functions))
        warning("Cannot compute Colless' index with non-binary tree")
    output <- lapply(functions, function(x) try(x(data, dist=dist, abundance.weighted=abundance.weighted, which.eigen=which.eigen), silent=TRUE))
    
    ## Clean up output and return
    output <- Filter(function(x) !inherits(x, "try-error"), output)
    output <- do.call(cbind, output)
    return(as.data.frame(output))
}
