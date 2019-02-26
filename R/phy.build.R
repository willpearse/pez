#' Build a novel phylogeny from existing data
#' 
#' @param tree \code{\link[ape:phylo]{phylo}} phylogeny to have those
#' species inserted into it
#' @param species vector of species names to be bound into the tree if
#' missing from it
#' @param split the character that splits genus and species names in
#' your phylogeny. Default is \code{_}, i.e. Quercus_robur.
#' @param ... ignored
#' @details \code{congeneric.merge} Binds missing species into a
#' phylogeny by replacing all members of the clade it belongs to with
#' a polytomy. Assumes the \code{tip.labels} represent Latin
#' binomials, split by the \code{split} argument. This code was
#' originally shipped with phyloGenerator - this is the \code{merge}
#' method in that program.
#' @note Thank you to Josep Padulles Cubino, who found that the genus
#'     name splitting in a previous version of this function could
#'     cause incorrect placement in oddly named cases. As with all
#'     phylogenetic construction tools, the output from
#'     \code{congeneric.merge} should be checked for accuracy before
#'     use.
#' @return \code{\link[ape:phylo]{phylo}} phylogeny
#' @author Will Pearse
#' @importFrom ape drop.tip
#' @rdname phy.build
#' @name phy.build
#' @references Pearse W.D. & Purvis A. phyloGenerator: an automated phylogeny generation tool for ecologists. Methods in Ecology and Evolution 4(7): 692--698.
#' @export
#' @examples
#' tree <- read.tree(text="((a_a:1,b_b:1):1, c_c:2):1;")
#' tree <- congeneric.merge(tree, c("a_nother", "a_gain", "b_sharp"))
congeneric.merge <- function(tree, species, split="_", ...){
    if(!inherits(tree, "phylo"))
        stop("'tree' must be a phylogeny")
    before <- sum(species %in% tree$tip.label)
    for(i in seq_along(species)){
        if(!is.null(species[i]) & !species[i] %in% tree$tip.label){
            genus <- paste0("^",strsplit(species[i], split, fixed=TRUE)[[1]][1],split)
            matches <- unique(grep(genus, tree$tip.label, value=TRUE))
            if(length(matches) > 0){
                tree <- drop.tip(tree, matches[-1])
                tip.length <- .find.unique.branch.length(tree, matches[1])
                polytomy <- .make.polytomy(unique(c(matches, species[i])), (tip.length/2))
                tree <- bind.replace(tree, polytomy, matches[1], tip.length)
            }
        }
    }
    cat("\nNumber of species in tree before:", before)
    cat("\nNumber of species in tree now:   ", sum(species %in% tree$tip.label), "\n")
    return(tree)
}

#' @param backbone backbone phylogeny (\code{\link[ape:phylo]{phylo}})
#' into which the donor is to be bound
#' @param donor phylogeny (\code{\link[ape:phylo]{phylo}}) to bound
#' into the backbone phylogeny
#' @param replacing.tip.label the species in the donor phylogeny
#' that's being replaced by the donor phylogeny
#' @param donor.length how deep the donor phylogeny should be cut into
#' the backbone phylogeny. If NA (default), then the bladj algorithm
#' is followed (or, in plain English, it's put half-way along the
#' branch)
#' @details \code{bind.replace} Binds a phylogeny (donor) into a
#' bigger phylogeny ('backbone'); useful if you're building a
#' phylogeny a la Phylomatic. A version of this R code was shipped
#' with phyloGenerator (Pearse & Purvis 2013). This is really an
#' internal function for \code{congeneric.merge}, but hopefully it's
#' of some use to you!
#' @return phylogeny (\code{\link[ape:phylo]{phylo}})
#' @author Will Pearse
#' @importFrom ape bind.tree
#' @rdname phy.build
bind.replace <- function(backbone, donor, replacing.tip.label, donor.length=NA){	
    bind.point <- which(backbone$tip.label == replacing.tip.label)
    backbone <- bind.tree(backbone, donor, where=bind.point)
    which.tip <- which(backbone$tip.label == donor$tip.label[1])
    which.node <- backbone$edge[which(backbone$edge[,2] == which.tip),1]
    which.edge <- which(backbone$edge[,2] == which.node)
    tip.length <- backbone$edge.length[which.edge]
    if(is.na(donor.length)){
        backbone$edge.length[which.edge] <- tip.length/2
    } else {
        backbone$edge.length[which.edge] <- tip.length - donor.length/2
    }
    return(backbone)	
}

#' @importFrom ape as.phylo.formula
.make.polytomy <- function(species, tip.length=NA){	
    d.f <- data.frame(spp=factor(species))	
    polytomy <- as.phylo.formula(~spp, data=d.f)	
    if(!is.na(tip.length)) polytomy$edge.length <- rep(tip.length, length(species))	
    return(polytomy)	
}

.find.unique.branch.length <- function(tree, tip){	
    which.tip <- which(tree$tip.label == tip)
    which.edge <- which(tree$edge[,2] == which.tip)
    tip.length <- tree$edge.length[which.edge]
    return(tip.length)	
}

#' \code{congeneric.impute} sequentially add species to a phylogeny to
#' form an _imputed_ bifurcating tree. Makes use of a result from
#' Steel & Mooers (2010) that gives the expected branch-length under a
#' Yule model whether the rate of diversification has been
#' estimated. The intention of this is to approximate the method by
#' which phylogenetic structure is sampled from the prior in BEAST;
#' i.e., to approximate the standard Kuhn et al. (2011) method for
#' imputing a phylogeny. When using \code{congeneric.impute} you
#' should (1) repeat your analyses across many (if in doubt,
#' thousands) of separate runs (see Kuhn et al. 2011) and (2) check
#' for yourself that your trees are unbiased for your purpose - I make
#' no guarantee this is appropriate, and in many cases I think it
#' would not be. See also 'notes' below.
#' 
#' @note Caveats for \code{congeneric.impute}: something I noticed is
#'     that BEAST randomly picks an edge to break when adding species
#'     (starting from a null tree), and this is the behaviour I have
#'     (attempted to) replicate here. It is not clear to me that this
#'     is unbiased, since a clade undergoing rapid diversification
#'     will have many edges but these will be short (and so cannot
#'     have an edge inserted into them following the method below). My
#'     understanding is this is a known problem, and I simply cannot
#'     think of a better way of doing this that doesn't incorporate
#'     what I consider to be worse pathology. Thus this method, even
#'     if it works (which I can't guarantee), it should tend to break
#'     long branches.
#' @param max.iter Sometimes the random draw for the new branch length
#'     to be added will be too large to allow it to be added to the
#'     tree. In such cases, \code{congeneric.imput} will randomly draw
#'     another branch length, and it will repeat this process
#'     \code{max.iter} times. See 'notes' for more on this.
#' @references Steel, M., & Mooers, A. (2010). The expected length of
#'     pendant and interior edges of a Yule tree. Applied Mathematics
#'     Letters, 23(11), 1315-1319.
#' @references Kuhn, T. S., Mooers, A. O., & Thomas, G. H. (2011). A
#'     simple polytomy resolver for dated phylogenies. Methods in
#'     Ecology and Evolution, 2(5), 427-436.
#' @importFrom ape node.depth.edgelength bind.tree is.ultrametric getMRCA 
#' @importFrom phytools getDescendants
#' @rdname phy.build
#' @name phy.build
#' @export
#' @examples
#' tree <- read.tree(text="((a_a:1,b_b:1):1, c_c:2):1;")
#' tree <- congeneric.impute(tree, c("a_nother", "a_gain", "b_sharp"))
congeneric.impute <- function(tree, species, split="_", max.iter=1000, ...){
    tree$missing.species <- character(0)
    if(!is.ultrametric(tree))
        stop("Cannot safely bind into a non-ultrametric tree")
    div.rate <- (length(tree$tip.label)-2) / sum(tree$edge.length)
    
    for(i in seq_along(species)){
        heights <- node.depth.edgelength(tree)
        heights <- max(heights) - heights
        if(!is.null(species[i]) & !species[i] %in% tree$tip.label){
            genus <- strsplit(species[i], split, fixed=TRUE)[[1]][1]
            matches <- unique(grep(genus, tree$tip.label))
            if(length(matches) > 0){
                if(length(matches) > 1){
                node <- getMRCA(tree, matches)
                edges <- getDescendants(tree, node)
                orig.edge <- sample(edges, 1)
                } else {
                    # Only one member of clade in tree - can't get MRCA of a species
                    edges <- orig.edge <- matches
                }
                failed <- TRUE
                for(j in seq_len(max.iter)){
                    rnd.br <- rexp(1, div.rate)
                    old.br <- heights[orig.edge]
                    if(rnd.br > tree$edge.length[which(tree$edge[,2]==orig.edge)])
                        next
                    to.bind <- list(edge=matrix(2:1, 1, 2), tip.label=species[i], edge.length=old.br+rnd.br, Nnode=1)
                    class(to.bind) <- "phylo"
                    tree <- bind.tree(tree, to.bind, where=orig.edge, position=rnd.br)
                    failed <- FALSE
                    break
                }
                if(failed)
                    tree$missing.species <- append(tree$missing.species, species[i])
            } else {
                tree$missing.species <- append(tree$missing.species, species[i])
            }
        }
    }
    return(tree)
}

