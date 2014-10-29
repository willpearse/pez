#' Bind missing species into a phylogeny based on taxonomy
#' 
#' @param tree \code{\link{ape::phylo}} phylogeny to have those
#' species inserted into it
#' @param lookup vector of species names to be bound into the tree if
#' missing from it
#' @param split the character that splits genus and species names in
#' your phylogeny. Default is \code{_}, i.e. Quercus_robur.
#' @details Binds missing species into a phylogeny by replacing all
#' members of the clade it belongs to with a polytomy. Assumes the
#' \code{tip.labels} represent Latin binomials, split by the
#' \code{split} argument. This code was originally shipped with
#' phyloGenerator - this is the \code{merge} method in that program.
#' @return \code{\link{ape::phylo}} phylogeny
#' @author Will Pearse
#' @importFrom ape drop.tip
#' @rdname phy.build
#' @name phy.build
#' @references Pearse W.D. & Purvis A. phyloGenerator: an automated phylogeny generation tool for ecologists. Methods in Ecology and Evolution 4(7): 692--698.
#' @export
merge.phylo <- function(tree, species, split="_"){
    if(!inherits(phy, "phylo"))
        stop("'tree' must be a phylogeny")
    before <- sum(lookup$clean %in% tree$tip.label)
    for(i in seq(nrow(lookup))){
        if(!is.null(lookup$clean[i]) & !lookup$clean[i] %in% tree$tip.label){
            genus <- strsplit(lookup$clean[i], split, fixed=TRUE)[[1]][1]
            matches <- unique(grep(genus, tree$tip.label, value=TRUE))
            if(length(matches) > 0){
                tree <- drop.tip(tree, matches[-1])
                tip.length <- .find.unique.branch.length(tree, matches[1])
                polytomy <- .make.polytomy(unique(c(matches, lookup$clean[i])), (tip.length/2))
                tree <- bind.replace(tree, polytomy, matches[1], tip.length)
            }
        }
    }
    cat("\nNumber of species in tree before:", before)
    cat("\nNumber of species in tree now:   ", sum(lookup$clean %in% tree$tip.label), "\n")
    return(tree)
}

#' Bind one phylogeny into another, giving an ultrametric phylogeny
#' 
#' @param backbone backbone phylogeny (\code{\link{ape::phylo}}) into
#' which the donor is to be bound
#' @param donor phylogeny (\code{\link{ape::phylo}}) to bound into the
#' backbone phylogeny
#' @param replacing.tip.label the species in the donor phylogeny
#' that's being replaced by the donor phylogeny
#' @param donor.length how deep the donor phylogeny should be cut into
#' the backbone phylogeny. If NA (default), then the bladj algorithm
#' is followed (or, in plain English, it's put half-way along the
#' branch)
#' @details Binds a phylogeny (donor) into a bigger phylogeny
#' ('backbone'); useful if you're building a phylogeny a la
#' Phylomatic. This code was originally shipped with phyloGenerator -
#' this is the \code{merge} method in that program.
#' @return phylogeny (\code{\link{ape::phylo}})
#' @author Will Pearse
#' @importFrom ape bind.tree
#' @rdname phy.build
#' @references Pearse W.D. & Purvis A. phyloGenerator: an automated phylogeny generation tool for ecologists. Methods in Ecology and Evolution 4(7): 692--698.
#' @export
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
