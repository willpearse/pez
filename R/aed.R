# ED and related methods
# - taken from ecoPD, WIP

ed <- function(data, na.rm=TRUE) {
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")

    #Internal function
    .edi <- function(tree) {
        # set length of root edge to zero
        edgeLength(tree)[edgeId(tree, "root")] <- 0

        all.nodes <- nodeId(tree, type = "all")
        des <- descendants(tree, all.nodes, type="tips")
        nv <- edgeLength(tree, all.nodes) / sapply(des, length)
        names(nv) <- all.nodes

        tip.nodes <- nodeId(tree, "tip")
        anc <- ancestors(tree, tip.nodes, "ALL")

        res <- sapply(anc, function(n) sum(nv[as.character(n)], na.rm=TRUE))
        names(res) <- tipLabels(tree)
        return(res)
    }

    #Do the work
    subtrees <- lapply(assemblage.phylogenies(data), as, "phylo4")
    res <- lapply(subtrees, .edi)[as.character(comms)]
    names(res) <- names(comms)
    return(res)

}

hed <- function(data, na.rm=TRUE) {
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")

    ED <- ed(data)
    PD <- pd(data)
    res <- sapply(seq_along(PD), function(x) {scaledED <- ED[[x]] / PD[[x]]; -sum(scaledED * log(scaledED))})
    names(res) <- names(PD)
    return(res)
}

eed <- function(data, na.rm=TRUE) {
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    subtrees <- lapply(assemblage.phylogenies(data), as, "phylo4")
    output <- hed(data) / log(sapply(subtrees, function(x) length(x$tip.label)))
    return(output)
}


# TODO: This function includes its own code for not counting root edge
# length. Maybe this should maybe be done at a higher level?
aed <- function(data, na.rm=TRUE) {
    #Argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    
    #Internal functions
    .isD <- function(tree) {
        # get all node IDs, but excluding root node
        nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
        # Create logical matrix indicating which tips (in columns) are
        # descendants of each node (in rows), self-inclusive
        t(sapply(ancestors(tree, tipLabels(tree), "ALL"),
                 function(n) nonroot.nodes %in% n))
    }
    .elen <- function(tree) {
        nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
        edgeLength(tree, nonroot.nodes)
    }
    .aed <- function(x){
        spp <- row.names(isDescendant[[i]])        
        dAbund <- N[spp, i] * isDescendant[[i]]
        # Calculate individual-based AED of each species
        AED <- colSums(edge.length[[i]] * t(prop.table(dAbund, margin=2)))
        AED/N[spp, i]
    }

    subtrees <- lapply(assemblage.phylogenies(data), as, "phylo4")
    isDescendant <- lapply(subtrees, .isD)[as.character(comms)]
    edge.length <- lapply(subtrees, .elen)[as.character(comms)]
    N <- abundance(x)

    res <- lapply(seq_along(N), .aed)
    names(res) <- names(comms)
    return(res)

}

haed <- function(x, na.rm=TRUE) {
    # Recast AED in terms of individuals
    AED <- aed(x)
    PD <- pd(x)
    N <- abundance(x)
    scaled.AED <- lapply(seq_along(N), function(i) {
        spp <- names(AED[[i]])        
        rep(unname(AED[[i]]), N[spp, i]) / PD[[i]]
    })
    res <- sapply(scaled.AED, function(x) -sum(x * log(x)))
    names(res) <- names(AED)
    return(res)
}
