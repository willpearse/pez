iac <- function(data, na.rm=TRUE) {
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    
    subtrees <- assemblage.phylogenies(data)
    .denom <-  function(tree, template) {
        # Count number of lineages originating at each internal node
        # (i.e. number of splits)
        int.nodes <- nodeId(data$phy, "internal")
        nSplits <- sapply(int.nodes, function(x) length(children(tree,
            x)))
        names(nSplits) <- int.nodes
        # For each tip, take the product of the number of splits across
        # all of its ancestral nodes
        res <- sapply(ancestors(tree, tipLabels(tree)), function(x)
            prod(nSplits[as.character(x)]))
        template[match(names(res), names(template))] <- res
        return(template)
    }

    # now for each subtree...
    tmp <- setNames(rep(NA, length(data$phy$tip.label)), data$phy$tip.label)
    denom <- lapply(subtrees, .denom, tmp)
    denom <- do.call("cbind", denom)
    nnodes <- sapply(subtrees, function(x) x$Nnode)

    # Calculate expected number of individuals under null hypothesis
    # of equal allocation to each lineage at each (node) split 
    N <- rowSums(data$comm)
    expected <- t(colSums(N, na.rm=na.rm) / t(denom))

    # IAC: summed absolute difference between expected and observed
    # abundances, divided by number of nodes
    return(colSums(abs(expected - N), na.rm=na.rm) / nnodes)
}

