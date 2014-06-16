simpson <- function(x, method=c("phylogenetic", "traditional")) {
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    method <- match.arg(method)
    N.relative <- prop.table(t(rowsums(data$comm), 1)
    if (method=="phylogenetic") {
        dmat <- pair.dist(data$phy)
    } else {
        dmat <- matrix(1, length(data$phy$tip.label), length(data$phy$tip.label))
        diag(dmat) <- 0
    }   
    out <- apply(N.relative, 1, function(n) sum((n %o% n)*dmat))
    return(out) 
}
