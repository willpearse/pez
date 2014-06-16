##
## Shannon's index
##

setGeneric("shannon",
  function(x) {
    standardGeneric("shannon")
})

setMethod("shannon", signature(x="phylo4d"),
  function(x) {
    phyc <- phylo4com(x)
    shannon(phyc)
})

setMethod("shannon", signature(x="phylo4com"),
  function(x) {
    N.relative <- prop.table(t(abundance(x)), 1)
    out <- apply(N.relative, 1, function(n) -sum(n * log(n)))
    return(out) 
})

