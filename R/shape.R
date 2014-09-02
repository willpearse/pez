#' Calculate shape phylogenetic biodiversity metrics across communities
#'
#' \code{shape} calculates phylogenetic biodiversity metrics
#'
#' @param data a \code{comparative.comm} object
#' @param metric specify particular metrics to calculate, default is \code{all}
#' @param which.eigen The eigen vector to calculate for the PhyloEigen metric
#' @param ... Additional arguments to passed to metric functions (unlikely you will want this!)
#' @details Calculates various metrics of phylogenetic biodiversity that are categorized as \emph{shape} metrics by Pearse \emph{et al.} (2014)
#' @return a \code{phy.structure} list object of metric values
#' @author M.R. Helmus, Will Pearse
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus M.R. (2014). Metrics and Models of Community Phylogenetics. In: Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Springer Berlin Heidelberg, pp. 451-464.
#' @references \code{PSV,PSR} Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007). Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83.
#' @references \code{PD} Faith D.P. (1992). Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references \code{colless} Colless D.H. (1982). Review of phylogenetics: the theory and practice of phylogenetic systematics. Systematic Zoology, 31, 100-104.
#' @references \code{gamma}
#' @references \code{taxon} Clarke K.R. & Warwick R.M. (1998). A taxonomic distinctness index and its statistical properties. J. Appl. Ecol., 35, 523-531.
#' @references \code{eigen.sum} Diniz-Filho J.A.F., Cianciaruso M.V., Rangel T.F. & Bini L.M. (2011). Eigenvector estimation of phylogenetic and functional diversity. Functional Ecology, 25, 735-744.
#' @references \code{cadotte.pd} (i.e., \emph{Eed, Hed}) Cadotte M.W., Davies T.J., Regetz J., Kembel S.W., Cleland E. & Oakley T.H. (2010). Phylogenetic diversity metrics for ecological communities: integrating species richness, abundance and evolutionary history. Ecology Letters, 13, 96-105.
#' @examples \dontrun{
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' (output<-shape(data))
#' str(output)
#' shape(data, "colless")
#' shape(data, "eigen.sum", 2)
#' }
#' @importFrom picante psd mpd pd
#' @importFrom vegan taxondive
#' @importFrom PVR PVRdecomp
#' @importFrom apTreeshape as.treeshape as.treeshape.phylo colless tipsubtree
#' @importFrom ape gammaStat cophenetic.phylo drop.tip
#' @export
shape <- function(data,metric=c("all", "psv", "psr", "mpd", "pd","colless", "gamma", "taxon", "eigen.sum","cadotte.pd", "mfpd"),which.eigen=1, phylogeneticWeight = 0.5, removeErrors = TRUE, ...)
{
# vecnums chooses the eigenvector to calculate sumvar in Diniz-Filho
# J.A.F., Cianciaruso M.V., Rangel T.F. & Bini
# L.M. (2011). Eigenvector estimation of phylogenetic and functional
# diversity. Functional Ecology, 25, 735-744.

#-------------------- Unused (old?) parameters --------------------
# @param pa If TRUE (default), all metrics are calculated across
# presence-absence matrices, and species abundances are ignored
# @param permute Number of permutations for metric (currently only for PCD)
    
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  
  #Setup
  nspp <- ncol(data$comm)
  nsite <- nrow(data$comm)
  SR <- rowSums(data$comm > 0)
  data$comm[data$comm > 0] <- 1
  coefs <- data.frame(row.names=rownames(data$comm))
  if(is.null(data$vcv))
    data$vcv <- cophenetic(data$phy)
  output <- list(psv=NULL, psr=NULL, mpd=NULL, pd=NULL, pd.ivs=NULL, colless=NULL, gamma=NULL, taxon=NULL, eigen.sum=NULL, cadotte.pd=NULL)
  
  #Calculate measures
  if(metric == "psv" | metric == "all")
    output$psv <- coefs$psv <- try(psd(data$comm, data$phy, ...)[,1], silent = TRUE)
  
  if(metric == "psr" | metric == "all")
    output$psr <- coefs$psr <- try(psd(data$comm, data$phy, ...)[,4], silent = TRUE)
  
  if(metric == "mpd" | metric == "all")
    output$mpd <- coefs$mpd <- try(mpd(data$comm, data$vcv, abundance.weighted=FALSE, ...), silent = TRUE)
  
  if(metric == "pd" | metric == "all")
    try({
      output$pd <- coefs$pd <- pd(data$comm, data$phy, ...)[,1]
      output$pd.ivs <- coefs$pd.ivs <- unname(resid(lm(coefs$pd ~ rowSums(data$comm))))
    }, silent = TRUE)
  
  if(metric == "colless" | metric == "all")
      try({
          tree.shape <- as.treeshape(data$phy)
          output$colless <- coefs$colless <- apply(data$comm, 1, .colless, tree.shape, tree.shape$names)
      }, silent = TRUE)
  
  
  if(metric == "gamma" | metric == "all")
      try({
          tree.shape <- as.treeshape(data$phy)
          nams <- tree.shape$names
          output$gamma <- coefs$gamma <- apply(data$comm, 1, .gamma, data$phy, nams)
      }, silent = TRUE)
  
  #Note - I've cut out some here because I think the simplification
  #can happen in the summary output
  if(metric == "taxon" | metric == "all")
      try({
        output$taxon <- taxondive(data$comm, data$vcv, ...)
        t <- data.frame(output$taxon$D, output$taxon$Dstar, output$taxon$Lambda, output$taxon$Dplus, output$taxon$SDplus)
        names(t) <- c("Delta", "DeltaStar", "LambdaPlus", "DeltaPlus", "S.DeltaPlus")
        coefs <- cbind(coefs, t)
      }, silent = TRUE)
  
  if(metric == "eigen.sum" | metric == "all")
      try({
          evc <- PVRdecomp(data$phy, ...)@Eigen$vectors
          output$eigen.sum <- coefs$eigen.sum <- apply(data$comm, 1, .eigen.sum, evc, which.eigen)
      }, silent = TRUE)
  
  if(metric == "cadotte.pd" | metric == "all")
    try({
      output$cadotte.pd <- data.frame(Eed=.eed(data), Hed=.hed(data))
      coefs$EED <- output$cadotte.pd$Eed
      coefs$HED <- output$cadotte.pd$Hed
    }, silent = TRUE)

  if(removeErrors) output <- lapply(output, .removeErrors)
  
  #Prepare output
  output$type <- "shape"
  output$coefs <- coefs
  class(output) <- "phy.structure"
  return(output)
}

#Internal Colless function
#' @importFrom apTreeshape colless tipsubtree
.colless <- function(pa.vec,tree,nams)
{
  if(sum(pa.vec)<3)
  {
    return(NA)
  } else {
    return(colless(tipsubtree(tree, nams[pa.vec!=0])))
  }
}

.gamma<-function(pa.vec,tree,nams){
    if(sum(pa.vec)<3){
        return(NA)
    } else {
        if(length(setdiff(tree$tip.label, nams)) != 0){
            tree <- drop.tip(setdiff(tree$tip.label, ))
        }
        return(gammaStat(drop.tip(tree,nams[pa.vec==0])))
    }
}

#' @importFrom picante pd
#' @importFrom picante evol.distinct
.hed <- function(data){
    #Argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")

    #Setup
    ed <- evol.distinct(data$phy, "fair.proportion")$w
    pd <- pd(data$comm, data$phy)$PD

    #Internal assemblage calc.
    ..hed <- function(ed, pd.comm){
        hed <- ed / pd.comm
        return(-sum(hed * log(hed)))
    }

    #Calculate, clean, and return
    output <- numeric(nrow(data$comm))
    names(output) <- rownames(data$comm)
    for(i in seq(nrow(data$comm)))
        output[i] <- ..hed(ed, pd[i])
    return(output)
}

.eed <- function(data, na.rm=TRUE) {
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    output <- .hed(data) / log(apply(data$comm, 1, function(x) sum(x != 0)))
    names(output) <- rownames(data)
    return(output)
}

.eigen.sum <- function(x, evc, vecnums) {
    if(sum(x>0)) {
        return(sum(apply(as.matrix(evc[x>0,vecnums]),2,var)))
    } else {
        return(NA)
    }
}
