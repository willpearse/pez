#TODO: I'm not sure PD should be standardised by the sum of species,
# but rather the sum of individuals. In general, we need to have a
# better discussion about individual vs. species metrics and scaling
# it's going to be hard for people to install ecoPD, so we need some
# sensible 'require'-foo to get everything in there I don't like
# having to 'require' ape when using phylobase, but either I'm doing
# something really stupid with the dependencies or phylobase doesn't
# have a properly setup NAMESPACE. Grrrr.

#' Calculate dissimilarity metrics across communities
#' 
#' \code{dissimilarity} calculates dissimilarity metrics in comparative.comm communities
#' 
#' @param data a comparative community ecology object
#' @param metric specify (a) particular metric(s) to calculate (unifrac, pcd, phylosor), or the default 'all'
#' @param pa If TRUE (default), all metrics are calculated across presence-absence matrices, and species abundances are ignored
#' @param permute Number of permutations for metric (currently only for PCD)
#' @param removeErrors If \code{FALSE}, metrics that failed to be computed are returned as error messages
#' @details This calculates UniFrac, PCD (and its subcomponents), and the phylogenetic Sorenson's index, all defined as dissimilarity metrics in Pearse et al.
#' @note This function uses XXX's version of the PCD function, not that included in picante
#' @return cc.dissimilarity object (a named list with the output from each metric)
#' @author Matt Helmus, Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample)
#' shape(data)
#' shape(data, "colless")
#' shape(data, "eigen.sum", 2)
#' }
#' @importFrom picante psd mpd pd
#' @importFrom vegan taxondive
#' @importFrom PVR PVRdecomp
#' @importFrom apTreeshape as.treeshape as.treeshape.phylo colless tipsubtree
#' @importFrom ape gammaStat cophenetic.phylo drop.tip
#' @export
shape <- function(data,
                  metric=c("all", "psv", "psr", "mpd", "pd",
                      "colless", "gamma", "taxon", "eigen.sum", "cadotte.pd"),
                  which.eigen=1, removeErrors = TRUE)
{
# vecnums chooses the eigenvector to calculate sumvar in Diniz-Filho
# J.A.F., Cianciaruso M.V., Rangel T.F. & Bini
# L.M. (2011). Eigenvector estimation of phylogenetic and functional
# diversity. Functional Ecology, 25, 735-744.
    
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  
  #Setup
  nspp <- ncol(data$comm)
  nsite <- nrow(data$comm)
  SR <- rowSums(data$comm > 0)
  pa <- data$comm > 0
  if(is.null(data$vcv))
    data$vcv <- cophenetic(data$phy)
  output <- list(psv=NULL, psr=NULL, mpd=NULL, pd=NULL, pd.ivs=NULL, colless=NULL, gamma=NULL, taxon=NULL, eigen.sum=NULL, cadotte.pd=NULL)
  
  #Calculate measures
  if(metric == "psv" | metric == "all")
    output$psv <- try(psd(pa, data$phy)[,1], silent = TRUE)
  
  if(metric == "psr" | metric == "all")
    output$psr <- try(psd(pa, data$phy)[,4], silent = TRUE)
  
  if(metric == "mpd" | metric == "all")
    output$mpd <- try(mpd(pa, data$vcv), silent = TRUE)
  
  if(metric == "pd" | metric == "all"){
      output$pd.ivs <- .trywith(output, {
          pd <- pd(pa, data$phy)[,1]
          resid(lm(pd ~ rowSums(pa)))
      })
  }
  
  if(metric == "colless" | metric == "all"){
      output$colless <- .trywith(data, {
          tree.shape <- as.treeshape(phy)
          nams <- tree.shape$names
          apply(pa, 1, .colless, tree.shape, nams)
      })
  }
  
  if(metric == "gamma" | metric == "all"){
      output$gamma <- .trywith(data, {
          tree.shape <- as.treeshape(phy)
          nams <- tree.shape$names
          apply(pa, 1, .gamma, phy, nams)
      })
  }
  
  #Note - I've cut out some here because I think the simplification
  #can happen in the summary output
  if(metric == "taxon" | metric == "all")
      output$taxon <- try(taxondive(pa, data$vcv), silent = TRUE)
  
  if(metric == "eigen.sum" | metric == "all"){
      output$eigen.sum <- .trywith(data, {
          evc <- PVRdecomp(phy)@Eigen$vectors
          apply(pa, 1, 
                function(x, evc, vecnums) {
                    if(sum(x>0)){
                        return(sum(apply(as.matrix(evc[x>0,vecnums]),2,var)))
                    } else {
                        return(NA)
                    }
                }, evc, which.eigen)
      })
  }
  
  if(metric == "cadotte.pd" | metric == "all")
    output$cadotte.pd <- try(data.frame(Eed=.eed(data), Hed=.hed(data)), silent = TRUE)

  if(removeErrors) output <- lapply(output, .removeErrors)
  
  #Prepare output
  output$type <- "shape"
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
#' @importFrom caper ed.calc
.hed <- function(data, na.rm=TRUE) {
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")

    ED <- lapply(assemblage.phylogenies(data), function(x) setNames(ed.calc(x)$spp$ED, x$tip.label))
    PD <- pd(data$comm, data$phy)$PD
    res <- sapply(seq_along(PD), function(x) {scaledED <- ED[[x]] / PD[x]; -sum(scaledED * log(scaledED))})
    names(res) <- rownames(data$comm)
    return(res)
}

.eed <- function(data, na.rm=TRUE) {
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    output <- .hed(data) / log(apply(data$comm, 1, function(x) sum(x != 0)))
    names(output) <- rownames(data)
    return(output)
}


.trywith <- function(data, expr) try(with(data, expr), silent = TRUE)

.removeErrors <- function(object) {
    if(inherits(object, "try-error")) return(NULL)
    object
}
