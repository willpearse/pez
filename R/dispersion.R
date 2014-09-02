#TODO:
# Carry through arguments

#' Calculate dispersion phylogenetic biodiversity metrics across communities
#' 
#' \code{dispersion} calculates phylogenetic biodiversity metrics
#' 
#' @param data a \code{comparative.comm} object
#' @param permute the number of null permutations to perform (default 1000)
#' @param metric specify particular metrics to calculate, default is \code{all}
#' @param null.model One of "taxa.labels", "richness", "frequency",
#' "sample.pool", "phylogeny.pool", "independentswap", or
#' "independentswap". These correspond to the null models available in
#' picante; only D does not use these null models currently
#' @param abundance Whether to use abundance-weighted forms of these
#' metrics (default: FALSE)
#' @param ... additional parameters to be passed to metrics (unlikely
#' you will want to use this!)
#' @details Calculates various metrics of phylogenetic biodiversity
#' that are categorized as \emph{dispersion} metrics by Pearse
#' \emph{et al.} (2014). All these are defined as dispersion metrics
#' in Pearse et al., Dc=0 is the Brownian expectation, Dc=1 is the
#' random expectation.
#' @return a \code{phy.structure} list object of metric values
#' @author M.R. Helmus, Will Pearse
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus
#' M.R. (2014). Metrics and Models of Community Phylogenetics. In:
#' Modern Phylogenetic Comparative Methods and Their Application in
#' Evolutionary Biology. Springer Berlin Heidelberg, pp. 451-464.
#' @references \code{sesmpd,sesmntd} Webb C.O. (2000). Exploring the
#' phylogenetic structure of ecological communities: An example for
#' rain forest trees. American Naturalist, 156, 145-155.
#' @references \code{sespd} Webb C.O., Ackerly D.D. & Kembel
#' S.W. (2008). Phylocom: software for the analysis of phylogenetic
#' community structure and trait evolution. Bioinformatics
#' Applications Note, 24, 2098-2100.
#' @references \code{innd} Ness J.H., Rollinson E.J. & Whitney
#' K.D. (2011). Phylogenetic distance can predict susceptibility to
#' attack by natural enemies. Oikos, 120, 1327-1334.
#' @references \code{d} Fritz S.A. & Purvis A. (2010). Selectivity in
#' Mammalian Extinction Risk and Threat Types: a New Measure of
#' Phylogenetic Signal Strength in Binary Traits. Conservation
#' Biology, 24, 1042-1051.
#' @examples \dontrun{
#' data(phylocom)
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' dispersion(data)
#' dispersion(data, 100, "sesmpd")
#' }
#' @importFrom caper phylo.d
#' @importFrom picante ses.mntd ses.mpd ses.pd
#' @export
dispersion <- function(data, metric=c("all", "sesmpd", "sesmntd", "sespd", "innd", "d"), permute=1000, null.model=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), abundance=FALSE, ...)
{
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  if(permute < 0) stop("Can't have negative null permutations!")
  null.model <- match.arg(null.model)
  coefs <- data.frame(row.names=rownames(data$comm))
  
  #Setup
  tree.dist <- cophenetic(data$phy)
  output <- list(sesmpd=NULL, sesmntd=NULL, sespd=NULL, innd=NULL, d=NULL)
  
  #Caculate measures
  if(metric == "sesmpd" | metric == "all"){
    output$sesmpd <- ses.mpd(data$comm,tree.dist, null.model=null.model, abundance.weighted=FALSE, ...)
    coefs$sesmpd <- output$sesmpd$mpd.obs.z
  }
  
  if(metric == "sesmntd" | metric == "all"){
    output$sesmntd <- ses.mntd(data$comm,tree.dist, null.model=null.model, abundance.weighted=FALSE, ...)
    coefs$sesmntd <- output$sesmntd$mntd.obs.z
  }
  
  if(metric == "sespd" | metric == "all"){
    output$sespd <- ses.pd(data$comm,data$phy, null.model=null.model, abundance.weighted=FALSE, ...)
    coefs$sespd <- output$sespd$pd.obs.z
  }
  
  if(metric == "innd" | metric == "all"){
    output$innd <- ses.mpd(data$comm,1/tree.dist, null.model=null.model, abundance.weighted=FALSE, ...)
    coefs$innd <- output$innd$mpd.obs.z
  }
  
  if(metric == "d" | metric == "all"){
    t <- data
    t$comm[t$comm > 0] <- 1
    output$d <- D.c(t, traits=FALSE)
    coefs$d <- output$d[,1]
  }
  
  #Prepare output
  output$type <- "dispersion"
  output$coefs <- coefs
  class(output) <- "phy.structure"
  return(output)
}
