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
#' metrics (default: FALSE). Doesn't apply to D, which is
#' presence/absence only.
#' @param ... additional parameters to be passed to metrics (unlikely
#' you will want to use this!)
#' @details Calculates various metrics of phylogenetic biodiversity
#' that are categorized as \emph{dispersion} metrics by Pearse
#' \emph{et al.} (2014).
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
    output$sesmpd <- ses.mpd(data$comm,tree.dist, null.model=null.model, abundance.weighted=abundance, ...)
    coefs$sesmpd <- output$sesmpd$mpd.obs.z
  }
  
  if(metric == "sesmntd" | metric == "all"){
    output$sesmntd <- ses.mntd(data$comm,tree.dist, null.model=null.model, abundance.weighted=abundance, ...)
    coefs$sesmntd <- output$sesmntd$mntd.obs.z
  }
  
  if(metric == "sespd" | metric == "all"){
    output$sespd <- ses.pd(data$comm,data$phy, null.model=null.model, abundance.weighted=abundance, ...)
    coefs$sespd <- output$sespd$pd.obs.z
  }
  
  if(metric == "innd" | metric == "all"){
    output$innd <- ses.mpd(data$comm,1/tree.dist, null.model=null.model, abundance.weighted=abundance, ...)
    coefs$innd <- output$innd$mpd.obs.z
  }
  #WARNING: Altering commnity matrix in-place
  if(metric == "d" | metric == "all"){
    data$comm[data$comm > 0] <- 1
    output$d <- .d(data, traits=FALSE)
    coefs$d <- output$d[,1]
  }
  
  #Prepare output
  output$type <- "dispersion"
  output$coefs <- coefs
  class(output) <- "phy.structure"
  return(output)
}

#Interal D calculation. *Heavily* based on D from caper.
#' @importFrom caper contrCalc VCV.array
#' @importFrom mvtnorm rmvnorm
.d <- function(data, permute=1000, traits=TRUE) {
  #Checking
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  if (!is.numeric(permute)) (stop("'", permute, "' is not numeric."))
  data$comm[data$comm > 1] <- 1
  # check tree branch lengths
  el    <- data$phy$edge.length
  elTip <- data$phy$edge[,2] <= length(data$phy$tip.label)
  
  if(any(el[elTip] == 0)) 
    stop('Phylogeny contains pairs of tips on zero branch lengths, cannot currently simulate')
  if(any(el[! elTip] == 0)) 
    stop('Phylogeny contains zero length internal branches. Use di2multi.')


  #Internal D calculation
  ..d <- function(ds, vcv, permute, phy){
      dsSort <- sort(ds)
      
      ## Random Association model
      ds.ran <- replicate(permute, sample(ds))
      
      ## Brownian Threshold model random data
      ds.phy <- rmvnorm(permute, sigma=unclass(vcv)) # class of 'VCV.array' throws the method dispatch
      ds.phy <- as.data.frame(t(ds.phy))
      
                                        # turn those into rank values, then pull that rank's observed value
      ds.phy <- apply(ds.phy, 2, rank, ties="random")
      ds.phy <- apply(ds.phy, 2, function(x) as.numeric(dsSort[x]))
      
      ## Get change along edges
      ## insert observed and set dimnames for contrCalc
      ds.ran <- cbind(Obs=ds, ds.ran)
      ds.phy <- cbind(Obs=ds, ds.phy)
      dimnames(ds.ran) <- dimnames(ds.phy) <- list(data$phy$tip.label, c('Obs', paste('V',1:permute, sep='')))
      
      ## now run that through the contrast engine 
      ds.ran.cc <- contrCalc(vals=ds.ran, phy=phy, ref.var='V1', picMethod='phylo.d', crunch.brlen=0)
      ds.phy.cc <- contrCalc(vals=ds.phy, phy=phy, ref.var='V1', picMethod='phylo.d', crunch.brlen=0)
      
      ## get sums of change and distributions
      ransocc <- colSums(ds.ran.cc$contrMat)
      physocc <- colSums(ds.phy.cc$contrMat)
      # double check the observed, but only to six decimal places or you can get floating point errors
      if(round(ransocc[1], digits=6) != round(physocc[1], digits=6)) stop('Problem with character change calculation in phylo.d')
      obssocc <- ransocc[1]
      ransocc <- ransocc[-1]
      physocc <- physocc[-1]
      
      soccratio <- (obssocc - mean(physocc)) / (mean(ransocc) - mean(physocc))
      soccpval1 <- sum(ransocc < obssocc) / permute
      soccpval0 <- sum(physocc > obssocc) / permute
      
      return(c(soccratio, soccpval1, soccpval0))
  }

  ## being careful with the edge order - pre-reorder the phylogeny
  phy <- reorder(data$phy, 'pruningwise')
  
  vcv <- VCV.array(data$phy)
  vals <- matrix(ncol=3, nrow=nrow(data$comm))
  rownames(vals) <- rownames(data$comm)
  colnames(vals) <- c("D", "P(D=1)", "P(D=0)")
  for(i in seq(nrow(data$comm)))
      vals[i,] <- ..d(data$comm[i,], vcv, permute, data$phy)
  
  return(vals)
}
