#' Calculate (phylogenetic) dispersion: examine assemblages in the
#' context of a source pools
#'
#' As described in Pearse et al. (2014), an evenness metric is one the
#' examines the phylogenetic structure of species present in each
#' assemblage in the context of a source pool of potentially present
#' species. Unlike other metrics, the values of a dispersion metric is
#' *contingent* on the definition of source pool, and (often)
#' randomisations used to conduct that comparison. For completeness,
#' options are provided to calculate these metrics using species
#' traits.
#'
#' Most of these metrics do not involve comparison with some kind of
#' evolutionary-derived expectation for phylogenetic shape. Those that
#' do, however, such as D, make no sense unless applied to a
#' phylogenetic distance matrix - their null expectation *requires*
#' it. Using square-rooted distance matrices, or distance matrices
#' that incorporate trait information, can be an excellent thing to
#' do, but (for the above reasons), \code{pez} won't give you an
#' answer for metrics for which WDP thinks it makes no sense. SESpd
#' can (...up to you whether it should!...) be used with a
#' square-rooted distance matrix, but the results *will always be
#' wrong* if you do not have an ultrametric tree (branch lengths
#' proportional to time) and you will be warned about this. WDP
#' strongly feels you should only be using ultrametric phylogenies in
#' any case, but code to fix this bug is welcome.
#' 
#' @param data \code{\link{comparative.comm}} object
#' @param permute number of null permutations to perform (default
#' 1000)
#' @param metric default (\code{all}) calculates everything;
#' individually call-able metrics are: \code{SESmpd}, \code{SESmntd},
#' \code{SESpd}, \code{innd}, \code{d}.
#' @param null.model one of "taxa.labels", "richness", "frequency",
#' "sample.pool", "phylogeny.pool", "independentswap", or
#' "independentswap". These correspond to the null models available in
#' \code{\link{picante}}; only \code{d} does not use these null models
#' @param abundance Whether to use abundance-weighted forms of these
#' metrics (default: FALSE). D, which is presence/absence only, and so
#' will not be calculated when \code{TRUE}.
#' @param sqrt.phy If TRUE (default is FALSE) your phylogenetic
#' distance matrix will be square-rooted; specifying TRUE will force
#' the square-root transformation on phylogenetic distance matrices
#' (in the spirit of Leitten and Cornwell, 2014). See `details' for
#' details about different metric calculations when a distance matrix
#' is used.
#' @param traitgram If not NULL (default), a number to be passed to
#' \code{funct.phylo.dist} (\code{phyloWeight}; the `a' parameter),
#' causing analysis on a distance matrix reflecting both traits and
#' phylogeny (0 --> only phylogeny, 1 --> only traits; see
#' \code{funct.phylo.dist}). If a vector of numbers is given,
#' \code{shape} iterates across them and returns a \code{data.frame}
#' with coefficients from each iteration. See `details' for details
#' about different metric calculations when a distance matrix is used.
#' @param ext.dist Supply an external species-level distance matrix
#' for use in calculations. See `details' for comments on the use of
#' distance matrices in different metric calculations.
#' @param traitgram.p A value for `p' to be used in conjunction with
#' \code{traitgram} when calling \code{funct.phylo.dist}.
#' @param ... additional parameters to be passed to metrics (unlikely
#' you will want to use this!)
#' @return a \code{phy.structure} list object of metric values
#' @author M.R. Helmus, Will Pearse
#' @seealso \code{\link{shape}} \code{\link{evenness}} \code{\link{dissimilarity}}
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
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' \dontrun{dispersion(data)}
#' dispersion(data, metric = "sesmpd", permute = 100)
#' @importFrom caper phylo.d
#' @importFrom picante ses.mntd ses.mpd ses.pd
#' @importFrom ape is.ultrametric as.phylo
#' @export
dispersion <- function(data, metric=c("all", "sesmpd", "sesmntd", "sespd", "innd", "d"), permute=1000, null.model=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), abundance=FALSE, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, ...)
{
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  if(permute < 0) stop("Can't have negative null permutations!")
  null.model <- match.arg(null.model)
  coefs <- data.frame(row.names=rownames(data$comm))

  if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
      stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
  if(!is.null(traitgram)){
      if(length(traitgram) > 1){
          output <- vector("list", length(traitgram))
          for(i in seq_along(output))
              output[[i]] <- cbind(coef(Recall(data, metric=metric, permute=permute, null.model=null.model, abundance=abundance, sqrt.phy=sqrt.phy, traitgram=traitgram[i], traitgram.p=traitgram.p, ext.dist=ext.dist)), traitgram[i], sites(data))
          output <- do.call(rbind, output)
          names(output)[ncol(output)-1] <- "traitgram"
          names(output)[ncol(output)] <- "site"
          rownames(output) <- NULL
          return(output)
      } else {
          dist <- as.matrix(funct.phylo.dist(data, traitgram, traitgram.p))
          traitgram <- TRUE
      }      
  } else traitgram <- FALSE
  
  if(!is.null(ext.dist)){
      dist <- .check.ext.dist(ext.dist, species(data), ncol(data$comm))
      ext.dist <- TRUE
  } else ext.dist <- FALSE
  
  #Setup
  if(traitgram==FALSE & ext.dist==FALSE)
      dist <- cophenetic(data$phy)
  if(sqrt.phy){
        if(!is.ultrametric(data$phy))
            warning("Phylogeny is not ultrametric; see function details")
      dist <- sqrt(dist)
      data$phy <- as.phylo(hclust(as.dist(dist)))
  }
  dist <- cophenetic(data$phy)
  output <- list(sesmpd=NULL, sesmntd=NULL, sespd=NULL, innd=NULL, d=NULL)
  
  #Caculate measures
  if(metric == "sesmpd" | metric == "all"){
    output$sesmpd <- ses.mpd(data$comm, dist, null.model=null.model, abundance.weighted=abundance, ...)
    coefs$sesmpd <- output$sesmpd$mpd.obs.z
  }
  
  if(metric == "sesmntd" | metric == "all"){
    output$sesmntd <- ses.mntd(data$comm, dist, null.model=null.model, abundance.weighted=abundance, ...)
    coefs$sesmntd <- output$sesmntd$mntd.obs.z
  }
  
  if((metric == "sespd" | metric == "all") & ext.dist==FALSE){
    output$sespd <- ses.pd(data$comm, data$phy, null.model=null.model, ...)
    coefs$sespd <- output$sespd$pd.obs.z
  }
  
  if(metric == "innd" | metric == "all"){
    output$innd <- ses.mpd(data$comm,1/dist, null.model=null.model, abundance.weighted=abundance, ...)
    coefs$innd <- output$innd$mpd.obs.z
  }
  #WARNING: Altering commnity matrix in-place
  if((metric == "d" | metric == "all") & (ext.dist==FALSE & sqrt.phy==FALSE & abundance==FALSE & traitgram==FALSE)){
    data$comm[data$comm > 0] <- 1
    output$d <- .d(data)
    coefs$d <- output$d[,1]
  }
  
  #Prepare output
  output$type <- "dispersion"
  output$coefs <- coefs
  class(output) <- "phy.structure"
  return(output)
}

#Interal D calculation. *Heavily* based on D from caper.
# - placed in here for future compatibility with Dc
#' @importFrom caper contrCalc VCV.array
#' @importFrom mvtnorm rmvnorm
.d <- function(data, permute=1000) {
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
