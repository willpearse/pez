#' Calculate D(c) across communities
#' 
#' \code{D.c} runs the extension of D (Dc) across a comparative community object
#' 
#' @param data a comparative community ecology object
#' @param permute the number of null permutations to perform
#' @param traits if TRUE (default) calculate phylogenetic signal of species' traits; otherwise calculate it of community abundance
#' @details This is heavily based on David Orme's implementation of 'phylo.d', which is in turn heavily based on
#' Susanne Fritz's implementation of 'phylo.d'. There is no 'D' function, because Dc calculated across binary
#' data should be equivalent to D. Please don't use this without contacting Will Pearse first - this isn't published, and I'd like to be given some credit for it!
#' Dc=0 is the Brownian expectation, Dc=1 is the random expectation.
#' @return Matrix with columns for each trait or community, and rows for Dc value, and P(Dc=0) and P(Dc=1).
#' @author Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample)
#' D.c(data, permute=1000)
#' D.c(data, traits=FALSE)
#' }
#' @importFrom caper contrCalc VCV.array
#' @importFrom mvtnorm rmvnorm
#' @export
D.c <- function(data, permute=1000, traits=TRUE) {
  #Checking
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  if (!is.numeric(permute)) (stop("'", permute, "' is not numeric.")) 
  # check tree branch lengths
  el    <- data$phy$edge.length
  elTip <- data$phy$edge[,2] <= length(data$phy$tip.label)
  
  if(any(el[elTip] == 0)) 
    stop('Phylogeny contains pairs of tips on zero branch lengths, cannot currently simulate')
  if(any(el[! elTip] == 0)) 
    stop('Phylogeny contains zero length internal branches. Use di2multi.')


  #Internal Dc calculation
  .d.c <- function(ds, vcv, permute, phy){
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
  
  if(traits){
    if(is.null(data$traits)) stop("Need trait data to calculate D.c of traits!")
    vals <- matrix(ncol=3, nrow=ncol(data$traits))
    rownames(vals) <- colnames(data$traits)
    colnames(vals) <- c("D.c", "P(D.c=1)", "P(D.c=0)")
    for(i in seq(ncol(data$traits)))
      vals[i,] <- .d.c(data$traits[,i], vcv, permute, data$phy)
  } else {
    vals <- matrix(ncol=3, nrow=nrow(data$comm))
    rownames(vals) <- rownames(data$comm)
    colnames(vals) <- c("D.c", "P(D.c=1)", "P(D.c=0)")
    for(i in seq(nrow(data$comm)))
      vals[i,] <- .d.c(data$comm[i,], vcv, permute, data$phy)
  }
  return(vals)
}
