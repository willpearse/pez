#' Calculate a lambda, delta, or kappa for a comparative community ecology dataset
#' 
#' \code{classic.phy.signal} calculates Pagel's lambda, delta, or kappa for all community abundances or trait values. Essentially a wrapper around 
#' 
#' @param data a comparative community ecology object on which you want to calculate phylogenetic signal
#' @param traits if TRUE (default) calculate phylogenetic signal of species' traits; otherwise calculate it of community abundance
#' @details This may not necessarily be a good thing to do!
#' You may not have trait or community data where these measures make sense!
#' @return Named numeric vector, where each element is a trait or community.
#' @author Will Pearse, Jeannine Cavender-Bares
#' @references R. P. Freckleton, P. H. Harvey, and M. Pagel. Phylogenetic analysis and comparative data: A test and review of evidence. American Naturalist, 160:712-726, 2002.
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
#' classic.phy.signal(data)
#' classic.phy.signal(data, traits=FALSE, method="kappa")
#' }
#' @importFrom caper comparative.data pgls
#' @export
phy.signal <- function(data, traits=TRUE, method=c("lambda", "delta", "kappa")){
  #Assertions and argument handling
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  method <- match.arg(method)
  
  #Traits
  if(traits){
    if(is.null(data$traits)) stop("Need traits to compute phylogenetic signal of traits!")
    data$traits$this.breaks.pez <- rownames(data$traits)
    ## FIXME:  are we supposed to look for this.breaks.pez in data$traits?
    comparative.data <- comparative.data(data$phy, data$traits, data$traits$this.breaks.pez)
    traits <- numeric(ncol(comparative.data$data))
    names(traits) <- names(comparative.data$data)
    for(i in seq(ncol(comparative.data$data))){
      #I know, this should be a case statement...
      if(method == "lambda"){
        model <- pgls(comparative.data$data[,1] ~ 1, data=comparative.data, lambda="ML")
        traits[i] <- summary(model)$param.CI$lambda$opt
      }
      if(method == "delta"){
        model <- pgls(comparative.data$data[,1] ~ 1, data=comparative.data, delta="ML")
        traits[i] <- summary(model)$param.CI$delta$opt
      }
      if(method == "kappa"){
        model <- pgls(comparative.data$data[,1] ~ 1, data=comparative.data, kappa="ML")
        traits[i] <- summary(model)$param.CI$kappa$opt
      }
    }
    return(traits)
  } else {
    #Community composition
    comm <- data.frame(t(data$comm))
    comm$this.breaks.pez <- rownames(comm)
    comparative.data <- comparative.data(data$phy, comm, data$traits$this.breaks.pez)
    communities <- numeric(ncol(comparative.data$data))
    names(communities) <- names(comparative.data$data)
    for(i in seq(ncol(comparative.data$data))){
      if(method == "lambda"){
        model <- pgls(comparative.data$data[,1] ~ 1, data=comparative.data, lambda="ML")
        communities[i] <- summary(model)$param.CI$lambda$opt
      }
      if(method == "delta"){
        model <- pgls(comparative.data$data[,1] ~ 1, data=comparative.data, delta="ML")
        communities[i] <- summary(model)$param.CI$delta$opt
      }
      if(method == "kappa"){
        model <- pgls(comparative.data$data[,1] ~ 1, data=comparative.data, kappa="ML")
        communities[i] <- summary(model)$param.CI$kappa$opt
      }
    }
    return(communities)
  }
}
