## TODO:
##
## 1. Carry through arguments
##
## 2. Make is so that the user can change the arguements in the
## ses. functions such as different null models. MRH does not know how
## to do this, is it with the "..." argument?

## SCW: the list(...) idiom might help with issue #2.  But I wonder
## why not just add these arguments to the list, and then check to see
## whether or not they are being ignored (so the user can be usefully
## informed about whether they've made an honest mistake).  If you
## only need to pass arguments on to one function, then:
##   dispersion <- function( , ...){
##      some code
##      onlyFunctionThatNeedsDots( , ...)
##      some more code
##   }
## would work.

#' Calculate dispersion phylogenetic biodiversity metrics across communities
#' 
#' \code{dispersion} calculates phylogenetic biodiversity metrics
#' 
#' @param data a \code{comparative.comm} object
#' @param permute the number of null permutations to perform
#' @param metric specify particular metrics to calculate, default is
#' \code{all}
#' @param ... additional arguments for modifying the behaviour of
#' various metrics (see details)
#'
#' @details Calculates various metrics of phylogenetic biodiversity
#' that are categorized as \emph{dispersion} metrics by Pearse
#' \emph{et al.} (2014) This calculates $SES_{MPD}$, $SES_{MNTD}$,
#' $SES_{PD}$, INND, and D. All these are defined as dispersion
#' metrics in Pearse et al., Dc=0 is the Brownian expectation, Dc=1 is
#' the random expectation.
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
#
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample)
#' dispersion(data)
#' dispersion(data, 100, "sesmpd")
#' }
#' @importFrom caper phylo.d
#' @importFrom picante ses.mntd ses.mpd
#' @export
dispersion <- function(data, metric = c("all", "sesmpd", "sesmntd", "sespd", "innd", "d"),
                       permute=1000, ...)
{

    ## MRH, is this what you mean?
    l... <- list(...)
    Dc <- l...$Dc
    hasDc <- !is.null(reGenerators) # in case you need to know whether its available
    
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  if(permute < 0) stop("Can't have negative null permutations!")
  
  #Setup
  tree.dist <- cophenetic(data$phy)
  output <- list(sesmpd=NULL, sesmntd=NULL, sespd=NULL, innd=NULL, d=NULL)
  
  #Caculate measures
  if(metric == "sesmpd" | metric == "all")
    output$sesmpd <- ses.mpd(data$comm,tree.dist)
  
  if(metric == "sesmntd" | metric == "all")
    output$sesmntd <- ses.mntd(data$comm,tree.dist)
  
  if(metric == "sespd" | metric == "all")
    output$sespd <- ses.pd(data$comm,data$phy)
  
  if(metric == "innd" | metric == "all")
    output$innd <- ses.mpd(data$comm,1/tree.dist)
  
  if(metric == "d" | metric == "all")
    output$d <- D.c(data, traits=FALSE)
  
  #Prepare output
  output$type <- "dispersion"
  class(output) <- "phy.structure"
  return(output)
}
