#TODO:
# Carry through arguments

#' Calculate dispersion metrics across communities
#' 
#' \code{dispersion} calculates dispersion metrics in comparative.comm communities
#' 
#' @param data a comparative community ecology object
#' @param permute the number of null permutations to perform
#' @param metric specify (a) particular metric(s) to calculate (sesmpd, sesmntd, sespd, innd, d), or the default 'all'
#' @details This calculates $SES_{MPD}$, $SES_{MNTD}$, $SES_{PD}$, INND, and D. All these are defined as dispersion metrics in Pearse et al., 
#' Dc=0 is the Brownian expectation, Dc=1 is the random expectation.
#' @return cc.dispersion object (a named list with the output from each metric)
#' @author Matt Helmus, Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample)
#' dispersion(data)
#' dispersion(data, 100, "sesmpd")
#' }
#' @import caper picante
#' @export
dispersion <- function(data, metric=c("all", "sesmpd", "sesmntd", "sespd", "innd", "d"), permute=1000)
{
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
