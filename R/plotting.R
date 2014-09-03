#' Plots dot-plots of community presence/absence or abundance
#' @param x a comparative community ecology object (\code{comparative.comm})
#' @param sites names of sites to plot (default: all); see examples
#' @param abundance make size proportional to species abundance (default: FALSE)
#' @param dot.cex function to determine point size; see examples, this isn't as terrible-sounding as it seems.
#' @param site.col colours to use when plotting sites; if not same length as number of sites, only the first element is used (no recycling)
#' @param fraction fraction of plot window to be taken up with phylogeny; e.g., 3 (default) means phylogeny is 1/3 of plot
#' @param x.increment specify exact spacing of points along plot; see examples
#' @details Getting the right spacing of dots on the phylogeny may take some playing around with the fraction and x.increment arguments, see below. It may seem a little strange to set point size using a function; however, this gives you much more flexibility when dealing with abundance data, and allows you to transform the data - see the examples
#' @param ... additional arguments passed to plotting functions
#' @return List containing plot.phylo information, as well as the used x.adj values
#' @author Will Pearse, Matt Helmus
#' @method plot comparative.comm
#' @export
plot.comparative.comm <- function(x, sites=NULL, abundance=FALSE, dot.cex=NULL, site.col="black", fraction=3, x.increment=NULL, ...){
  #Assertions and argument checking
  data <- x
  if (!inherits(data, "comparative.comm"))
    stop("ERROR:", deparse(substitute(data)), "is not of class 'comparative.comm'")
  if(is.null(sites)){
    sites <- rownames(data$comm)
  } else {
    if(any(!sites %in% rownames(data$comm))){
      missing <- paste(sites[!sites %in% rownames(data$comm)], collapse=" ")
      stop("ERROR:", missing, "not site(s) in", deparse(substitute(data)))
    }
  }
  if(!is.null(dot.cex) & !is.function(dot.cex))
    stop("ERROR:", deparse(substitute(dot.cex)), "is not a function!")
  
  #Process data
  sites <- data$comm[sites, ]
  if(!abundance)
    sites <- sites > 0
  if(is.function(dot.cex))
    sites <- dot.cex(sites)
  if(length(site.col) != nrow(sites))
    site.col <- rep(site.col[1], nrow(sites))
  
  #Plot phylogeny and dots
  display <- plot(data$phy, show.tip.label=FALSE, plot=FALSE, no.margin=TRUE, ...)
  plot(data$phy, show.tip.label=FALSE, x.lim=c(0, display$x.lim[2] * fraction), no.margin=TRUE, ...)
  if(!is.null(x.increment)){
    if(length(x.increment) != nrow(sites))
      stop("ERROR:", deparse(substitute(x.increment)), "'s length does not match the number of sites to be plotted") else x.adj <- x.increment
  } else {
    x.adj <- (display$x.lim[2] * fraction) / (nrow(sites) + 2)
    x.adj <- seq(from=x.adj/2, by=x.adj, length.out=nrow(sites))
  }
  for(i in seq(nrow(sites)))
    tiplabels(tip=match(colnames(sites), data$phy$tip.label), col=site.col[i], adj=x.adj[i], cex=sites[i,], ...)
  
  #Invisibly return
  output <- display
  output$x.increment <- x.adj
  invisible(output)
}
