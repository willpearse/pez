#' Plots a barplot of a trait in a comparative community object
#' 
#' \code{cc.barplot} plots a barplot of a trait in a comparative community object
#' 
#' @param data a comparative community ecology object
#' @param trait name of trait you wish to plot (only one!)
#' @details This may not necessarily be a good thing to do!
#' @return Named numeric vector, where each element is a trait or community.
#' @author Matt Helmus, Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
#' cc,barplotc(data, c("traitA", "traitB"))
#' }
#' @importFrom ade4 newick2phylog
#' @importFrom ape plot.phylo tiplabels
#' @export
cc.barplot <- function(data, traits, scaling = FALSE, ranging = TRUE, yranging = NULL, tipsontop=TRUE,ceti = 1, csub = 0, f.phylog = 1/(1 + ncol(values)), clabel.leaves = 0, grid.col = gray(0.85), connect.col = "black",vert.line.col = grey(0.5), vert.line.lty = 1, bar.col = heat.colors(dim(values)[1], alpha = .8), spacing = 2.5,border.col = "black",draw.connecting.lines = FALSE, draw.grid = FALSE, draw.vert.line = TRUE, phylog.lwd=1, ...)
  
{
  #Assertions and argument handling
  if (!inherits(data, "comparative.comm"))
    stop("ERROR:", deparse(substitute(data)), "is not of class 'comparative.comm'")
  phylog <- newick2phylog(write.tree(data$phy))
  if(is.null(data$traits))
    stop("ERROR:", deparse(substitute(data)), "does not contain traits to plot!")
  if(any(!traits %in% names(data$traits))){
    missing <- paste(traits[!traits %in% names(data$traits)], collapse=" ")
    stop("ERROR:", missing, "not (a) trait(s) in", deparse(substitute(data)))
  }
  values <- data$traits[,names(data$traits) %in% traits]
  if (!is.numeric(as.matrix(values)))
    stop("ERROR: specified trait(s) not numeric")
  
  names.var <- colnames(values)
  if (scaling == TRUE) {
    values <- scalewt(values)
    values <- as.data.frame(values)
    names(values) <- names.var
  }
  values<-as.data.frame(values)
  lwd.old <- par("lwd")
  par(lwd=phylog.lwd)
  w <- plot(x = phylog, y = NULL, f.phylog = f.phylog, clabel.leaves = clabel.leaves)
  par(lwd=lwd.old)
  mar.old <- par("mar")
  on.exit(par(mar = mar.old))
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  par(usr = c(0, 1, -0.05, 1))
  x1 <- w$xbase
  space <- (1 - w$xbase - (w$xbase - max(w$xy$x))/2 * ncol(values))/ncol(values)
  x2 <- x1 + space
  fun1 <- function(x) {
    x1 + (x2 - x1) * (x - x1.use)/(x2.use - x1.use)
  }
  ret <- cbind.data.frame(values, w$xy[, "y"])
  for (i in seq(ncol(values))) {
  if (ranging == TRUE) {
    if (is.null(yranging)) {val.ref <- pretty(range(values), 4)} else {val.ref <- pretty(yranging, 4)}
  } else {val.ref <- pretty(values[, i], 4)}
  x1.use <- min(val.ref)
  x2.use <- max(val.ref)
  xleg <- fun1(val.ref)
  miny <- 0
  maxy <- max(w$xy$y)
  nleg <- length(xleg)
  if(draw.connecting.lines){segments(w$xy$x, w$xy$y, rep(max(w$xy$x), nrow(values)), w$xy$y,col = connect.col)}
  if(draw.grid)
  {
    segments(xleg, rep(miny, nleg), xleg, rep(maxy, nleg), col = grid.col)
    segments(rep(xleg[1], nrow(values)), w$xy$y, rep(max(xleg), nrow(values)), w$xy$y, col = grid.col)
  }
  if (ceti > 0) {
    if (trunc(i/2) < (i/2))
    {
      y.a<-((miny - 0.05) * 2/3)
      text(xleg, rep(y.a,nleg), as.character(val.ref), cex = par("cex") * ceti)
    } else {
      y.a<-((miny - 0.05) * 1/3)
      text(xleg, rep(y.a,nleg) , as.character(val.ref), cex = par("cex") * ceti)
    }
    lines(c(xleg[1],xleg[length(xleg)]),c(miny,miny))
  }
  
  #draw the bars
  if(scaling){origin<-0} else {origin<-min(val.ref)}
  
  if(draw.vert.line) {segments(fun1(origin), miny, fun1(origin), maxy, lty = vert.line.lty, col = vert.line.col)}
  spacin<-(w$xy$y[1]-w$xy$y[2])/spacing
  for(jj in seq(nrow(values)))
  {
    
    polygon(c(fun1(values[jj, i]),fun1(values[jj, i]),fun1(origin),fun1(origin)),
            c(w$xy$y[jj]-spacin, w$xy$y[jj]+spacin, w$xy$y[jj]+spacin, w$xy$y[jj]-spacin), border = border.col, col = bar.col[jj])
    if(tipsontop){text(w$xy$x[jj]+spacin/2,w$xy$y[jj],rownames(w$xy)[jj],pos=4)}
    
  }
  
  if (csub > 0){text(xleg[3], 1 - (1 - max(w$xy$y))/3, names(values)[i],cex = par("cex") * csub)}
  ret[, i] <- fun1(values[, i])
  x1 <- x1 + space + (w$xbase - max(w$xy$x))/2
  x2 <- x2 + space + (w$xbase - max(w$xy$x))/2
  }
  return(invisible(ret))
}

#' Plots dot-plots of community presence/absence or abundance
#' 
#' \code{plot.comparative.comm} plots a barplot of a trait in a comparative community object
#' 
#' @param data a comparative community ecology object
#' @param sites names of sites to plot (default: all); see examples
#' @param abundance make size proportional to species abundance (default: FALSE)
#' @param dot.cex function to determine point size; see examples, this isn't as terrible-sounding as it seems.
#' @param site.col colours to use when plotting sites; if not same length as number of sites, only the first element is used (no recycling)
#' @param fraction fraction of plot window to be taken up with phylogeny; e.g., 3 (default) means phylogeny is 1/3 of plot
#' @param x.increment specify exact spacing of points along plot; see examples
#' @details Getting the right spacing of dots on the phylogeny may take some playing around with the fraction and x.increment arguments, see below. It may seem a little strange to set point size using a function; however, this gives you much more flexibility when dealing with abundance data, and allows you to transform the data - see the examples
#' @return List containing plot.phylo information, as well as the used x.adj values
#' @author Will Pearse, Matt Helmus
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
#' cc.dotplot(data)
#' cc.dotplot(data, sites=c("clump1", "clump2a"), fraction=1.5)
#' settings <- cc.dotplot(data, sites=c("clump1", "clump2a"), site.col=rainbow(2), fraction=1.5, cex=3)
#' cc.dotplot(data, sites=c("clump1", "clump2a"), site.col=rainbow(2), fraction=1.2, x.increment=settings$x.increment/2)
#' cc.dotplot(data, sites=c("clump1", "clump2a"), site.col=rainbow(2), fraction=1.2, x.increment=settings$x.increment/2, abundance=TRUE, dot.cex=function(x) 3*x)
#' abund.sqrt <- function(x) ifelse(x>0, sqrt(x), 0)
#' cc.dotplot(data, sites=c("clump1", "clump2a"), site.col=rainbow(2), fraction=1.2, x.increment=settings$x.increment/2, abundance=TRUE, dot.cex=abund.sqrt)
#' }
#' @export
plot.comparative.comm <- function(x, sites=NULL, abundance=FALSE, dot.cex=NULL, site.col="black", pch=20, fraction=3, x.increment=NULL, ...){
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
    tiplabels(tip=match(colnames(sites), data$phy$tip.label), pch=pch, col=site.col[i], adj=x.adj[i], cex=sites[i,], ...)
  
  #Invisibly return
  output <- display
  output$x.increment <- x.adj
  invisible(output)
}
