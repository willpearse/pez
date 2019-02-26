#' \code{fibre.plot} (fibrously) plots a phylogeny
#' 
#' @param tree a phylogeny (of class phylo) you wish to plot
#' @param gif name of GIF you would like to create. This should *not*
#'     including a folder name (this is due to the use of
#'     \code{\link[animation]{saveGIF}}); "plot.gif" would be fine,
#'     but "work/plot.gif" would not
#' @param focal species numbers or clade numbers to plot differently
#'     (see examples). Note that specifying a clade will highlight the
#'     clade *before* it arises; this is by design. If not specified
#'     (the default) there will be no focal species; this is fine.
#' @param frames number of frames for animation; this will also
#'     determine the time internals for the plot
#' @param colours a function that will return a colour ramp for use in
#'     plotting of species on the fiber plot itself as well as the
#'     standard phylogeny to the right (e.g., \code{\link{rainbow}})
#' @param f.colours as \code{colours} but for the \code{focal} species
#' @param pca PCA (of class \code{\link[stats]{prcomp}}) of
#'     phylogenetic dissimilarity matrix; NULL calculates one, I
#'     recommend you use the output from a previous run to speed
#'     things up
#' @param clade.mat clade matrix (from
#'     \code{\link[caper]{clade.matrix}}$clade.matrix) of phylogeny; NULL
#'     calculates one, I recommend you use the output from a previous
#'     run to speed things up
#' @param delay the delay between each slice's frame in the output
#'     GIF; default 0.2 seconds
#' @param side.tree whether to plot a standard phylogeny to the right
#'     of the plot to aid with interpretation (default: TRUE). You
#'     almost certainly want this option
#' @param width width of animation
#' @param height height of animation
#' @details Probably best to just plot it out and see what happens, to
#'     be honest.
#' @return The data that were plotted last, the PCA and clade.matrix
#'     to speed later plots, and the colours used.
#' @note I would be grateful if you could cite the article this code
#'     was released in when using this code. I maintain this code in
#'     the package "willeerd" on GitHub. I give an example of how to
#'     install this code from there below. Updates will be released
#'     through that, and I welcome code improvements!
#' @author Will Pearse
#' @examples \dontrun{
#' fibre.plot(rlineage(0.1,0), "Yule_fibre.gif")
#' }
#' @importFrom caper clade.matrix
#' @importFrom ape branching.times
#' @importFrom animation saveGIF
#' @importFrom stats prcomp rexp
#' @importFrom graphics layout rect
#' @importFrom grDevices colorRampPalette
#' @export
fibre.plot <- function(tree, gif, focal, frames=60, colours=colorRampPalette(c("blue","black","red")), f.colours=colorRampPalette(c("darkgreen","lightgreen")), pca=NULL, clade.mat=NULL, delay=0.2, side.tree=TRUE, width=600, height=600){
    # Internal function to wrap fibre positions around the plotting matrix
    .snake <- function(x){
        rev.x <- x[,ncol(x):1]
        dim <- dim(x)
        x <- unname(matrix(x, byrow=TRUE))
        rev.x <- unname(matrix(rev.x, byrow=TRUE))
        return(matrix(ifelse(rep(c(FALSE,TRUE),each=dim[2],length.out=length(x)), x, rev.x), dim[1], byrow=TRUE))
    }
    
    ## Assertions and argument checking
    if(!inherits(tree, "phylo"))
        stop("'", deparse(substitute(tree)), "' must be of class 'phylo'")
    if(!is.numeric(frames))
        stop("'", deparse(substitute(frames)), "' must be a numeric!")
    if(!is.function(colours) | !is.function(f.colours))
        stop("'", deparse(substitute(colours)), "' and '", deparse(substitute(f.colours)), "' must be functions (that return colours)")
    if(!is.null(pca) & !inherits(pca, "prcomp"))
        stop("'", deparse(substitute(pca)), "' must be of class 'prcomp'")
    if(!is.character(gif))
        stop("'", deparse(substitute(gif)), "' needs to be a filename!")
    
    ## Setup
    timing <- branching.times(tree)
    colours <- colours(frames)
    f.colours <- f.colours(frames)
    slices <- seq(from=max(timing)+(.Machine$double.eps*5), to=min(c(0,min(timing)))-(.Machine$double.eps*5), length.out=frames)
    if(is.null(pca))
        pca <- prcomp(cophenetic(tree), scale=TRUE, center=TRUE)
    if(is.null(clade.mat))
        clade.mat <- clade.matrix(tree)$clade.matrix

    spp.val <- pca$x[,1]
    dimension <- floor(sqrt(length(spp.val))) + 1
    
    # Setup loop start
    curr.colours <- rep(1, length(tree$tip.label))
    x <- 2
    curr.colours <- rep(NA, dimension^2)
    curr.colours[seq_along(spp.val)] <- 1

    # Setup colours (aside from focals)
    t <- c(timing, setNames(rep(0,length(tree$tip.label)),seq_along(tree$tip.label)))
    t <- t[match(tree$edge[,1], names(t))]
    branch.colours <- rev(colours)[as.numeric(cut(t, slices))]
    
    # Setup focal species
    if(!missing(focal)){
        mask <- rep(FALSE, length(curr.colours))
        if(is.character(focal)){
            mask[which(tree$tip.label %in% focal)] <- TRUE
        } else {
            mask[which(colSums(clade.mat[focal,,drop=FALSE]) > 0)] <- TRUE
        }
        focal.edges <- unique(c(
            which(rowSums(clade.mat[focal,,drop=FALSE])>0),
            which(rowSums(clade.mat[,Filter(function(x) x<=length(tree$tip.label), focal),drop=FALSE])>0)
        ))
        f.branches <- Filter(Negate(is.na), match(tree$edge[,2], focal.edges))
        branch.colours[f.branches] <- rev(f.colours)[as.numeric(cut(t[f.branches], slices))]
    } else {
        mask <- rep(FALSE, length(curr.colours))
    }
    
    ## Loop over all slices and print
    saveGIF({
        if(side.tree)
            layout(matrix(1:2, 1, 2), c(.7,.3), rep(1,2))
        image(t(.snake(matrix(ifelse(mask,NA,curr.colours), nrow=dimension,byrow=TRUE))), col=colours, main=round(slices[1]), zlim=c(1, length(colours)), bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
        image(t(.snake(matrix(ifelse(mask, curr.colours, NA), nrow=dimension,byrow=TRUE))), col=f.colours, main=round(slices[1]), zlim=c(1, length(colours)), bty="n", xaxt="n", yaxt="n", ylab="", xlab="", add=TRUE)
        for(i in seq_along(slices)[-1]){
            to.update <- names(timing)[timing > slices[i] & timing <= slices[i-1]]
            for(j in seq_along(to.update))
                curr.colours[which(clade.mat[to.update[j],] == 1)] <- i
            if(side.tree)
                layout(matrix(1:2, 1, 2), c(.7,.3), rep(1,2))
            image(t(.snake(matrix(ifelse(mask,NA,curr.colours), nrow=dimension,byrow=TRUE))), col=colours, main=paste(round(slices[i],2),"mya"), zlim=c(1, length(colours)), bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
            image(t(.snake(matrix(ifelse(mask,curr.colours,NA), nrow=dimension,byrow=TRUE))), col=f.colours, main=round(slices[1]), zlim=c(1, length(colours)), bty="n", xaxt="n", yaxt="n", ylab="", xlab="", add=TRUE)
            if(side.tree){
                t <- plot(tree, direction="down", show.tip.label=FALSE, edge.color=branch.colours)
                rect(t$x.lim[1]-10, 0, t$x.lim[2]+10, slices[i]-min(slices), col="white", border=NA)
                # weirdness in y-definition needed to account for non-ultrametric trees (with negative branching times)
            }
        }
    }, interval = delay, movie.name = gif, ani.width = width, ani.height = height)
    
    ## Invisibly return
    invisible(list(plot=matrix(curr.colours, nrow=dimension), colours=curr.colours, clade.mat=clade.mat, pca=pca))
}
