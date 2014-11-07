#' Calculate (phylogenetic) shape: examine assemblage composition
#'
#' @param data \code{\link{comparative.comm}} object
#' @param metric metrics to calculate. Default (\code{all-quick})
#' calculates everything bar \code{fd.dist}, \code{all} calculates
#' everything. Individually call-able metrics are: \code{psv},
#' \code{psr}, \code{mpd}, \code{mntd}, \code{pd}, \code{colless},
#' \code{gamma}, \code{taxon}, \code{eigen.sum}, \code{eed}, \code{hed}, &
#' \code{dist.fd}. Note that \code{dist.fd} is a trait distance
#' metric, but will be calculated on the phylogenetic distance matrix
#' if not traits are supplied, and the external/square-rooted
#' phylogenetic distance matrix if either is specified.
#' @param which.eigen The eigen vector to calculate for the PhyloEigen
#' metric (\code{eigen.sum})
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
#' @param traitgram.p A value for `p' to be used in conjunction with
#' \code{traitgram} when calling \code{funct.phylo.dist}.
#' @param remove.errors suppress errors about metrics that failed to
#' run (default: TRUE). This will happen for some metrics if you have,
#' for example, only one species in a community.
#' @param ext.dist Supply an external species-level distance matrix
#' for use in calculations. See `details' for comments on the use of
#' distance matrices in different metric calculations.
#' @param ... Additional arguments to passed to metric functions
#' (unlikely you will want this!)
#' @details Calculates various metrics of phylogenetic biodiversity
#' that are categorized as \emph{shape} metrics by Pearse \emph{et
#' al.} (2014).
#' @details As described in Pearse et al. (2014), a shape metric is
#' one the examines the phylogenetic structure of species present in
#' each assemblage, ignoring abundances entirely. For completeness,
#' options are provided to calculate these metrics using species
#' traits.
#' @details Most of these metrics do not involve comparison with some
#' kind of evolutionary-derived expectation for phylogenetic
#' shape. Those that do, however, such as PSV or Colless' index, make
#' no sense unless applied to a phylogenetic distance matrix - their
#' null expectation *requires* it. Using square-rooted distance
#' matrices, or distance matrices that incorporate trait information,
#' can be an excellent thing to do, but (for the above reasons),
#' \code{pez} won't give you an answer for metrics for which WDP
#' thinks it makes no sense. \code{pd}, \code{eed} & \code{hed} can
#' (...up to you whether you should!...) be used with a square-rooted
#' distance matrix, but the results *will always be wrong* if you do
#' not have an ultrametric tree (branch lengths proportional to time)
#' and you will be warned about this. WDP strongly feels you should
#' only be using ultrametric phylogenies in any case, but code to fix
#' this bug is welcome.
#' @note As mentioned above, \code{dist.fd} is calculated using a
#' phylogenetic distance matrix if no trait data are available, or if
#' you specify \code{sqrt.phy}. It is not calculated by default
#' because it generates warning messsages (which WDP is loathe to
#' suppress) which are related to the general tendency for a low rank
#' of phylogenetic distance matrices. Much ink has been written about
#' this, and in par this problem is why the \code{eigen.sum} measure
#' came to be suggested.
#' @return \code{phy.structure} list object of metric values. Use
#' \code{coefs} to extract a summary metric table, or examine each
#' individual metric (which gives more details for each) by calling
#' \code{print} on the output (i.e., type \code{output} in the example
#' below).
#' @author M.R. Helmus, Will Pearse
#' @seealso \code{\link{evenness}} \code{\link{dispersion}} \code{\link{dissimilarity}}
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus M.R. (2014). Metrics and Models of Community Phylogenetics. In: Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Springer Berlin Heidelberg, pp. 451-464.
#' @references \code{PSV,PSR} Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007). Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83.
#' @references \code{PD} Faith D.P. (1992). Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references \code{colless} Colless D.H. (1982). Review of phylogenetics: the theory and practice of phylogenetic systematics. Systematic Zoology, 31, 100-104.
#' @references \code{gamma} Pybus O.G. & Harvey P.H. (2000) Testing macro-evolutionary models using incomplete molecular phylogenies. _Proceedings of the Royal Society of London. Series B. Biological Sciences 267: 2267--2272.
#' @references \code{taxon} Clarke K.R. & Warwick R.M. (1998). A taxonomic distinctness index and its statistical properties. J. Appl. Ecol., 35, 523-531.
#' @references \code{eigen.sum} Diniz-Filho J.A.F., Cianciaruso M.V., Rangel T.F. & Bini L.M. (2011). Eigenvector estimation of phylogenetic and functional diversity. Functional Ecology, 25, 735-744.
#' @references \code{eed,hed} (i.e., \emph{Eed, Hed}) Cadotte M.W., Davies T.J., Regetz J., Kembel S.W., Cleland E. & Oakley T.H. (2010). Phylogenetic diversity metrics for ecological communities: integrating species richness, abundance and evolutionary history. Ecology Letters, 13, 96-105.
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' (output<-shape(data))
#' str(output)
#' shape(data, "colless")
#' shape(data, "eigen.sum", which.eigen=2)
#' @importFrom picante psd mpd pd mntd
#' @importFrom vegan taxondive
#' @importFrom PVR PVRdecomp
#' @importFrom apTreeshape as.treeshape as.treeshape.phylo colless tipsubtree
#' @importFrom ape gammaStat cophenetic.phylo drop.tip is.ultrametric as.phylo
#' @importFrom FD dbFD
#' @export
shape <- function(data,metric=c("all-quick", "all", "psv", "psr", "mpd", "mntd", "pd","colless", "gamma", "taxon", "eigen.sum","eed", "hed", "dist.fd"), sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, which.eigen=1, remove.errors = TRUE,  ...)
{
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
      stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
  if(!is.null(traitgram)){
      if(length(traitgram) > 1){
          output <- vector("list", length(traitgram))
          for(i in seq_along(output))
              output[[i]] <- cbind(coef(Recall(data, metric, sqrt.phy, traitgram=traitgram[i], which.eigen=which.eigen, remove.errors=remove.errors, traitgram.p=traitgram.p, ...)), traitgram[i], sites(data))
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
      if(!inherits(ext.dist, "dist"))
          stop("'ext.dist' must be a distance matrix")
      if(attr(ext.dist, "Size") != ncol(data$comm))
          stop("'ext.dist' must have dimensions matching comparative.comm object's species'")
      if(!identical(attr(ext.dist, "Labels"), species(data)))
          warning("'ext.dist' names do not match species data; continuing regardless")
      dist <- as.matrix(ext.dist)
      ext.dist <- TRUE
  } else ext.dist <- FALSE
  
  #Setup
  nspp <- ncol(data$comm)
  nsite <- nrow(data$comm)
  SR <- rowSums(data$comm > 0)
  data$comm[data$comm > 0] <- 1
  coefs <- data.frame(row.names=rownames(data$comm))
  if(traitgram==FALSE & ext.dist==FALSE)
      dist <- cophenetic(data$phy)
  if(sqrt.phy){
        if(!is.ultrametric(data$phy))
            warning("Phylogeny is not ultrametric; see function details")
      dist <- sqrt(dist)
      data$phy <- as.phylo(hclust(as.dist(dist)))
  }
  #Remove missing species
  dist <- dist[colSums(data$comm)>0, colSums(data$comm)>0]
  data <- data[,colSums(data$comm)>0]
  output <- list(psv=NULL, psr=NULL, mpd=NULL, mntd=NULL, pd=NULL, pd.ivs=NULL, colless=NULL, gamma=NULL, taxon=NULL, eigen.sum=NULL, Eed=NULL, Hed=NULL, dist.fd=NULL)
  
  #Calculate measures
  if((metric == "psv" | metric == "all" | metric == "all-quick") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE))
    output$psv <- coefs$psv <- try(psd(data$comm, data$phy, ...)[,1], silent = TRUE)
  
  if((metric == "psr" | metric == "all" | metric == "all-quick") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE))
    output$psr <- coefs$psr <- try(psd(data$comm, data$phy, ...)[,4], silent = TRUE)
  
  if(metric == "mpd" | metric == "all" | metric == "all-quick")
    output$mpd <- coefs$mpd <- try(mpd(data$comm, dist, abundance.weighted=FALSE, ...), silent = TRUE)
  
  if((metric == "pd" | metric == "all" | metric == "all-quick") & traitgram==FALSE & ext.dist==FALSE)
    try({
      output$pd <- coefs$pd <- pd(data$comm, data$phy, ...)[,1]
      output$pd.ivs <- coefs$pd.ivs <- unname(resid(lm(coefs$pd ~ rowSums(data$comm))))
    }, silent = TRUE)

  if(metric == "mntd" | metric == "all" | metric == "all-quick")
    try(output$mntd <- coefs$mntd <- mntd(data$comm, dist, abundance.weighted=FALSE, ...), silent = TRUE)
  
  if((metric == "colless" | metric == "all" | metric == "all-quick") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE))
      try(output$colless <- coefs$colless <- .colless(data), silent = TRUE)
  
  
  if((metric == "gamma" | metric == "all" | metric == "all-quick") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE))
      try({
          tree.shape <- as.treeshape(data$phy)
          nams <- tree.shape$names
          output$gamma <- coefs$gamma <- apply(data$comm, 1, .gamma, data$phy, nams)
      }, silent = TRUE)
  
  #Note - I've cut out some here because I think the simplification
  #can happen in the summary output
  if(metric == "taxon" | metric == "all" | metric == "all-quick")
      try({
        output$taxon <- taxondive(data$comm, dist, ...)
        t <- data.frame(output$taxon$D, output$taxon$Dstar, output$taxon$Lambda, output$taxon$Dplus, output$taxon$SDplus)
        names(t) <- c("Delta", "DeltaStar", "LambdaPlus", "DeltaPlus", "S.DeltaPlus")
        coefs <- cbind(coefs, t)
      }, silent = TRUE)
  
  if(metric == "eigen.sum" | metric == "all" | metric == "all-quick")
      try({
              #PVRdecomp ignores a phylogeny if given a distance matrix; suppress its warning
          t <- options("warn")
          options(warn=-10)
          evc <- PVRdecomp(data$phy, dist=dist, ...)@Eigen$vectors
          options(warn=as.numeric(t))
          output$eigen.sum <- coefs$eigen.sum <- apply(data$comm, 1, .eigen.sum, evc, which.eigen)
      }, silent = TRUE)
  
  if((metric == "eed" | metric == "all" | metric == "all-quick") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE))
      try(output$Eed <- coefs$Eed <- .eed(data), silent=TRUE)

  if((metric == "hed" | metric == "all" | metric == "all-quick") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE))
      try(output$Hed <- coefs$Hed <- .hed(data), silent=TRUE)

  if(metric == "dist.fd" | metric == "all")
      try({
          if(!is.null(data$data) | traitgram==FALSE | ext.dist==FALSE){
              t <- dist
          } else t <- data$data
          output$dist.fd$output <- capture.output(output$dist.fd <- dbFD(t, data$comm, w.abun=FALSE, messages=TRUE, ...), file=NULL)
          coefs <- with(output$dist.fd, cbind(coefs, cbind(FRic, FEve, FDiv, FDis, RaoQ)))
          #Only bother getting CWMs if we have trait data
          if(!is.null(data$data)){
              t <- output$dist.fd$CWM
              colnames(t) <- paste(colnames(t), "cmw", sep=".")
              coefs$dist.fd <- rbind(coef$dist.fd, t)
          }
      }, silent=TRUE)
    
  if(remove.errors) output <- lapply(output, .removeErrors)
  
  #Prepare output
  output$type <- "shape"
  output$coefs <- coefs
  class(output) <- "phy.structure"
  return(output)
}

#Internal Colless function
#' @importFrom apTreeshape colless tipsubtree
.colless <- function(data)
{
    output <- numeric(nrow(data$comm))
    for(i in seq(nrow(data$comm)))
        output[i] <- colless(as.treeshape(drop_tip(data$phy, colnames(data$comm)[data$comm[i,]==0])))
    names(output) <- rownames(data$comm)
    return(output)
}

#' @importFrom ape gammaStat
.gamma<-function(pa.vec,tree,nams){
    if(sum(pa.vec)<3){
        return(NA)
    } else {
        if(length(setdiff(tree$tip.label, nams)) != 0){
            tree <- drop.tip(setdiff(tree$tip.label, ))
        }
        return(gammaStat(drop.tip(tree,nams[pa.vec==0])))
    }
}

#' @importFrom picante pd
#' @importFrom picante evol.distinct
.hed <- function(data){
    #Argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")

    #Setup
    ed <- evol.distinct(data$phy, "fair.proportion")$w
    pd <- pd(data$comm, data$phy)$PD

    #Internal assemblage calc.
    ..hed <- function(ed, pd.comm){
        hed <- ed / pd.comm
        return(-sum(hed * log(hed)))
    }

    #Calculate, clean, and return
    output <- numeric(nrow(data$comm))
    names(output) <- rownames(data$comm)
    for(i in seq(nrow(data$comm)))
        output[i] <- ..hed(ed, pd[i])
    return(output)
}

.eed <- function(data, na.rm=TRUE) {
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    output <- .hed(data) / log(apply(data$comm, 1, function(x) sum(x != 0)))
    names(output) <- rownames(data)
    return(output)
}

.eigen.sum <- function(x, evc, vecnums) {
    if(sum(x>0)) {
        return(sum(apply(as.matrix(evc[x>0,vecnums]),2,var)))
    } else {
        return(NA)
    }
}
