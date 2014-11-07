#' Calculate (phylogenetic) dissimilarity: compare assemblages to
#' one-another
#'
#' As described in Pearse et al. (2014), a dissimilarity metric
#' compares diversity between communities. WARNING: Phylosor is
#' presented as a distance matrix here, i.e. it is *not* the fraction
#' of shared branch length among communities, but rather '1 - shared
#' branch length'. This means \code{dissimilarity} always returns a
#' *distance* object, not a similarity object; this is a different
#' convention from other packages.
#'
#' Using square-rooted distance matrices, or distance matrices that
#' incorporate trait information, can be an excellent thing to do, but
#' (for the above reasons), \code{pez} won't give you an answer for
#' metrics for which WDP thinks it makes no sense. All results from
#' this other than \code{comdist} *will always be wrong* if you do not
#' have an ultrametric tree and square-root (branch lengths
#' proportional to time) and you will be warned about this. WDP
#' strongly feels you should only be using ultrametric phylogenies in
#' any case, but code to fix this bug is welcome.
#' 
#' @param data \code{comparative.comm} object
#' @param metric default (\code{all}) calculates everything;
#' individually call-able metrics are: \code{unifrac}, \code{pcd},
#' \code{phylosor}, \code{comdist}.
#' @param abundance If TRUE (default) metrics are calculated
#' incorporating species abundances (currently only comdist)
#' @param permute Number of permutations for metric (currently only
#' for \code{pcd})
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
#' @param ext.dist Supply an external species-level distance matrix
#' for use in calculations. See `details' for comments on the use of
#' distance matrices in different metric calculations.
#' @param ... additional parameters to be passed to `metric
#' function(s) you are calling
#' @return list object of metric values. A \code{coef} method does not
#' exist for this function, because there's no nice way to simplify
#' all the distance matrices. Sorry!
#' @author M.R. Helmus, Will Pearse
#' @seealso \code{\link{shape}} \code{\link{evenness}} \code{\link{dispersion}}
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus
#' M.R. (2014). Metrics and Models of Community Phylogenetics. In:
#' Modern Phylogenetic Comparative Methods and Their Application in
#' Evolutionary Biology. Springer Berlin Heidelberg, pp. 451-464.
#' @references \code{unifrac} Lozupone C.A. & Knight
#' R. (2005). UniFrac: a new phylogenetic method for comparing
#' microbial communities. Applied and Environmental Microbiology, 71,
#' 8228-8235.
#' @references \code{pcd} Ives A.R. & Helmus M.R. (2010). Phylogenetic
#' metrics of community similarity. The American Naturalist, 176,
#' E128-E142.
#' @references \code{phylosor} Bryant J.A., Lamanna C., Morlon H.,
#' Kerkhoff A.J., Enquist B.J. & Green J.L. (2008). Microbes on
#' mountainsides: Contrasting elevational patterns of bacterial and
#' plant diversity. Proceedings of the National Academy of Sciences of
#' the United States of America, 105, 11505-11511.
#' @references \code{comdist} C.O. Webb, D.D. Ackerly, and
#' S.W. Kembel. 2008. Phylocom: software for the analysis of
#' phylogenetic community structure and trait
#' evolution. Bioinformatics 18:2098-2100.
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' \dontrun{dissimilarity(data)}
#' dissimilarity(data, "unifrac")
#' @importFrom picante unifrac phylosor pcd comdist
#' @importFrom ape is.ultrametric as.phylo
#' @export
dissimilarity <- function(data, metric=c("all", "unifrac", "pcd", "phylosor", "comdist"), abundance=TRUE, permute=100, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, ...)
{   
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  if(permute < 0) stop("Can't have negative null permutations!")

  if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
      stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
  if(!is.null(traitgram)){
      if(length(traitgram) > 1)
          stop("Cannot calculate dissimilarities en mass for traitgram. See help file.") else {
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
  output <- list(unifrac=NULL, pcd=NULL, phylosor=NULL, comdist=NULL)
  if(traitgram==FALSE & ext.dist==FALSE)
      dist <- cophenetic(data$phy)
  if(sqrt.phy){
        if(!is.ultrametric(data$phy))
            warning("Phylogeny is not ultrametric; see function details")
      dist <- sqrt(dist)
      data$phy <- as.phylo(hclust(as.dist(dist)))
  }
  if(abundance == FALSE)
      data$comm[data$comm > 1] <- 1
  #Caculate measures
  if((metric == "unifrac" | metric == "all") & (ext.dist==FALSE & traitgram==FALSE))
    output$unifrac <- unifrac(data$comm, data$phy, ...)
  
  if((metric == "pcd" | metric == "all") & (ext.dist==FALSE & traitgram==FALSE))
    output$pcd <- pcd(data$comm, data$phy, reps=permute, ...)

  #NOTE: I'm flipping phylosor to be a distance matrix
  if((metric == "phylosor" | metric == "all")  & (ext.dist==FALSE & traitgram==FALSE)){
    output$phylosor <- phylosor(data$comm, data$phy, ...)
    output$phylosor <- as.dist(1 - as.matrix(output$phylosor))
  }

  if(metric == "comdist" | metric == "all")
    output$comdist <- comdist(data$comm, dist, abundance.weighted=abundance, ...)
  
  #Prepare output
  output$type <- "dissimilarity"
  class(output) <- "phy.structure"
  return(output)
}
