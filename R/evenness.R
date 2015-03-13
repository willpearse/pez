#' Calculate (phylogenetic) evenness: examine assemblage composition
#' and abundance
#'
#'
#' As described in Pearse et al. (2014), an evenness metric is one the
#' examines the phylogenetic structure of species present in each
#' assemblage, taking into account their abundances. For completeness,
#' options are provided to calculate these metrics using species
#' traits.
#'
#' Most of these metrics do not involve comparison with some kind of
#' evolutionary-derived expectation for phylogenetic shape. Those that
#' do, however, such as PSE, make no sense unless applied to a
#' phylogenetic distance matrix - their null expectation *requires*
#' it. Using square-rooted distance matrices, or distance matrices
#' that incorporate trait information, can be an excellent thing to
#' do, but (for the above reasons), \code{pez} won't give you an
#' answer for metrics for which WDP thinks it makes no
#' sense. \code{pae}, \code{iac}, \code{haead} & \code{eaed} can
#' (...up to you whether you should!...)  be used with a square-rooted
#' distance matrix, but the results *will always be wrong* if you do
#' not have an ultrametric tree (branch lengths proportional to time)
#' and you will be warned about this. WDP strongly feels you should
#' only be using ultrametric phylogenies in any case, but code to fix
#' this bug is welcome.
#' @param data \code{\link{comparative.comm}} object
#' @param metric default (\code{all-quick}) calculates everything bar
#' \code{fd.dist} and the Pagel transformations
#' (\eqn{$\lambda$}{lambda}, \eqn{$\delta$}{delta},
#' \eqn{$\kappa$}{kappa}). \code{all} calculates
#' everything. Individually call-able metrics are: \code{rao},
#' \code{taxon}, \code{entropy}, \code{pae}, \code{iac}, \code{haed},
#' \code{eaed}, \code{lambda}, \code{kappa}, \code{delta}, \code{mpd},
#' \code{mntd}, \code{pse}, & \code{dist.fd}. Note that \code{dist.fd}
#' is a trait distance metric, but will be calculated on the
#' phylogenetic distance matrix if not traits are supplied, and the
#' external/square-rooted phylogenetic distance matrix if either is
#' specified.
#' @param sqrt.phy If TRUE (default is FALSE) your phylogenetic
#' distance matrix will be square-rooted; specifying TRUE will force
#' the square-root transformation on phylogenetic distance matrices
#' (in the spirit of Leitten and Cornwell, 2014). See `details' for
#' details about different metric calculations when a distance matrix
#' is used.
#' @param traitgram If not NULL (default), a number to be passed to
#' \code{funct.phylo.dist} (\code{phyloWeight}; the `a' parameter),
#' causing analysis on a distance matrix reflecting both traits and
#' phylogeny (0-->only phylogeny, 1--> only traits; see
#' \code{funct.phylo.dist}). If a vector of numbers is given,
#' \code{shape} iterates across them and returns a \code{data.frame}
#' with coefficients from each iteration. See `details' for details
#' about different metric calculations when a distance matrix is used.
#' @param traitgram.p A value for `p' to be used in conjunction with
#' \code{traitgram} when calling \code{funct.phylo.dist}.
#' @param ext.dist Supply an external species-level distance matrix
#' for use in calculations. See `details' for comments on the use of
#' distance matrices in different metric calculations.
#' @param ... Additional arguments to passed to metric functions
#' (unlikely you will want this!)
#' @note As mentioned above, \code{dist.fd} is calculated using a
#' phylogenetic distance matrix if no trait data are available, or if
#' you specify \code{sqrt.phy}. It is not calculated by default
#' because it generates warning messsages (which WDP is loathe to
#' suppress) which are related to the general tendency for a low rank
#' of phylogenetic distance matrices. Much ink has been written about
#' this, and in par this problem is why the \code{eigen.sum} measure
#' came to be suggested.
#'
#' Some of these metrics can cause (inconsequential) warnings if given
#' assemblages with only one species/individual in them, and return
#' NA/NaN values depending on the metric. I consider these `features',
#' not bugs.
#' @return \code{phy.structure} list object of metric values. Use
#' \code{coefs} to extract a summary metric table, or examine each
#' individual metric (which gives more details for each) by calling
#' \code{print} on the output (i.e., type \code{output} in the example
#' below).
#' @author M.R. Helmus, Will Pearse
#' @seealso shape dispersion dissimilarity
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus
#' M.R. (2014). Metrics and Models of Community Phylogenetics. In:
#' Modern Phylogenetic Comparative Methods and Their Application in
#' Evolutionary Biology. Springer Berlin Heidelberg, pp. 451-464.
#' @references \code{pse} Helmus M.R., Bland T.J., Williams C.K. &
#' Ives A.R. (2007). Phylogenetic measures of biodiversity. American
#' Naturalist, 169, E68-E83.
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus M.R. (2014). Metrics and Models of Community Phylogenetics. In: Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Springer Berlin Heidelberg, pp. 451-464.
#' @references \code{pse} Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007). Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83.
#' @references \code{rao} Webb C.O. (2000). Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. American Naturalist, 156, 145-155.
#' @references \code{taxon} Clarke K.R. & Warwick R.M. (1998). A taxonomic distinctness index and its statistical properties. J. Appl. Ecol., 35, 523-531.
#' @references \code{entropy} Allen B., Kon M. & Bar-Yam Y. (2009). A New Phylogenetic Diversity Measure Generalizing the Shannon Index and Its Application to Phyllostomid Bats. The American Naturalist, 174, 236-243.
#' @references \code{pae,iac,haed,eaed} Cadotte M.W., Davies T.J., Regetz J., Kembel S.W., Cleland E. & Oakley T.H. (2010). Phylogenetic diversity metrics for ecological communities: integrating species richness, abundance and evolutionary history. Ecology Letters, 13, 96-105.
#' @references \code{lambda,delta,kappa} Mark Pagel (1999) Inferring the historical patterns of biological evolution. Nature 6756(401): 877--884.
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' evenness(data)
#' evenness(data, "rao")
#' @importFrom ape cophenetic.phylo drop.tip is.ultrametric as.phylo
#' @importFrom picante pse raoD
#' @importFrom vegan taxondive
#' @importFrom caper comparative.data pgls summary.pgls coef.pgls
#' @importFrom ade4 newick2phylog
#' @importFrom FD dbFD
#' @export
evenness <- function(data, metric=c("all-quick", "all", "rao", "taxon", "entropy", "pae", "iac", "haed", "eaed", "lambda", "delta", "kappa", "mpd", "mntd", "pse", "dist.fd"), sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, ...)
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
              output[[i]] <- cbind(coef(Recall(data, metric, sqrt.phy, traitgram=traitgram[i], traitgram.p=traitgram.p, ...)), traitgram[i], sites(data))
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
  SR <- rowSums(data$comm>0)
  nsite <- nrow(data$comm)
  nspp <- ncol(data$comm)
  coefs <- data.frame(row.names=rownames(data$comm))
  if(sqrt.phy)
      data <- .sqrt.phy(data)
  if(traitgram==FALSE & ext.dist==FALSE)
      dist <- cophenetic(data$phy)

  #Remove missing species
  dist <- dist[colSums(data$comm)>0, colSums(data$comm)>0]
  data <- data[,colSums(data$comm)>0]
  output <- list(pse=NULL, rao=NULL, taxon=NULL, entropy=NULL, cadotte=NULL, delta=NULL, mpd=NULL, mntd=NULL, dist.fd=NULL)
  
  #Filter metrics according to suitability and calculate
  functions <- setNames(c(.rao, .phylo.entropy, .pae, .iac, .haed, .eaed, .lambda, .delta, .kappa, .mpd, .mntd, .pse, .dist.fd), c("rao", "entropy", "pae", "iac", "haed", "eaed", "lambda", "delta", "kappa", "mpd", "mntd", "pse", "dist.fd"))
  if(metric == "all-quick")
      functions <- functions[!names(functions) %in% c("dist.fd","lambda","delta","kappa")]
  if(sqrt.phy == TRUE)
      functions <- functions[!names(functions) %in% c("pse","lambda","delta","kappa")]
  if(traitgram == TRUE)
      functions <- functions[!names(functions) %in% c("pse", "entropy", "iac", "haed", "lambda", "delta", "kappa", "rao")]
  if(ext.dist == TRUE)
      functions <- functions[!names(functions) %in% c("pse", "rao", "entropy", "pae", "iac", "haed", "eaed", "lambda", "delta", "kappa")]
  output <- lapply(functions, function(x) try(x(data, dist=dist, abundance.weighted=TRUE), silent=TRUE))

  #Clean up output and return
  output <- Filter(function(x) !inherits(x, "try-error"), output)
  output <- do.call(cbind, output)
  return(output)
}
