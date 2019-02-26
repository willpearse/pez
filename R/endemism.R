#' Calculate (phylogenetic) endemism
#'
#' At present, only a small number of metrics, but we intend for this
#' to grow with time. Note that metrics that incorporate abundance are
#' mixed in with those that do not. Some of these metrics make sense
#' when used with probabilities, for example those derived from an
#' SDM; some do not. You will have to use your own judgement (as with
#' everything in science!).
#' @param data \code{\link{comparative.comm}} object
#' @param sqrt.phy If TRUE (default is FALSE) your phylogenetic
#'     distance matrix will be square-rooted; specifying TRUE will
#'     force the square-root transformation on phylogenetic distance
#'     matrices (in the spirit of Leitten and Cornwell, 2014). See
#'     `details' for details about different metric calculations when
#'     a distance matrix is used.
#' @return \code{data.frame} with metric values.
#' @author Will Pearse, Dan Rosauer
#' @seealso \code{\link{pez.shape}} \code{\link{pez.evenness}}
#'     \code{\link{pez.dispersion}} \code{\link{pez.dissimilarity}}
#' @references \code{BED} Cadotte, M. W., & Jonathan Davies,
#'     T. (2010). Rarest of the rare: advances in combining
#'     evolutionary distinctiveness and scarcity to inform
#'     conservation at biogeographical scales. Diversity and
#'     Distributions, 16(3), 376-385.
#' @references \code{PE} Rosauer, D. A. N., Laffan, S. W., Crisp,
#'     M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic
#'     endemism: a new approach for identifying geographical
#'     concentrations of evolutionary history. Molecular Ecology,
#'     18(19), 4061-4072.
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' (output<-pez.endemism(data))
#' @export
pez.endemism <- function(data, sqrt.phy=FALSE){
  #Assertions, argument handling, and setup
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  if(sqrt.phy)
      data <- .sqrt.phy(data)
  
  #Filter metrics according to suitability and calculate
  functions <- setNames(c(.pe, .bed), c("PE", "BED"))
  output <- lapply(functions, function(x) try(x(data)))
  
  #Clean up output and return
  output <- Filter(function(x) !inherits(x, "try-error"), output)
  output <- do.call(cbind, output)
  return(as.data.frame(output))
}
