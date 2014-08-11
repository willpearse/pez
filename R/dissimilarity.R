## TODO: Carry through arguments (?)  as.dist argument would allow
## hclust-ing, etc., of these measures more easily are there
## implications from pa=FALSE?  Will has altered UniFrac and PhyloSor
## to remove the rowMeans (is that OK?)
#
#
#
#' Calculate dissimilarity phylogenetic biodiversity metrics across communities
#' 
#' \code{dissimilarity} calculates phylogenetic biodiversity metrics
#' 
#' @param data a \code{comparative.comm} object
#' @param metric specify particular metrics to calculate, default is \code{all}
#' @param pa If TRUE (default), all metrics are calculated across presence-absence matrices, and species abundances are ignored
#' @param permute Number of permutations for metric (currently only for PCD)
#' @details Calculates various metrics of phylogenetic biodiversity that are categorized as \emph{dissimilarity} metrics by Pearse \emph{et al.} (2014)
#' @note This function uses a version of the PCD function, that is not included in \code{picante} and can be slow if \code{metric}="all"
#' @return a \code{phy.structure} list object of metric values
#' @author M.R. Helmus, Will Pearse
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
#' @examples \dontrun{
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' dissimilarity(data)
#' dispersion(data, "pcd")
#' }
#' @importFrom picante unifrac phylosor
#' @export
dissimilarity <- function(data, metric=c("all", "unifrac", "pcd", "phylosor"), pa=TRUE, permute=100)
{

    
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  if(permute < 0) stop("Can't have negative null permutations!")
  
  #Setup
  if(pa)
    data$comm <- data$comm[rowSums(data$comm)>0, ]
  output <- list(unifrac=NULL, pcd=NULL, phylosor=NULL)
  
  #Caculate measures
  if(metric == "unifrac" | metric == "all")
    output$unifrac <- unifrac(data$comm, data$phy)
  
  if(metric == "pcd" | metric == "all")
    output$pcd <- PCD(data$comm, data$phy, reps=permute)

  if(metric == "phylosor" | metric == "all")
    output$phylosor <- phylosor(data$comm, data$phy)
  
  #Prepare output
  output$type <- "dissimilarity"
  class(output) <- "phy.structure"
  return(output)
}
