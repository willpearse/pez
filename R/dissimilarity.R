#TODO:
# Carry through arguments (?)
# as.dist argument would allow hclust-ing, etc., of these measures more easily
# are there implications from pa=FALSE?
# Will has altered UniFrac and PhyloSor to remove the rowMeans (is that OK?)

#' Calculate dissimilarity metrics across communities
#' 
#' \code{dissimilarity} calculates dissimilarity metrics in comparative.comm communities
#' 
#' @param data a comparative community ecology object
#' @param metric specify (a) particular metric(s) to calculate (unifrac, pcd, phylosor), or the default 'all'
#' @param pa If TRUE (default), all metrics are calculated across presence-absence matrices, and species abundances are ignored
#' @param permute Number of permutations for metric (currently only for PCD)
#' @details This calculates UniFrac, PCD (and its subcomponents), and the phylogenetic Sorenson's index, all defined as dissimilarity metrics in Pearse et al.
#' @note This function uses XXX's version of the PCD function, not that included in picante
#' @return cc.dissimilarity object (a named list with the output from each metric)
#' @author Matt Helmus, Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample)
#' dissimilarity(data)
#' dispersion(data, "pcd")
#' }
#' @import picante
#' @export

dissimilarity <- function(data, metric=c("all", "unifrac", "pcd", "phylosor"), pa=TRUE, permute=1000)
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
  class(output) <- c("cc.dissimilarity", class(output))
  return(output)
}

#' Print a dissimilarity object (...by summarising it...)
#' @method print cc.dissimilarity
#' @S3method print cc.dissimilarity
#' @export
print.cc.dissimilarity <- function(x, ...){
  summary(x)
}

#' Summarise a dissimilarity object
#' @method summary cc.dissimilarity
#' @S3method summary cc.dissimilarity
#' @export
summary.cc.dissimilarity <- function(object, ...){
  cat("\nDissimilarity metrics in this object:\n")
  if(!is.null(object$unifrac)){
    cat("\tUniFrac (unifrac)\n")
  }
  if(!is.null(object$pcd)){
    cat("\tPCD (pcd)\n")
  }
  if(!is.null(object$phylosor)){
    cat("\tPhyloSor (phylosor)\n")
  }
  cat("Use something like 'output$unifrac' to work with each measure\n")
}
