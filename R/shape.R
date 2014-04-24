#TODO:
# I'm not sure PD should be standardised by the sum of species, but rather the sum of individuals. In general, we need to have a better discussion about individual vs. species metrics and scaling
# incorporate cat with null models by doing temp<-cat;cat<-function(...){} to hide all the output from colless.test!
# it's going to be hard for people to install ecoPD, so we need some sensible 'require'-foo to get everything in there
# I don't like having to 'require' ape when using phylobase, but either I'm doing something really stupid with the dependencies or phylobase doesn't have a properly setup NAMESPACE. Grrrr.

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
#' shape(data)
#' shape(data, "colless")
#' shape(data, "eigen.sum", 2)
#' }
#' @importFrom picante psd mpd pd
#' @importFrom vegan taxondive
#' @importFrom PVR PVRdecomp
#' @importFrom apTreeshape as.treeshape as.treeshape.phylo colless tipsubtree
#' @importFrom ecoPD eed hed
#' @importFrom phylobase phylo4d
#' @importFrom ape gammaStat cophenetic.phylo drop.tip
#' @export
shape <- function(data, metric=c("all", "psv", "psr", "mpd", "pd", "colless", "gamma", "taxon", "eigen.sum", "cadotte.pd"), which.eigen=1)  #vecnums chooses the eigenvector to calculate sumvar in Diniz-Filho J.A.F., Cianciaruso M.V., Rangel T.F. & Bini L.M. (2011). Eigenvector estimation of phylogenetic and functional diversity. Functional Ecology, 25, 735-744.
{
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  
  #Setup
  nspp <- ncol(data$comm)
  nsite <- nrow(data$comm)
  SR <- rowSums(data$comm > 0)
  pa <- data$comm > 0
  if(is.null(data$vcv))
    data$vcv <- cophenetic(data$phy)
  output <- list(psv=NULL, psr=NULL, mpd=NULL, pd=NULL, pd.ivs=NULL, colless=NULL, gamma=NULL, taxon=NULL, eigen.sum=NULL, cadotte.pd=NULL)
  
  #Caculate measures
  if(metric == "psv" | metric == "all")
    output$psv <- psd(pa, data$phy)[,1]
  
  if(metric == "psr" | metric == "all")
    output$psr <- psd(pa, data$phy)[,4]
  
  if(metric == "mpd" | metric == "all")
    output$mpd <- mpd(pa, data$vcv)
  
  if(metric == "pd" | metric == "all"){
    output$pd <- pd(pa, data$phy)[,1]
    output$pd.ivs <- resid(lm(output$pd ~ rowSums(pa)))
  }
  
  if(metric == "colless" | metric == "all"){
    tree.shape <- as.treeshape(data$phy)
    nams <- tree.shape$names
    output$colless <- apply(pa, 1, .colless, tree.shape, nams)
  }
  
  if(metric == "gamma" | metric == "all"){
    tree.shape <- as.treeshape(data$phy)
    nams <- tree.shape$names
    output$gamma <- apply(pa, 1, .gamma, data$phy, nams)
  }
  
  #Note - I've cut out some here because I think the simplification can happen in the summary output
  if(metric == "taxon" | metric == "all")
    output$taxon <- taxondive(pa, data$vcv)
  
  if(metric == "eigen.sum" | metric == "all"){
      evc <- PVRdecomp(data$phy)@Eigen$vectors  
      output$eigen.sum <- apply(pa, 1, 
                                function(x, evc, vecnums) if(sum(x>0)) return(sum(apply(as.matrix(evc[x>0,vecnums]),2,var))) else return(NA), 
                                evc,which.eigen)
  }
  
  if(metric == "cadotte.pd" | metric == "all"){
    .pa <- pa[rowSums(pa>0)>1,]
    .tree <- phylo4d(data$phy, t(.pa))
    temp <- data.frame(Eed=eed(.tree), Hed=hed(.tree))
    output$cadotte.pd <- temp[match(rownames(data$comm), rownames(temp)),]
  }
  
  #Prepare output
  output$type <- "shape"
  class(output) <- "phy.structure"
  return(output)
}

#Internal Colless function
.colless <- function(pa.vec,tree,nams)
{
  if(sum(pa.vec)<3)
  {
    return(NA)
  } else {
    return(colless(tipsubtree(tree, nams[pa.vec!=0])))
  }
}

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
