#TODO:
# Are Delta et al., when calculated here, based on abundances? Worth checking...
# There is an error message coming from phylo.entropy that is described in its script
# Could encapsulate phylogenetic signal measures by writing a wrapper and using expression to wrap the correct thing in (if you cared enough...)
# optimise (? and if can be bothered)
# mis-match between proportions and partial.tree$edge.length - not introduced by the changes below
# NEEDED lambda, delta, kappa reference

#' Calculate evenness phylogenetic biodiversity metrics across communities
#'
#' \code{evenness} calculates phylogenetic biodiversity metrics
#' 
#' @param data a comparative community ecology object
#' @param metric specify particular metrics to calculate, default is \code{all}
#' @details Calculates various metrics of phylogenetic biodiversity that are categorized as \emph{evenness} metrics by Pearse \emph{et al.} (2014)
#' @return a \code{phy.structure} list object of metric values
#' @author M.R. Helmus, Will Pearse
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
#' @references \code{cadotte} (i.e., \emph{PAE, IAC, Haed, Eaed}) Cadotte M.W., Davies T.J., Regetz J., Kembel S.W., Cleland E. & Oakley T.H. (2010). Phylogenetic diversity metrics for ecological communities: integrating species richness, abundance and evolutionary history. Ecology Letters, 13, 96-105.
#' @references \code{pst} Hardy O.J. & Senterre B. (2007). Characterizing the phylogenetic structure of communities by an additive partitioning of phylogenetic diversity. Journal of Ecology, 95, 493-506.
#' @references \code{lambda}
#' @references \code{delta}
#' @references \code{kappa}
#' @examples \dontrun{
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' evenness(data)
#' evenness(data, "rao")
#' }
#' @importFrom ape cophenetic.phylo drop.tip
#' @importFrom picante pse raoD
#' @importFrom vegan taxondive
#' @importFrom caper comparative.data pgls summary.pgls coef.pgls
#' @importFrom ade4 newick2phylog
#' @export
evenness <- function(data, metric=c("all", "rao", "taxon", "entropy", "cadotte", "pst", "lambda", "delta", "kappa"))
{
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  metric <- match.arg(metric)
  
  #Setup
  if(is.null(data$vcv))
    data$vcv <- cophenetic(data$phy)
  SR <- rowSums(data$comm>0)
  nsite <- nrow(data$comm)
  nspp <- ncol(data$comm)
  output <- list(pse=NULL, rao=NULL, taxon=NULL, entropy=NULL, cadotte=NULL, pst=NULL, delta=NULL)
  
  #Caculate measures
  if(metric == "pse" | metric == "all")
    output$pse <- pse(data$comm, data$phy)[,1]
  
  if(metric == "rao" | metric == "all")
    output$rao <- raoD(data$comm, data$phy)$Dkk
  
  if(metric == "taxon" | metric == "all"){
    output$taxon <- taxondive(data$comm, data$vcv)
  }
    
  
  if(metric == "entropy" | metric == "all")
    output$entropy <- .phylo.entropy(data)
  
  if(metric == "cadotte" | metric == "all"){
    temp <- data.frame(PAE=.pae(data), IAC=.iac(data), Haed=.haed(data), Eaed=.haed(data)/log(rowSums(data$comm)))
    output$cadotte <- temp[match(rownames(data$comm), rownames(temp)),]
  }
  
  if(metric == "pst" | metric == "all"){
    .abund <- data$comm[rowSums(data$comm>0)>1, ]
    if(length(setdiff(data$phy$tip.label, colnames(.abund))))
      tree <- drop.tip(data$phy, setdiff(data$phy$tip.label, colnames(.abund))) else tree <- data$phy
    temp <- .simpson.phylogenetic(data)
    output$pst <- temp[match(rownames(data$comm), names(temp))]
  }
  
  if(metric == "lambda" | metric == "all"){
    output$lambda <- list(models=vector("list", nrow(data$comm)), values=numeric(nrow(data$comm)))
    for(i in seq(nrow(data$comm))){
      c.data <- data$comm[i,][data$comm[i,] > 0]
      if(length(unique(c.data)) > 1){
        c.data <- data.frame(comm=c.data, names=names(c.data))
        c.data <- comparative.data(phy=drop_tip(data$phy, setdiff(data$phy$tip.label, c.data$names)), data=c.data, names.col=names)
        output$lambda$models[[i]] <- pgls(comm ~ 1, c.data, lambda="ML")
        output$lambda$values[i] <- summary(output$lambda$models[[i]])$param.CI$lambda$opt
      } else {
        output$lambda$models[[i]] <- NA
        output$lambda$values[i] <- NA
      }
    }
  }
  
  if(metric == "delta" | metric == "all"){
    output$delta <- list(models=vector("list", nrow(data$comm)), values=numeric(nrow(data$comm)))
    for(i in seq(nrow(data$comm))){
      c.data <- data$comm[i,][data$comm[i,] > 0]
      if(length(unique(c.data)) > 1){
        c.data <- data.frame(comm=c.data, names=names(c.data))
        c.data <- comparative.data(phy=drop_tip(data$phy, setdiff(data$phy$tip.label, c.data$names)), data=c.data, names.col=names)
        output$delta$models[[i]] <- pgls(comm ~ 1, c.data, delta="ML")
        output$delta$values[i] <- summary(output$delta$models[[i]])$param.CI$delta$opt
      } else {
        output$delta$models[[i]] <- NA
        output$delta$values[i] <- NA
      }
    }
  }
  
  if(metric == "kappa" | metric == "all"){
    output$kappa <- list(models=vector("list", nrow(data$comm)), values=numeric(nrow(data$comm)))
    for(i in seq(nrow(data$comm))){
      c.data <- data$comm[i,][data$comm[i,] > 0]
      if(length(unique(c.data)) > 1){
        c.data <- data.frame(comm=c.data, names=names(c.data))
        c.data <- comparative.data(phy=drop_tip(data$phy, setdiff(data$phy$tip.label, c.data$names)), data=c.data, names.col=names)
        output$kappa$models[[i]] <- pgls(comm ~ 1, c.data, kappa="ML")
        output$kappa$values[i] <- summary(output$kappa$models[[i]])$param.CI$kappa$opt
      } else {
        output$kappa$models[[i]] <- NA
        output$kappa$values[i] <- NA
      }
    }
  }
  
  #Prepare output
  output$type <- "evenness"
  class(output) <- "phy.structure"
  return(output)
}

#' @importFrom ape write.tree
#' @importFrom ade4 newick2phylog
.phylo.entropy <- function(data)
{
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  
  #Setup
  tree.phylog <- newick2phylog(write.tree(data$phy))
  species <- colnames(data$comm)
  hp.sites <- numeric(nrow(data$comm))
  
  ## beginning of the sites iteration
  for (i in seq(nrow(data$comm)))
  {
    ## species which occured at the site
    site.species <- species[which(data$comm[i,] > 0)]
    
    ## species which NOT occured at the site. They will be removed
    ## from the phylogenetic tree.
    other.species <- species[which(data$comm[i,] == 0)]
    
    ## proportions of occurrence of each species
    proportions <- data$comm[i,site.species] / sum(data$comm[i,])
    
    ## god, how I hate these tree conversions...
    ## also, this comparison is very ugly, it has to be a better
    ## way...
    ## TODO: Search for a better way to verify if an object is empty
    if (length(other.species) == 0) {
      partial.tree <- tree.phylog
    } else {
      if(length(site.species) == 1) {other.species<-c(other.species,site.species)}
      partial.tree <- drop.tip(data$phy, other.species)
      if (all(partial.tree$edge.length[1] == partial.tree$edge.length) | length(site.species) == 2)
      {
        #SOMETHING WRONG HERE
        hp.sites[rownames(data$comm)[i]] <- abs(sum(proportions * log(proportions) * partial.tree$edge.length))
        next
      }
      partial.tree <- newick2phylog(write.tree(drop.tip(data$phy, other.species)))
    }
    ## TODO: Some (I think) partial trees are not rooted. The results
    ## seem to be ok, but the paper states that Hp should be
    ## calculated from a rooted tree. Do not know why though. Check
    ## what can be done to keep partial trees rooted.
    if ("Root" %in% names(partial.tree$nodes))
    {
      partial.branches <-
        partial.tree$nodes[-c(length(partial.tree$nodes))]
    } else {
      partial.branches <- partial.tree$nodes
    }
    ## terminal branches sizes
    partial.leaves <- partial.tree$leaves
    
    ## first part of the calculations. Here we calculate the index for
    ## each terminal branch
    sum.leaves <- sum(partial.leaves *
                        proportions[names(partial.leaves)] *
                        log(proportions[names(partial.leaves)]))
    
    ## storing the first part of the calculation
    hp <- c(sum.leaves)
    
    ## initilizing the list that will hold the descending leaves for
    ## each branch
    descending.leaves <- list()
    
    ## determining the descending leaves for each branch
    for (j in names(partial.branches)) {
      if (all(partial.tree$parts[[j]] %in% names(partial.leaves))) {
        descending.leaves[[j]] <- partial.tree$parts[[j]]
      } else {
        branches <- partial.tree$parts[[j]][!partial.tree$parts[[j]]
                                            %in% names(partial.leaves)]
        leaves <- partial.tree$parts[[j]][partial.tree$parts[[j]] %in%
                                            names(partial.leaves)]
        for (k in branches) {
          leaves <- c(leaves, descending.leaves[[k]])
        }
        descending.leaves[[j]] <- leaves
      }
    }
    ## calculating the index for each internal branch
    for (j in names(partial.branches)) {
      sum.proportions.desc.leaves <-
        sum(proportions[descending.leaves[[j]]])
      hp <- c(hp, (partial.branches[[j]] * sum.proportions.desc.leaves
                   * log(sum.proportions.desc.leaves)))
    }
    ## putting it all together
    hp.sites[rownames(data$comm)[i]] <- abs(sum(hp))
  }
  #Make the 0 values NAs
  hp.sites[hp.sites==0]<-NA
  ## the end.
  return(hp.sites)
}

#' @importFrom caper clade.matrix
.aed <- function(data, na.rm=TRUE) {
    #Argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    
    #Internal functions
    .desc <- function(tree) {
        t <- clade.matrix(tree)
        mat <- t$clade.matrix
        mat <- mat[!apply(mat, 1, function(x) all(x==1)), ]
        mat <- apply(mat, 1, as.logical)
        rownames(mat) <- tree$tip.label
        return(list(mat=mat, edge=t$edge.length))
    }
    aed <- function(i){
        spp <- row.names(desc[[i]]$mat)
        dAbund <- data$comm[i, spp] * desc[[i]]$mat
        AED <- colSums(edge.length[[i]] * t(prop.table(dAbund, margin=2)), na.rm=TRUE)
        AED/data$comm[i, spp]
    }
    
    subtrees <- assemblage.phylogenies(data)
    desc <- lapply(subtrees, .desc)
    edge.length <- lapply(desc, function(x) x$edge[as.numeric(colnames(x$mat))])
    
    res <- lapply(seq(nrow(data$comm)), aed)
    names(res) <- rownames(data$comm)
    return(res)
}

#' @importFrom picante pd
.haed <- function(data, na.rm=TRUE) {
    #Argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    
    # Recast AED in terms of individuals
    AED <- .aed(data)
    PD <- pd(data$comm, data$phy)
    scaled.AED <- lapply(seq(nrow(data$comm)), function(i) {
        spp <- names(AED[[i]])        
        rep(unname(AED[[i]]), sum(data$comm[i, spp])) / PD$PD[i]
    })
    res <- sapply(scaled.AED, function(x) -sum(x * log(x)))
    names(res) <- names(AED)
    return(res)
}

#' @importFrom ape cophenetic.phylo 
.simpson.phylogenetic <- function(data) {
    N.relative <- prop.table(data$comm, 2)
    dmat <- cophenetic(data$phy)
    out <- apply(N.relative, 1, function(n) sum((n %o% n)*dmat))
    return(out) 
}

.iac <- function(data, na.rm=TRUE) {
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    
    subtrees <- assemblage.phylogenies(data)
    .ancestors <- function(tree, no.root=TRUE){
        mat <- clade.matrix(tree)$clade.matrix
        if(no.root)
            mat <- mat[!apply(mat, 1, function(x) all(x==1)), ]
        #Stop self-references
        diag(mat) <- 0
        apply(mat, 2, function(x) which(x == 1))
    }
        
    .denom <-  function(tree) {
        # Count number of lineages originating at each internal node
        # (i.e. number of splits)
        nSplits <- table(tree$edge[,1])
        # For each tip, take the product of the number of splits across
        # all of its ancestral nodes
        res <- sapply(.ancestors(tree), function(x)
            prod(nSplits[as.character(x)]))
        return(res * 2)
    }

    # now for each subtree...
    denom <- lapply(subtrees, .denom)
    denom <- do.call("cbind", denom)
    nnodes <- sapply(subtrees, function(x) x$Nnode)

    # Calculate expected number of individuals under null hypothesis
    # of equal allocation to each lineage at each (node) split 
    expected <- rowSums(data$comm, na.rm=na.rm) / t(denom)

    # IAC: summed absolute difference between expected and observed
    # abundances, divided by number of nodes
    return(rowSums(abs(expected - data$comm), na.rm=na.rm) / nnodes)
}

.pae <- function(data, na.rm=TRUE) {
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    
    subtrees <- assemblage.phylogenies(data)
    PD <- pd(data$comm, data$phy)
    tmp <- setNames(rep(0, ncol(data$comm)), colnames(data$comm))
    TL <- lapply(subtrees, function(tree) {
        #Get terminal edge length
        res <- data$phy$edge.length[data$phy$edge[,2] <= length(data$phy$tip.label)]
        tmp[match(data$phy$tip.label, names(tmp))] <- res
        tmp
    })
    TL <- do.call("cbind", TL)
    numer <- PD$PD + colSums(TL * (t(data$comm) - 1))
    denom <- PD$PD + (rowSums(data$comm, na.rm = na.rm) / rowSums(data$comm,
        na.rm=na.rm) - 1) * colSums(TL)
    res <- numer/denom
    names(res) <- rownames(data$comm)
    return(res)
}
