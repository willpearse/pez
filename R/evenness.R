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
  output <- list(pse=NULL, rao=NULL, taxon=NULL, entropy=NULL, cadotte=NULL, delta=NULL, mpd=NULL, mntd=NULL, dist.fd=NULL)
  
  #Caculate measures
  if((metric == "pse" | metric == "all" | metric == "all-quick") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE))
    try(output$pse <- coefs$pse <- pse(data$comm, data$phy)[,1], silent=TRUE)
  
  if((metric == "rao" | metric == "all" | metric == "all-quick") & (ext.dist==FALSE))
    try(output$rao <- coefs$rao <- raoD(data$comm, data$phy)$Dkk, silent=TRUE)
  
  if(metric == "taxon" | metric == "all" | metric == "all-quick")
    try({
      output$taxon <- taxondive(data$comm, dist, ...)
      t <- data.frame(output$taxon$D, output$taxon$Dstar, output$taxon$Lambda, output$taxon$Dplus, output$taxon$SDplus)
      names(t) <- c("taxon.Delta", "taxon.DeltaStar", "taxon.LambdaPlus", "taxon.DeltaPlus", "taxon.S.DeltaPlus")
      coefs <- cbind(coefs, t)
    }, silent = TRUE)
  
  if((metric == "entropy" | metric == "all" | metric == "all-quick") & (traitgram==FALSE & ext.dist==FALSE))
    try({
      output$entropy <- .phylo.entropy(data)
      coefs$entropy <- output$entropy[!is.na(output$entropy)]
    }, silent=TRUE)

  if((metric == "pae" | metric == "all" | metric == "all-quick") & (traitgram==FALSE & ext.dist==FALSE))
      try(output$pae <- coefs$pae <- .pae(data), silent=TRUE)

  if((metric == "iac" | metric == "all" | metric == "all-quick") & (traitgram==FALSE & ext.dist==FALSE))
      try(output$iac <- coefs$iac <- .iac(data), silent=TRUE)

  if((metric == "haed" | metric == "all" | metric == "all-quick") & (traitgram==FALSE & ext.dist==FALSE))
      try(output$Haed <- coefs$Haed <- .haed(data), silent=TRUE)

  if((metric == "eaed" | metric == "all" | metric == "all-quick") & (traitgram==FALSE & ext.dist==FALSE))
      try(output$Eaed <- coefs$Eaed <- .haed(data)/log(rowSums(data$comm)), silent=TRUE)
    
  if((metric == "lambda" | metric == "all") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE)){
    output$lambda <- list(models=vector("list", nrow(data$comm)), values=numeric(nrow(data$comm)))
    for(i in seq(nrow(data$comm))){
      c.data <- data$comm[i,][data$comm[i,] > 0]
      if(length(unique(c.data)) > 1){
        c.data <- data.frame(comm=c.data, names=names(c.data))
        c.data <- comparative.data(phy=drop_tip(data$phy, setdiff(data$phy$tip.label, c.data$names)), data=c.data, names.col=names)
        output$lambda$models[[i]] <- pgls(comm ~ 1, c.data, lambda="ML", ...)
        output$lambda$values[i] <- summary(output$lambda$models[[i]])$param.CI$lambda$opt
      } else {
        output$lambda$models[[i]] <- NA
        output$lambda$values[i] <- NA
      }
    }
    coefs$lambda <- output$lambda$values
  }
  
  if((metric == "delta" | metric == "all") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE)){
    output$delta <- list(models=vector("list", nrow(data$comm)), values=numeric(nrow(data$comm)))
    for(i in seq(nrow(data$comm))){
      c.data <- data$comm[i,][data$comm[i,] > 0]
      if(length(unique(c.data)) > 1){
        c.data <- data.frame(comm=c.data, names=names(c.data))
        c.data <- comparative.data(phy=drop_tip(data$phy, setdiff(data$phy$tip.label, c.data$names)), data=c.data, names.col=names)
        output$delta$models[[i]] <- pgls(comm ~ 1, c.data, delta="ML", ...)
        output$delta$values[i] <- summary(output$delta$models[[i]])$param.CI$delta$opt
      } else {
        output$delta$models[[i]] <- NA
        output$delta$values[i] <- NA
      }
    }
    coefs$delta <- output$delta$values
  }

  if((metric == "kappa" | metric == "all") & (sqrt.phy==FALSE & traitgram==FALSE & ext.dist==FALSE)){
    output$kappa <- list(models=vector("list", nrow(data$comm)), values=numeric(nrow(data$comm)))
    for(i in seq(nrow(data$comm))){
      c.data <- data$comm[i,][data$comm[i,] > 0]
      if(length(unique(c.data)) > 1){
        c.data <- data.frame(comm=c.data, names=names(c.data))
        c.data <- comparative.data(phy=drop_tip(data$phy, setdiff(data$phy$tip.label, c.data$names)), data=c.data, names.col=names)
        output$kappa$models[[i]] <- pgls(comm ~ 1, c.data, kappa="ML", ...)
        output$kappa$values[i] <- summary(output$kappa$models[[i]])$param.CI$kappa$opt
      } else {
        output$kappa$models[[i]] <- NA
        output$kappa$values[i] <- NA
      }
    }
    coefs$kappa <- output$kappa$values
  }

  if(metric == "mpd" | metric == "all" | metric == "all-quick")
    output$mpd <- coefs$mpd <- try(mpd(data$comm, dist, abundance.weighted=TRUE, ...), silent = TRUE)

  if(metric == "mntd" | metric == "all" | metric == "all-quick")
    output$mntd <- coefs$mntd <- try(mntd(data$comm, dist, abundance.weighted=TRUE, ...), silent = TRUE)

    if(metric == "dist.fd" | metric == "all")
      try({
          if(!is.null(data$data) | traitgram==FALSE | ext.dist==FALSE){
              t <- dist
          } else t <- data$data
          output$dist.fd$output <- capture.output(output$dist.fd <- dbFD(t, data$comm, w.abun=TRUE, messages=TRUE, ...), file=NULL)
          coefs <- with(output$dist.fd, cbind(coefs, cbind(FRic, FEve, FDiv, FDis, RaoQ)))
          #Only bother getting CWMs if we have trait data
          if(!is.null(data$data)){
              t <- output$dist.fd$CWM
              colnames(t) <- paste(colnames(t), "cmw", sep=".")
              coefs$dist.fd <- rbind(coef$dist.fd, t)
          }
      }, silent=TRUE)
  
  #Prepare output
  output$type <- "evenness"
  output$coefs <- coefs
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
        hp.sites[i] <- abs(sum(proportions * log(proportions) * partial.tree$edge.length[1]))
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
    hp.sites[i] <- abs(sum(hp))
  }
  #Make the 0 values NAs
  hp.sites[hp.sites==0]<-NA
  names(hp.sites) <- rownames(data$comm)
  ## the end.
  return(hp.sites)
}

#' @importFrom ape extract.clade
#' @importFrom caper clade.matrix
# Much of this is simplified from the original because we can assume
#   the order of the components of a comparative.comm
.aed <- function(data){
    #Setup
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    uni.nodes <- data$phy$edge[,2][!data$phy$edge[,2] %in% seq_along(data$phy$tip.label)]

    #Internal AED for each assemblage
    ..aed <- function(assemblage, tree, uni.nodes, clade.matrix){
        #Nodal values
        node.values <- numeric(length(uni.nodes))
        for(i in seq_along(uni.nodes)){
            t <- extract.clade(tree, uni.nodes[i])
            t.abund <- assemblage[names(assemblage) %in% t$tip.label]
            node.values[i] <- (tree$edge.length[which(tree$edge[,2]==uni.nodes[i])]) / sum(t.abund)
        }
        
        #AED
        aed <- numeric(length(assemblage))
        for(i in seq_along(tree$tip.label)){
            sp <- tree$tip.label[i]
            nodes <- rownames(clade.matrix)[clade.matrix[,i] == 1]
            splength <- tree$edge.length[tree$edge[,2] == i]
            t <- assemblage[i]
            aed[i] <- sum(node.values[which(uni.nodes %in% nodes)]) + unname(ifelse(t==0,0,splength/t))
        }
        return(aed)
    }

    #Calculate, neaten, and return
    aed <- apply(data$comm, 1, ..aed, data$phy, uni.nodes, clade.matrix(data$phy)$clade.matrix[-seq_along(data$phy$tip.label),])
    rownames(aed) <- data$phy$tip.label
    colnames(aed) <- rownames(data$comm)
    aed[aed == Inf | aed == -Inf] <- NA
    return(aed)
}


#' @importFrom picante pd
#' @importFrom picante evol.distinct
.haed <- function(data){
    #Argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")

    #Setup
    ed <- evol.distinct(data$phy, "fair.proportion")$w
    pd <- pd(data$comm, data$phy)$PD
    aed <- .aed(data)

    #Internal assemblage calc.
    ..haed <- function(ed, pd.comm, aed.comm, assemblage){
        s.aed <- rep(aed.comm, assemblage) / pd.comm      
        haed <- -sum(s.aed * log(s.aed), na.rm=TRUE)
        return(haed)
    }

    #Calculate, clean, and return
    output <- numeric(nrow(data$comm))
    names(output) <- rownames(data$comm)
    for(i in seq(nrow(data$comm)))
        output[i] <- ..haed(ed, pd[i], aed[,i], data$comm[i,])
    return(output)
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

#' @importFrom picante evol.distinct
.scheiner <- function(data, q=0, abund = TRUE){
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm")) stop("'data' must be a comparative community ecology object")
    
    #Setup
    ed <- evol.distinct(data$phy, "fair.proportion")$w
    pd <- pd(data$comm, data$phy)$PD
    if(!abund)
        data$comm <- as.numeric(data$comm > 0)
    
    #Calculate scheiner; beware dividing by zero inadvertantly
    output <- numeric(nrow(data$comm))
    for(i in seq(nrow(data$comm))){
        if(q==1)
            output[i] <- exp(-1*sum(((data$comm[i,]*ed[i])/(sum(data$comm[i,])*pd[i])) * log((data$comm[i,]*ed[i])/(sum(data$comm[i,])*pd[i]))))
        else
            output[i] <- sum(((data$comm[i,]*ed[i])/(sum(data$comm[i,])*pd[i]))^q)^(1/(1-q))
    }
    return(output)
}
