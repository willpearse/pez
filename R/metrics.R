#Internal Colless function
#' @importFrom apTreeshape colless tipsubtree
.colless <- function(data, ...)
{
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    output <- numeric(nrow(data$comm))
    for(i in seq(nrow(data$comm)))
        output[i] <- colless(as.treeshape(drop_tip(data$phy, colnames(data$comm)[data$comm[i,]==0])))
    names(output) <- rownames(data$comm)
    return(output)
}

#' @importFrom picante pd evol.distinct
.hed <- function(data, ...){
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

.eed <- function(data, na.rm=TRUE, ...) {
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    output <- .hed(data) / log(apply(data$comm, 1, function(x) sum(x != 0)))
    names(output) <- rownames(data)
    return(output)
}

#' @importFrom picante psv
.psv <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    psv(x$comm, x$phy)[,1]
}

#' @importFrom picante psd
.psr <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    psd(x$comm, x$phy)[,4]
}

#' @importFrom picante mpd
#' @importFrom ape cophenetic.phylo
.mpd <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    mpd(x$comm, dist, abundance.weighted=abundance.weighted)
}

#' @importFrom picante pd
#' @importFrom stats lm
.pd <- function(x, include.root=TRUE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    pd <- pd(x$comm, x$phy, include.root)[,1]
    pd.ivs <- unname(resid(lm(coefs$pd ~ rowSums(x$comm))))
    return(cbind(pd, pd.ivs))
}

#' @importFrom picante mntd
#' @importFrom ape cophenetic.phylo
.mntd <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(mntd(x$comm, dist, abundance.weighted=abundance.weighted))
}

#' @importFrom ape gammaStat
#' @importFrom apTreeshape as.treeshape
.gamma <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    
    ..gamma <- function(pa.vec,tree,nams){
        if(sum(pa.vec)<3){
            return(NA)
        } else {
            if(length(setdiff(tree$tip.label, nams)) != 0){
                tree <- drop_tip(setdiff(tree$tip.label, ))
            }
            return(gammaStat(drop_tip(tree,nams[pa.vec==0])))
        }
    }
    
    tree.shape <- as.treeshape(x$phy)
    nams <- tree.shape$names
    return(apply(x$comm, 1, ..gamma, x$phy, nams))
}

#' @importFrom ape cophenetic.phylo
#' @importFrom vegan taxondive 
.taxon <- function(x, dist=NULL, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    output <- taxondive(x$comm, dist)
    output <- with(output, data.frame(Delta=D, DeltaStar=Dstar, LambdaPlus=Lambda, DeltaPlus=Dplus, S.DeltaPlus=SDplus))
    return(output)
}

#' @importFrom ape cophenetic.phylo
.eigen.sum <- function(x, dist=NULL, which.eigen=1, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    ..eigen.sum <- function(x, evc, vecnums) {
        if(sum(x>0)) {
            return(sum(apply(as.matrix(evc[x>0,vecnums]),2,var)))
        } else {
            return(NA)
        }
    }
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    
    eigen <- -0.5 * dist
    l <- matrix(1/nrow(eigen), nrow=nrow(eigen), ncol=ncol(eigen))
    eigen <- (diag(nrow(eigen)) - l) %*% eigen %*% (diag(nrow(eigen)) - l)
    eigen <- eigen(eigen, symmetric=TRUE)$vectors
    return(apply(x$comm, 1, ..eigen.sum, eigen, which.eigen))
}

#' @importFrom FD dbFD
#' @importFrom ape cophenetic.phylo
.dist.fd <- function(x, method="phy", abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")

    if(method == "phy")
        data <- cophenetic(x$phy)
    if(method == "traits")
        data <- x$data
    if(is.matrix(method) | is.data.frame(method))
        data <- method
    output <- capture.output(dbFD(data, x$comm, w.abun=abundance.weighted, messages=TRUE), file=NULL)
    coefs <- with(output, cbind(coefs, cbind(FRic, FEve, FDiv, FDis, RaoQ)))
    
    #Only bother getting CWMs if we have trait data
    if(method=="traits" | is.data.frame(method)){
        t <- output$dist.fd$CWM
        colnames(t) <- paste(colnames(t), "cmw", sep=".")
        coefs$dist.fd <- rbind(coef$dist.fd, t)
    }
}

#' @importFrom ape is.ultrametric as.phylo cophenetic.phylo
.sqrt.phy <- function(x){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(!is.ultrametric(x$phy))
        stop("Phylogeny is not ultrametric; cannot square-root (known 'bug', see help)")
    dist <- sqrt(cophenetic(x$phy))
    x$phy <- as.phylo(hclust(as.dist(dist)))
    return(x)
}

#' @importFrom ape write.tree
#' @importFrom ade4 newick2phylog
.phylo.entropy <- function(data, ...)
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
.aed <- function(data, ...){
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
.haed <- function(data, ...){
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

.iac <- function(data, na.rm=TRUE, ...) {
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

.pae <- function(data, na.rm=TRUE, ...) {
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
.scheiner <- function(data, q=0, abundance.weighted = TRUE, ...){
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm")) stop("'data' must be a comparative community ecology object")
    
    #Setup
    ed <- evol.distinct(data$phy, "fair.proportion")$w
    pd <- pd(data$comm, data$phy)$PD
    if(!abundance.weighted)
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

#' @importFrom picante pse
.pse <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    pse(x$comm, x$phy)[,1]
}

#' @importFrom picante raoD
.rao <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    raoD(x$comm, x$phy)$Dkk
}

#' @importFrom caper pgls comparative.data
.lambda <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    
    output <- numeric(length(sites(x)))
    for(i in seq(nrow(x$comm))){
        c.data <- x$comm[i,][x$comm[i,] > 0]
        if(length(unique(c.x)) > 1){
            c.data <- data.frame(comm=c.data, names=names(c.data))
            c.data <- comparative.data(phy=drop_tip(x$phy, setdiff(x$phy$tip.label, c.data$names)), data=c.data, names.col=names)
            model <- pgls(comm ~ 1, c.data, lambda="ML", ...)
            output[i] <- summary(model)$param.CI$lambda$opt
        } else output[i] <- NA
    }
    return(output)
}

#' @importFrom caper pgls comparative.data
.delta <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    
    output <- numeric(length(sites(x)))
    for(i in seq(nrow(x$comm))){
        c.data <- x$comm[i,][x$comm[i,] > 0]
        if(length(unique(c.x)) > 1){
            c.data <- data.frame(comm=c.data, names=names(c.data))
            c.data <- comparative.data(phy=drop_tip(x$phy, setdiff(x$phy$tip.label, c.data$names)), data=c.data, names.col=names)
            model <- pgls(comm ~ 1, c.data, delta="ML", ...)
            output[i] <- summary(model)$param.CI$delta$opt
        } else output[i] <- NA
    }
    return(output)
}

#' @importFrom caper pgls comparative.data
.kappa <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    
    output <- numeric(length(sites(x)))
    for(i in seq(nrow(x$comm))){
        c.data <- x$comm[i,][x$comm[i,] > 0]
        if(length(unique(c.x)) > 1){
            c.data <- data.frame(comm=c.data, names=names(c.data))
            c.data <- comparative.data(phy=drop_tip(x$phy, setdiff(x$phy$tip.label, c.data$names)), data=c.data, names.col=names)
            model <- pgls(comm ~ 1, c.data, kappa="ML", ...)
            output[i] <- summary(model)$param.CI$kappa$opt
        } else output[i] <- NA
    }
    return(output)
}

.eaed <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(.haed(x)/log(rowSums(x$comm)))
}


#' @importFrom picante unifrac
.unifrac <- function(x, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(unifrac(x$comm, x$phy))
}

#' @importFrom picante ses.mpd
.pcd <- function(x, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(pcd(x$comm, x$phy, rep=permute))
}

#' @importFrom picante comdist
.comdist <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    return(comdist(x$comm, dist, abundance.weighted=abundance.weighted))
}

#' @importFrom picante phylosor
.phylosor <- function(x, dist=NULL, abundance.weighted=FALSE, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    if(is.null(dist))
        dist <- cophenetic(x$phy)
    output <- phylosor(data$comm, data$phy)
    output <- as.dist(1 - as.matrix(output))
    return(output)
}


#Interal D calculation. *Heavily* based on D from caper.
# - placed in here for future compatibility with Dc
#' @importFrom caper contrCalc VCV.array
#' @importFrom mvtnorm rmvnorm
.d <- function(data, permute=1000, ...) {
  #Checking
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  if (!is.numeric(permute)) (stop("'", permute, "' is not numeric."))
  data$comm[data$comm > 1] <- 1
  # check tree branch lengths
  el    <- data$phy$edge.length
  elTip <- data$phy$edge[,2] <= length(data$phy$tip.label)
  
  if(any(el[elTip] == 0)) 
    stop('Phylogeny contains pairs of tips on zero branch lengths, cannot currently simulate')
  if(any(el[! elTip] == 0)) 
    stop('Phylogeny contains zero length internal branches. Use di2multi.')


  #Internal D calculation
  ..d <- function(ds, vcv, permute, phy){
      dsSort <- sort(ds)
      
      ## Random Association model
      ds.ran <- replicate(permute, sample(ds))
      
      ## Brownian Threshold model random data
      ds.phy <- rmvnorm(permute, sigma=unclass(vcv)) # class of 'VCV.array' throws the method dispatch
      ds.phy <- as.data.frame(t(ds.phy))
      
                                        # turn those into rank values, then pull that rank's observed value
      ds.phy <- apply(ds.phy, 2, rank, ties="random")
      ds.phy <- apply(ds.phy, 2, function(x) as.numeric(dsSort[x]))
      
      ## Get change along edges
      ## insert observed and set dimnames for contrCalc
      ds.ran <- cbind(Obs=ds, ds.ran)
      ds.phy <- cbind(Obs=ds, ds.phy)
      dimnames(ds.ran) <- dimnames(ds.phy) <- list(data$phy$tip.label, c('Obs', paste('V',1:permute, sep='')))
      
      ## now run that through the contrast engine 
      ds.ran.cc <- contrCalc(vals=ds.ran, phy=phy, ref.var='V1', picMethod='phylo.d', crunch.brlen=0)
      ds.phy.cc <- contrCalc(vals=ds.phy, phy=phy, ref.var='V1', picMethod='phylo.d', crunch.brlen=0)
      
      ## get sums of change and distributions
      ransocc <- colSums(ds.ran.cc$contrMat)
      physocc <- colSums(ds.phy.cc$contrMat)
      # double check the observed, but only to six decimal places or you can get floating point errors
      if(round(ransocc[1], digits=6) != round(physocc[1], digits=6)) stop('Problem with character change calculation in phylo.d')
      obssocc <- ransocc[1]
      ransocc <- ransocc[-1]
      physocc <- physocc[-1]
      
      soccratio <- (obssocc - mean(physocc)) / (mean(ransocc) - mean(physocc))
      soccpval1 <- sum(ransocc < obssocc) / permute
      soccpval0 <- sum(physocc > obssocc) / permute
      
      return(c(soccratio, soccpval1, soccpval0))
  }

  ## being careful with the edge order - pre-reorder the phylogeny
  phy <- reorder(data$phy, 'pruningwise')
  
  vcv <- VCV.array(data$phy)
  vals <- matrix(ncol=3, nrow=nrow(data$comm))
  rownames(vals) <- rownames(data$comm)
  colnames(vals) <- c("D", "P(D=1)", "P(D=0)")
  for(i in seq(nrow(data$comm)))
      vals[i,] <- ..d(data$comm[i,], vcv, permute, data$phy)
  
  return(vals)
}

#' @importFrom picante ses.mpd
.ses.mpd <- function(x, dist=NULL, null.model="taxa.labels", abundance.weighted=FALSE, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(ses.mpd(x$comm, dis=dist, null.model=null.model, abundance.weighted=abundance.weighted, runs=permute))
}

#' @importFrom picante ses.mntd
.ses.mntd <- function(x, dist=NULL, null.model="taxa.labels", abundance.weighted=FALSE, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(ses.mntd(x$comm, dis=dist, null.model=null.model, abundance.weighted=abundance.weighted, runs=permute))
}

#' @importFrom picante ses.mpd
.inmd <- function(x, dist=NULL, null.model="taxa.labels", abundance.weighted=FALSE, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(ses.mpd(x$comm, dis=1/dist, null.model=null.model, abundance.weighted=abundance.weighted, runs=permute))
}

#' @importFrom picante ses.mntd
.innd <- function(x, dist=NULL, null.model="taxa.labels", abundance.weighted=FALSE, permute=1000, ...){
    if(!inherits(x, "comparative.comm"))
        stop("'x' must be a comparative.comm object")
    return(ses.mntd(x$comm, dis=1/dist, null.model=null.model, abundance.weighted=abundance.weighted, runs=permute))
}
