#TODO:
# Are Delta et al., when calculated here, based on abundances? Worth checking...
# There is an error message coming from phylo.entropy that is described in its script
# Could encapsulate phylogenetic signal measures by writing a wrapper and using expression to wrap the correct thing in (if you cared enough...)
# optimise (? and if can be bothered)
# mis-match between proportions and partial.tree$edge.length - not introduced by the changes below

#' Calculate evenness metrics across communities
#' 
#' \code{evenness} calculates evenness metrics in comparative.comm communities
#' 
#' @param data a comparative community ecology object
#' @param metric specify (a) particular metric(s) to calculate (all, rao, delta, entropy, cadotte, pst), or the default 'all'
#' @details This calculates Rao's quadratic entropy, taxonomic diversity, phylogenetic entropy, Cadotte's abundance measures, and Pst, all defined as evenness metrics in Pearse et al.
#' @return cc.evenness object (a named list with the output from each metric)
#' @author Matt Helmus, Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample)
#' evenness(data)
#' evenness(data, "rao")
#' }
#' @import picante ecoPD spacodiR ade4
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
    .abund <- data$comm[rowSums(data$comm>0)>1, ]
    if(length(setdiff(data$phy$tip.label, colnames(.abund))))
      tree <- drop.tip(data$phy, setdiff(data$phy$tip.label, colnames(.abund))) else tree <- data$phy
    temp <- phylo4d(tree, t(.abund))
    temp <- data.frame(PAE=pae(temp), IAC=iac(temp), Haed=haed(temp), Eaed=eaed(temp))
    output$cadotte <- temp[match(rownames(data$comm), rownames(temp)),]
  }
  
  if(metric == "pst" | metric == "all"){
    .abund <- data$comm[rowSums(data$comm>0)>1, ]
    if(length(setdiff(data$phy$tip.label, colnames(.abund))))
      tree <- drop.tip(data$phy, setdiff(data$phy$tip.label, colnames(.abund))) else tree <- data$phy
    temp <- phylo4d(tree, t(.abund))
    temp <- simpson(temp, "phylogenetic")
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
  class(output) <- c("cc.evenness", class(output))
  return(output)
}

#' Print an evenness object (...by summarising it...)
#' @method print cc.evenness
#' @S3method print cc.evenness
#' @export
print.cc.evenness <- function(x, ...){
  summary(x)
}

#' Summarise an evenness object
#' @method summary cc.evenness
#' @S3method summary cc.evenness
#' @export
summary.cc.evenness <- function(x, ...){
  cat("\nEvenness metrics in this object:\n")
  if(!is.null(x$rao)){
    cat("\tRao's quadratic entropy (rao)\n")
  }
  if(!is.null(x$taxon)){
    cat("\tTaxonomic diversity (taxon)\n")
  }
  if(!is.null(x$entropy)){
    cat("\tPhylogenetic entropy (entropy)\n")
  }
  if(!is.null(x$cadotte)){
    cat("\tCadotte's abundance evenness measures (cadotte)\n")
  }
  if(!is.null(x$pst)){
    cat("\tPst - Simpson's diversity (pst)\n")
  }
  if(!is.null(x$lambda)){
    cat("\tLambda transformation (lambda)\n")
  }
  if(!is.null(x$delta)){
    cat("\t Delta transformation (delta)\n")
  }
  if(!is.null(x$kappa)){
    cat("\tKappa transformation (kappa)\n")
  }
  cat("Use something like 'output$rao' to work with each measure\n")
}

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