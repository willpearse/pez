#To-DO:
# optimise (? and if can be bothered)
# mis-match between proportions and partial.tree$edge.length - not introduced by the changes below
require(ade4)

phylo.entropy <- function(data)
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

