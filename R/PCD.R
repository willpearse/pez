##' @importFrom ape compute.brlen
PCD <- function(samp, tree, PSVmncd=NULL, PSVpool=NULL, reps=10^4)
{
  SSii<-PSVmncd
  SCii<-PSVpool
  
  # Make samp matrix a pa matrix
  samp[samp>0]<-1
  
  # convert trees to VCV format
  if (is(tree)[1] == "phylo")
  {
    if (is.null(tree$edge.length)) {tree <- compute.brlen(tree, 1)}    #If phylo has no given branch lengths
    tree <- prune.sample(samp, tree)
    V <- vcv.phylo(tree, corr = TRUE)
    samp <- samp[, tree$tip.label]
  } else {
    V <- tree
    species<-colnames(samp)
    preval<-colSums(samp)/sum(samp)
    species<-species[preval>0]
    V<-V[species,species]
    samp<-samp[,colnames(V)]
  }
  
  if(any(rowSums(samp)==0))  #it is possible that some species are not included in the phylogeny and as such their removal creates communities of zero spp
  {
    warning("Some of the communities have zero species.")
  }
  samp.orig<-samp
  samp<-samp[rowSums(samp)>0,]
  m <- dim(samp)[1]
  n <- dim(samp)[2]
  
  
  if (!is.null(SSii) & length(SSii)!=max(rowSums(samp))) {
    stop("The length of PSVmncd does not match the community with the highest species richness")
  }
  
  # m=number of communities; n=number of species; nsr=maximum sr value across all communities
  
  nsr <-max(rowSums(samp))
  if(is.null(SSii) & is.null(SCii))   #If the user already has calculated the mean conditional PSV values for all levels of SR
  {                                   #and the PSV of the species pool
    SSii <- array(0,nsr)
    n1 <- 2
    for (n2 in 1:nsr)
    {
      temp <- array(0,reps)
      for (t in 1:reps)
      {
        rp <- sample(n)
        pick1 <- rp[1:n1]
        
        rp <- sample(n)
        pick2 <- rp[1:n2]
        
        C11 <- V[pick1,pick1]
        C22 <- V[pick2,pick2]
        C12 <- V[pick1,pick2]
        
        invC22 <- solve(C22)
        S11 <- C11 - C12%*%invC22%*%t(C12)
        SS11 <- (n1*sum(diag(S11))-sum(S11))/(n1*(n1-1))
        temp[t] <- SS11
      }
      SSii[n2] <- mean(temp)
    }
    SCii=1-(sum(V)-sum(diag(V)))/(n*(n-1))
  }
  
  # calculate PCD
  PCD <- array(NA,c(m,m))
  PCDc <- array(NA,c(m,m))
  PCDp <- array(NA,c(m,m))
  for (i in 1:(m-1))
  {
    for (j in (i+1):m)
    {
      pick1 <- (1:n)[samp[i,]==1]
      pick2 <- (1:n)[samp[j,]==1]
      
      n1 <- length(pick1)
      n2 <- length(pick2)
      
      C <- V[c(pick1, pick2),c(pick1, pick2)]
      
      C11 <- C[1:n1,1:n1]
      C22 <- C[(n1+1):(n1+n2),(n1+1):(n1+n2)]
      C12 <- C[1:n1,(n1+1):(n1+n2)]
      if(is.null(dim(C12)))
      {
        if(is.null(dim(C22))){C12<-as.matrix(C12)} else {C12<-t(as.matrix(C12))}
      }
      
      invC11 <- solve(C11)
      S22 <- C22 - t(C12)%*%invC11%*%C12
      
      invC22 <- solve(C22)
      S11 <- C11 - C12%*%invC22%*%t(C12)
      if(n1>1)
      {
        SC11 <- (n1*sum(diag(C11))-sum(C11))/(n1*(n1-1))
        SS11 <- (n1*sum(diag(S11))-sum(S11))/(n1*(n1-1))
      } else {          #Added to deal with communities of only one species
        SC11 <- 0
        SS11 <- S11
      }
      if(n2>1)
      {
        SC22 <- (n2*sum(diag(C22))-sum(C22))/(n2*(n2-1))
        SS22 <- (n2*sum(diag(S22))-sum(S22))/(n2*(n2-1))
      } else {          #Added to deal with communities of only one species
        SC22 <- 0
        SS22 <- S22
      }
      
      if((n1+n2)==2){ #both communities have only one species
        D=(n1*SS11 + n2*SS22) #we do not standardize by the unconditional PSVs, since the PSV of a 1 spp community is zero (or undefined)
      } else {
        D=(n1*SS11 + n2*SS22)/(n1*SC11 + n2*SC22)   #if one of the communities is of one species, then the unconditional PSV is 0
      }
      a <- length(unique(c(pick1, pick2)))
      b <- length(pick1)-a
      cc <- length(pick2)-a
      dsor <- 2*a/(2*a+b+cc) - 1
      
      pred.D <- (n1*SSii[n2]+n2*SSii[n1])/(n1*SCii+n2*SCii)
      pred.dsor <- 1 - 2*n1*n2/((n1+n2)*n)
      
      PCD[i,j] <- D/pred.D
      PCDc[i,j] <- dsor/pred.dsor
      PCDp[i,j] <- PCD[i,j]/PCDc[i,j]
    }
  }
  colnames(PCD)<-rownames(samp)
  rownames(PCD)<-rownames(samp)
  colnames(PCDc)<-rownames(samp)
  rownames(PCDc)<-rownames(samp)
  colnames(PCDp)<-rownames(samp)
  rownames(PCDp)<-rownames(samp)
  
  return(list(PCD=PCD,PCDc=PCDc,PCDp=PCDp,PSVmncd=SSii,PSVpool=SCii))
}
