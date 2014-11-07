# TO_DO:
# clean up Rd code, write nice example
#' Simulate phylogenetic community structure across a landscape
#' 
#' \code{scape} simulates communities that are phylogenetically structured
#' 
#' @param tree \code{\link{phylo}} object
#' @param scape.size edge dimension of square landscape
#' @param g.center strength of phylogenetic signal in species range centers
#' @param g.range strength of phylogenetic signal in species range sizes
#' @param g.repulse strength of phylogenetic repulsion
#' @param wd.all niche width, larger values simulate broader range sizes
#' @param signal.center simulate with phylosignal in range centers
#' @param signal.range simulate with phylosignal in range size
#' @param same.range make all range sizes equal
#' @param repulse include phylogenetic repulsion in range centers
#' @param center.scale adjust strength of phylogenetic attraction in range centers independent of \code{g.center}
#' @param range.scale adjust strength of phylogenetic signal in range size independent of \code{g.range}
#' @param repulse.scale adjust strength of phylogenetic repulsion independent of \code{g.repulse}
#' @param site.stoch.scale adjust strength of random variation in species richness across sites
#' @param sd.center sd in \code{\link{rnorm}} for the range centers, increase to get more variation in center values across species
#' @param sd.range sd \code{\link{rnorm}} for the range sizes, increase to get more variation in range sizes across gradients
#' @param rho Grafen branch adjustment of phylogenetic tree see \code{\link{corGrafen}}
#' @param th probability threshold 10^-th above which species are considered present at a site 
#' @details Simulates a landscape with species (i.e., tree tips)
#' distributions dependent on a supplied phylogenetic tree.  The
#' amount and type of structure is determened by the signal parameters
#' \code{g.center}, \code{g.range} and \code{g.repulse}. Parameters
#' are based on an Ornstein-Uhlenbeck model of evolution with
#' stabilizing selection. Values of g=1 indicate no stabilizing
#' selection and correspond to the Brownian motion model of evolution;
#' 0<g<1 represents stabilizing selection; and g>1 corresponds to
#' disruptive selection where phylogenetic signal for the supplied
#' tree is amplified. See \code{\link{corBlomberg}}.  Communities are
#' simulated along two gradients where the positions along those
#' gradients, \code{g.center} and range sizes \code{g.range}, can
#' exhibit phylogenetic signal. Phylogenetic attraction is simulated
#' in the \code{g.center} paramter, while repulsion in
#' \code{g.repulse}. Both can be exhibited such that closly related
#' species are generally found with similar range centers
#' (phylogenetic attraction) but just not at the same site
#' (phylogenetic repulsion).  The function then returns probabilities
#' of of each species across sites and the presence and absences of
#' species based a supplied threshold, \code{th}, which can be
#' increased to obtain more species at sites and thus increase average
#' site species richness.
#
#' @return \code{Y} presence/absence matrix
#' @return \code{index} spatial coordinates for X and Y (stacked columns)
#' @return \code{X.joint} full probabilities of species at sites, used to construct Y 
#' @return \code{X1} probabilities of species along gradient 1
#' @return \code{X2} probabilities of species along gradient 2
#' @return \code{sppXs} full probabilities of each species as an array arranged in a \code{scape.size}-by-\code{scape.size} matrix  
#' @return \code{V.phylo} initial phylogenetic covariance matrix from tree
#' @return \code{V.phylo.rho} phylogenetic covariance matrix from tree scaled by Grafen if rho is provided
#' @return \code{V.center} scaled by \code{g.center} phylo covariance matrix used in the simulations
#' @return \code{V.range}  scaled by \code{g.range} phylo covariance matrix used in the simulations
#' @return \code{V.repulse} scaled by \code{g.repulse} phylo covariance matrix used in the simulations
#' @return \code{bspp1} species optima for gradient 1
#' @return \code{bspp2} species optima for gradient 2                               
#' @return \code{u} the env gradients values for the two gradients
#' @return \code{wd} the denominator for species ranges 
#' @author M.R. Helmus
#' @references Helmus M.R. & Ives A.R. (2012). Phylogenetic diversity area curves. Ecology, 93, S31-S43.
#' @seealso \code{\link{sim.phy}}
#' @importFrom ape vcv corBlomberg compute.brlen vcv.phylo
#' @export
#' @examples
#' require(ape)
#' tree<-stree(8,type="balanced")       #signal in centers
#' kk<-scape(tree, scape.size=100, g.center=100, g.range=1, g.repulse=1, wd.all=150,
#' signal.center=TRUE, signal.range=FALSE, same.range=FALSE, repulse=FALSE,center.scale = 1,
#' range.scale = 1, repulse.scale = 1, site.stoch.scale = 0, sd.center=3, sd.range=1,rho=NULL, th=20)
#'
#' #Make some plots
#' require(plotrix)
#' par(mfrow=c(1,Ntip(tree)),mar=c(.1,.1,.1,.1))
#' for(j in 1:Ntip(tree)){color2D.matplot(1 - kk$sppXs[,,j]/max(kk$sppXs[,,j]), xlab = "",
#' ylab = "",main = "",border=NA,do.hex=FALSE,axes=FALSE)}
#' 
#' par(mfrow=c(2,1))
#' matplot((kk$X1), type = "l", xlab="gradient",ylab = "probability",
#' main = "Gradient 1",col=rainbow(dim(kk$X1)[2]),lty=1)
#' matplot((kk$X2), type = "l", xlab="gradient",ylab = "probability",
#' main = "Gradient 2",col=rainbow(dim(kk$X2)[2]),lty=1)
#' 
#' plot(x=1:dim(kk$Y)[1],y = rowSums(kk$Y), main = "SR",type = "l")
#' cor(kk$X1)
scape<-function(tree, scape.size=10, g.center=1, g.range=1, g.repulse=1, wd.all=150, signal.center=TRUE, signal.range=TRUE, same.range=TRUE,repulse=TRUE,center.scale = 1, range.scale = 1, repulse.scale = 1, site.stoch.scale = .5, sd.center=1, sd.range=1,rho=NULL, th=8)
{
  #deal with the tree
    if (is(tree)[1] == "phylo")
    {
        if (is.null(tree$edge.length))
        {
            tree <- compute.brlen(tree, 1)    #Note this assigns arbitrary branch-lengths
        }
        V <- vcv.phylo(tree, corr = TRUE)     #Note this will mess with trees that are not ultrametric, such as those with argitrary branch-lengths
    } else {
        V <- tree
    }
    Vinit<-V       
  #initialize
    nspp <- dim(V)[1]
    bspp2 <- NULL
    Vcomp <- NULL

        Xscale <- 1          #scale the strength of the probability matrix X
        Mscale <- site.stoch.scale        #scale stochasticity in niche distributions 
        Vscale1 <- center.scale         #scale the strength of the optimal values on axis one
        Vscale2 <- center.scale         #scale the strength of the optimal values on axis two
    
    # Grafen's rho adjust strength of phylogenetic signal overall.
    if(!is.null(rho)){
     V <- 1-(1-Vinit)^rho
     V <- V/max(V)
    }
 
    #########################################################################################################################
    #SIMULATION
    nsites<-scape.size                                                #number of sites for the square landscape
    mx <- t(as.matrix((-(nsites)/2):(nsites/2)))  #env gradient
    
    m <- length(mx)                                                   #new number of sites (equal to nsites + 1)
      
    ############
    #ATTRACTION (RANGE CENTERS,NICHE OPTIMA) AND RANGE WIDTHS
    if(signal.center){
      g<-abs(g.center)
      V.a<-vcv(corBlomberg(g, tree),corr=T)           #adjust phylogenetic signal
      iD <- t(chol(V.a))
    } else {
      V.a<-V
      iD <- diag(nspp)
    }
    if(signal.range){
      g<-abs(g.range)
      V.w<-vcv(corBlomberg(g, tree),corr=T)           #adjust phylogenetic signal
      iD.w <- t(chol(V.w))
    } else {
      V.w<-V
      iD.w <- diag(nspp)
    }
  ##environmental/geographical gradient 1   

    e <- iD %*% rnorm(nspp,sd=sd.center)                                                               #assign optimal values 
    #e <- iD %*% runif(nspp)                                                               #assign optimal values 
    #e <- iD %*% rep(1,nspp)                                                               #assign optimal values 
    e <- Vscale1 * (e - mean(e))/apply(e, 2, sd)                                           #z-scores and optional scaling of the optimal values
    bspp1 <- e
    if(!same.range){
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          wd <- range.scale*iD.w %*% rnorm(nspp,sd=sd.range)
          wd<-wd+(abs(min(wd)))
          wd<-wd/max(wd)

          dif<-sort(wd)[-1]-sort(wd)[-length(wd)]
          rat<-mean(dif/sort(wd)[-1])
          wd[wd==0]<-sort(wd)[2]-sort(wd)[2]*rat                                #Assign the zero with the mean ratio of nearest neighbor distances over the larger item
          wd<-wd.all*wd
          X <- exp(-((spmx - mxsp)^2)/t(matrix(wd,nspp,m))) #Niche distributions
        } else {
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          X <- exp(-((spmx - mxsp)^2)/wd.all) #Niche distributions
    }       
    X <- Xscale * X                                                                       #Scales this initial species x site probability matrix 
    #Xsmooth <- X                                                                         #Distributions without random variation
    X1 <- diag(1 - Mscale * runif(m)) %*% X                                               #Scale and include random variation into the niche distributions
    #}
    
    ##environmental/geographical gradient 2
    #if (envirogradflag2 == 1) {
        e <- iD %*% rnorm(nspp,sd=sd.center)
        e <- Vscale2 * (e - mean(e))/apply(e, 2, sd)
        bspp2 <- e

     if(!same.range){
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          wd <- range.scale*iD.w %*% rnorm(nspp,sd=sd.range)
          wd<-wd+(abs(min(wd)))
          wd<-wd/max(wd)
              #Assign the zero to the nonzero minimum
          dif<-sort(wd)[-1]-sort(wd)[-length(wd)]
          rat<-mean(dif/sort(wd)[-1])
          wd[wd==0]<-sort(wd)[2]-sort(wd)[2]*rat                                #Assign the zero with the mean rato of nearist neighbor distances over the larger item
          wd<-wd.all*wd
          X <- exp(-((spmx - mxsp)^2)/t(matrix(wd,nspp,m))) #Niche distributions     
     } else {
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          X <- exp(-((spmx - mxsp)^2)/wd.all) #Niche distributions
     }
     X <- Xscale * X
     X2 <- diag(1 - Mscale * runif(m)) %*% X

     ##################################        
     #REPULSION
     X.repulse <- NULL
      if (repulse) {
        compscale <- repulse.scale
        b0scale <- 0
        g<-abs(g.repulse)
        V.r<-vcv(corBlomberg(g, tree),corr=T)          #adjust phylogenetic signal 
      #calculate the repulsion matrix
        Vcomp <- solve(V.r, diag(nspp))
        Vcomp <- Vcomp/max(Vcomp)
        Vcomp <- compscale * Vcomp
        iDcomp <- t(chol(Vcomp))
        colnames(Vcomp) <- rownames(Vcomp)
        bcomp <- NULL
        for (i in 1:m) {
          bcomp <- cbind(bcomp, iDcomp %*% rnorm(nspp))
        }
        bcomp0 <- 0
        Xcomp <- exp(bcomp0 + bcomp)/(1 + exp(bcomp0 + bcomp))
        #X <- X * t(Xcomp)
        X1<-X1 * t(Xcomp)
        X2<-X2 * t(Xcomp)
        X.repulse<-t(Xcomp)
      }
      
##################################
#JOINT PROBABILITY MATRIX
  X.<-NULL
  spp.Xs<-array(NA,dim=c(m,m,nspp))
  for(i in 1:nspp){
    sppX<-matrix((X1[,i]) %*% t(X2[,i]))
    spp.Xs[,,i]<-sppX
    X.<-cbind(X.,matrix(sppX))
   }
   colnames(X.)<-colnames(X2)
  ######################
  #PA matrix
  m.<- dim(X.)[1]
  Y <- matrix(0, ncol = nspp, nrow = m.)
  Y[10^-th < X.] <- 1     #could also use a hard threshold
  colnames(Y) <- colnames(X.)
  index<-NULL
  
  index<-cbind(matrix(sapply(1:m,rep,times=m)),matrix(rep(1:m,times=m)))
  colnames(index)<-c("X1","X2")
  
  ########### OUTPUT  
    return(list(Y = Y, index=index, 
                X.joint=X., X1=X1, X2=X2, sppXs=spp.Xs, 
                V.phylo=Vinit, V.phylo.rho = V, V.center = V.a, V.range = V.w, V.repulse = Vcomp, 
                bspp1 = bspp1, bspp2 = bspp2, u = mx, wd=wd.all))
  
}  #function end


