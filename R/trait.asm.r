#TODO: NULL
#
#' Produces simulated communities based on species attributes
#' 
#' \code{trait.asm} calculates phylogenetic biodiversity metrics
#' 
#' @param a species attributes (e.g., traits like body size)
#' @param m number of communities to be simulated
#' @param meanSR target mean species richness across simulated communities
#' @param interval adjust to obtain \code{meanSR}
#' @param exponential use the exponential distribution when simulating communities
#' @param Pscale adjust this value when not using the exponential distribution in order to scale the species richnesses in the simulated communities
#' @details Simulates a set of communties based on the supplied attribute (trait) where larger values make it more likely for species to be in the communities.
#
#' @return \code{Y} presence/absence matrix
#' @return \code{P} probabilities
#' @return \code{a} the supplied trait
#' @return \code{exponential} if the exponential distribution was used
#' @return \code{meanSR} supplied \code{meanSR} value
#' @return \code{std} estimated sd
#' @author M.R. Helmus
#' @references Helmus M., Mercado-Silva N. & Vander Zanden M.J. (2013). Subsidies to predators, apparent competition and the phylogenetic structure of prey communities. Oecologia, 173, 997-1007.
#' @examples
#' \dontrun{
#'  data(laja)
#'  trait.asm(laja$fish.pref)
#' }
trait.asm<-function(a, m=1000,meanSR=NULL,interval=c(.001,10),exponential=TRUE,Pscale=1)
{
    a<-a[!is.na(a)]
	  n<-length(a)
    e<-as.matrix(a)
    
  if(!is.null(meanSR)){
    zz<-function(std,e,m,n,meanSR)
    {
      P <- exp(array(e,c(n,m))/std)/exp(max(e)/std)
      p <- array(P,c(n*m,1))
      Y <- rbinom(n*m, 1, p)
      Y <- t(array(Y, c(n,m)))
      abs(meanSR-mean(rowSums(Y)))
    }
    std<-unlist(optimize(zz,interval=interval,n=n,m=m,e=e,meanSR=meanSR))[1]    
  } else {std<-1}

  #exponential distribution
  if(exponential)
  {
    P <- (exp(array(e,c(n,m))/std)/exp(max(e)/std))
  } else {
    #do not use exponential distribution
    P <- Pscale*array(e,c(n,m))
  }
  p <- array(P,c(n*m,1))
	# convert distribution to presence/absence
	Y <- rbinom(n*m, 1, p)
	Y <- t(array(Y, c(n,m)))
  colnames(Y)<-names(a)
  return(list(Y=Y,P=P,a=a,exponential=exponential,meanSR=meanSR,std=std))
}