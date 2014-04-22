#TODO:
# Unit tests
# Documentation
#' Calculate dispersion metrics across communities
#' 
#' \code{dispersion} calculates dispersion metrics in comparative.comm communities
#' 
#' @param data a comparative community ecology object
#' @param permute the number of null permutations to perform
#' @param metric specify (a) particular metric(s) to calculate (sesmpd, sesmntd, sespd, innd, d), or the default 'all'
#' @details This calculates $SES_{MPD}$, $SES_{MNTD}$, $SES_{PD}$, INND, and D. All these are defined as dispersion metrics in Pearse et al., 
#' Dc=0 is the Brownian expectation, Dc=1 is the random expectation.
#' @return cc.dispersion object (a named list with the output from each metric)
#' @author Matt Helmus, Will Pearse
#' @examples \dontrun{
#' data(phylocom)
#' data <- comparative.comm(phylocom$phy, phylocom$sample)
#' dispersion(data)
#' dispersion(data, 100, "sesmpd")
#' }
#' @import caper picante
#' @export

pblm.traits<-function(assocs,tree1=NULL,tree2=NULL,covars1=NULL,covars2=NULL,bootstrap=FALSE,nreps=10,maxit=10000,pstart=c(.5,.5)){

  # Make a vector of associations
  A<-as.matrix(as.vector(as.matrix(assocs)))
  data.vecs<-A

  #numbers of species and interactions
  nassocs<-length(A)
  nspp1<-dim(assocs)[1]
	nspp2<-dim(assocs)[2]
  sppnames1<-rownames(assocs)
  sppnames2<-colnames(assocs)
  #make names of species pairs
  pairnames=NULL  # make a vector of pairwise comparison names
  for (o in 1:(nspp2))
  {
    for (u in 1:nspp1)
    {
      pairnames<-c(pairnames,paste(sppnames2[o],sppnames1[u],sep="-"))
    }
  }

  #Clean Covariates
  #If the covariate applies to both sets, then it should be in the matrix of the longer set
  covnames<-NULL
  C1covs<-NULL
  if(is.null(covars1))
  {
    C1<-NULL
  } else {
    if(is.null(dim(covars1)))
    {
      C1<-matrix(covars1,nspp1,nspp2,byrow=FALSE)
      if(is.factor(covars1))
      {
        C1<-as.matrix(as.vector(C1))
        C1covs<-cbind(C1covs,C1)
        C1<-as.matrix(model.matrix(~as.factor(C1)-1)[,-1])
        colnames(C1)<-paste(rep("covar1",length(levels(covars1))-1),levels(covars1)[-1],sep="-")
      } else {
        C1<-as.matrix(as.vector(C1))
        C1covs<-cbind(C1covs,C1)
      }
      covnames<-c(covnames,"covar1")
    } else {
      C1<-NULL
      for(i in 1:dim(covars1)[2])
      {
        C1hold<-matrix(covars1[,i],nspp1,nspp2,byrow=FALSE)
        if(is.factor(covars1[,i]))
        {
          C1hold<-as.matrix(as.vector(C1hold))
          C1covs<-cbind(C1covs,C1hold)
          C1hold<-as.matrix(model.matrix(~as.factor(C1hold)-1)[,-1])
          colnames(C1hold)<-paste(rep(colnames(covars1)[i],length(levels(covars1[,i]))-1),levels(covars1[,i])[-1],sep="-")
          C1<-cbind(C1,C1hold)
        } else {
          C1hold<-as.matrix(as.vector(C1hold))
          C1covs<-cbind(C1covs,C1hold)
          colnames(C1hold)<-colnames(covars1)[i]
          C1<-cbind(C1,C1hold)
        }
        covnames<-c(covnames,colnames(covars1)[i])
      }
    }

  data.vecs<-cbind(data.vecs,C1covs)
  }


  C2covs<-NULL
  if(is.null(covars2))
  {
    C2<-NULL
  } else {
    if(is.null(dim(covars2)))
    {
      C2<-matrix(covars2,nspp1,nspp2,byrow=TRUE)
      if(is.factor(covars2))
      {
        C2<-as.matrix(as.vector(C2))
        C2covs<-cbind(C2covs,C2)
        C2<-as.matrix(model.matrix(~as.factor(C2)-1)[,-1])
        colnames(C2)<-paste(rep("covar2",length(levels(covars2))-1),levels(covars2)[-1],sep="-")
      } else {
        C2<-as.matrix(as.vector(C2))
        C2covs<-cbind(C2covs,C2)
      }
      covnames<-c(covnames,"covar2")
    } else {
      C2<-NULL
      for(i in 1:dim(covars2)[2])
      {
        C2hold<-matrix(covars2[,i],nspp1,nspp2,byrow=TRUE)
        if(is.factor(covars2[,i]))
        {
          C2hold<-as.matrix(as.vector(C2hold))
          C2covs<-cbind(C2covs,C2hold)
          C2hold<-as.matrix(model.matrix(~as.factor(C2hold)-1)[,-1])
          colnames(C2hold)<-paste(rep(colnames(covars2)[i],length(levels(covars2[,i]))-1),levels(covars2[,i])[-1],sep="-")
          C2<-cbind(C2,C2hold)
        } else {
          C2hold<-as.matrix(as.vector(C2hold))
          C2covs<-cbind(C2covs,C2hold)
          colnames(C2hold)<-colnames(covars2)[i]
          C2<-cbind(C2,C2hold)
        }
        covnames<-c(covnames,colnames(covars2)[i])
      }
    }

  data.vecs<-cbind(data.vecs,C2covs)
  }


# Make U, the combined matrix of covariates
  U<-NULL
  if(is.null(C1) & is.null(C2))
  {
    U<-rep(1,length(A))
  } else {

    if(is.null(C1))
    {
      U<-rep(1,length(A))
    } else {
      U<-cbind(rep(1,length(A)),C1)
    }

    if(is.null(C2))
    {
      U<-U
    } else {
      U<-cbind(U,C2)
    }
  }

  # Begin to organize output
  if(is.null(dim(U)))
  {
    data.vecs<-data.frame(A)
    colnames(data.vecs)<-"associations"
  } else {
    colnames(data.vecs)<-c("associations", covnames)
  }
  rownames(data.vecs)<-pairnames

  ######
  # Calculate Star Regression Coefficients
  #calculate for the star (assuming no phylogenetic correlation)
	astar<-solve((t(U)%*%U),(t(U)%*%A))
	MSETotal<-cov(A)
	s2aStar<-as.vector(MSETotal)*qr.solve((t(U)%*%U))
	sdaStar<-t(diag(s2aStar)^(.5))
	approxCFstar<-rbind(t(astar)-1.96%*%sdaStar, t(astar), t(astar)+1.96%*%sdaStar)
  Pstar<-U%*%astar
  Estar<-A-Pstar
  MSEStar<-cov(matrix(Estar))

  #######
  if(is.null(tree1) | is.null(tree2))
  {
    coefs<-approxCFstar
    rownames(coefs)<-c("lower CI 95%","estimate","upper CI 95%")
    colnames(coefs)<-paste("star",c("intercept",colnames(U)[-1]),sep="-")
    MSEs<-cbind(data.frame(MSETotal),data.frame(MSEStar))
    Pstar<-data.frame(Pstar)
    colnames(Pstar)<-"star"
    Estar<-data.frame(Estar)
    colnames(Estar)<-"star"
    output<-list(MSE=MSEs,signal.strength=NULL,coefficients=data.frame(t(coefs)),CI.boot=NULL,variates=data.frame(data.vecs),residuals=Estar,predicted=Pstar,bootvalues=NULL,Vfull=NULL)
    class(output)<-"pblm"
    return(output)

  } else {

    #tree1 is the phylogeny for the rows
    #tree2 is the phylogeny for the columns

    #Clean Trees
    if(is(tree1)[1]=="phylo")
    {
      if(is.null(tree1$edge.length)){tree1<-compute.brlen(tree1, 1)}  #If phylo has no given branch lengths
      V1<-vcv.phylo(tree1,corr=TRUE)
      V1<-V1[rownames(assocs),rownames(assocs)]
    } else {
      V1<-tree1[rownames(assocs),rownames(assocs)]
    }

    if(is(tree2)[1]=="phylo")
    {
    if(is.null(tree2$edge.length)){tree2<-compute.brlen(tree2, 1)}  #If phylo has no given branch lengths
      V2<-vcv.phylo(tree2,corr=TRUE)
      V2<-V2[colnames(assocs),colnames(assocs)]
    } else {
      V2<-tree2[colnames(assocs),colnames(assocs)]
    }

    #Calculate Regression Coefficents for the base (assuming strict brownian motion evolution, ds=1)
    V1<-as.matrix(V1)
    V2<-as.matrix(V2)
  	V1.orig<-V1
  	V2.orig<-V2

    V1<-V1/det(V1)^(1/nspp1)   # scale covariance matrices (this reduces numerical problems caused by
  	V2<-V2/det(V2)^(1/nspp2)   # determinants going to infinity or zero)
  	V<-kronecker(V2,V1)
    invV<-qr.solve(V)

    abase<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
    MSEBase<-(t(A-U%*%abase)%*%invV%*%(A-U%*%abase))/(nassocs-1)
    s2abase<-as.vector(MSEBase)*qr.solve(t(U)%*%invV%*%U)
  	sdabase<-t(diag(s2abase)^(.5))
    approxCFbase<-rbind(t(abase)-1.96%*%sdabase, t(abase), t(abase)+1.96%*%sdabase)
    Pbase<-t(t(U%*%abase)%*%invV)
    Ebase<-A-Pbase

    ###################
    # Full EGLS estimates of phylogenetic signal
    ##################
  	initV1<-V1
  	initV2<-V2
    init.invV<-invV
    initV<-V
  	# tau = tau_i + tau_j where tau_i equals the node to tip distance
  	tau1<-matrix(diag(initV1),nspp1,nspp1) + matrix(diag(initV1),nspp1,nspp1)-2*initV1
  	tau2<-matrix(diag(initV2),nspp2,nspp2) + matrix(diag(initV2),nspp2,nspp2)-2*initV2

    # The workhorse function to estimate ds
    pegls<-function(parameters)
    {
      d1<-abs(parameters[1])
      d2<-abs(parameters[2])

      V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
      V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)

      V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
      V2<-V2/det(V2)^(1/nspp2)
      V<-kronecker(V2,V1)
      invV<-qr.solve(V)

      a<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
      E<-(A-U%*%a)
      #MSE
      t(E)%*%invV%*%E/(nassocs-1)
    }
    # estimate d1 and d2 via Nelder-Mead method same as fminsearch in Matlab, by minimizing MSE
    est<-optim(pstart,pegls,control=list(maxit=maxit))
    MSEFull<-est$value
  	d1<-abs(est$par[1])
  	d2<-abs(est$par[2])

    # Calculate EGLS coef w estimated ds
  	V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
    V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
    V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
    V2<-V2/det(V2)^(1/nspp2)
    V<-kronecker(V2,V1)
    invV<-qr.solve(V)
    aFull<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)
    s2aFull<-as.vector(MSEFull)*qr.solve(t(U)%*%invV%*%U)
  	sdaFull<-t(diag(s2aFull)^(.5))
  	approxCFfull<-rbind(t(aFull)-1.96%*%sdaFull, t(aFull), t(aFull)+1.96%*%sdaFull)
    Pfull<-t(t(U%*%aFull)%*%invV)
    Efull<-A-Pfull

    ########################################

    #organize output
    coefs<-cbind(approxCFfull,approxCFstar,approxCFbase)
    rownames(coefs)<-c("approx lower CI 95%","estimate","approx upper CI 95%")
    colnames(coefs)<-c(paste("full",c("intercept",colnames(U)[-1]),sep="-"),paste("star",c("intercept",colnames(U)[-1]),sep="-"),paste("base",c("intercept",colnames(U)[-1]),sep="-"))
    coefs<-t(coefs)
    CI.boot<-NULL
    MSEs<-cbind(data.frame(MSETotal),data.frame(MSEFull), data.frame(MSEStar), data.frame(MSEBase))
    residuals<-cbind(data.frame(Efull),data.frame(Estar),data.frame(Ebase))
    predicted<-cbind(data.frame(Pfull),data.frame(Pstar),data.frame(Pbase))
    rownames(residuals)<-pairnames
    rownames(predicted)<-pairnames
    colnames(predicted)<-c("full","star","base")
    colnames(residuals)<-c("full","star","base")
    phylocovs=list(V1=V1,V2=V2,V1.orig=V1.orig,V2.orig=V2.orig,init.invV=init.invV,initV=initV,V=V)

    ################
    #bootstrap CIs
    if(bootstrap)
    {
      Vtrue<-V
		  Atrue<-A
		  atrue<-aFull
		  dtrue<-c(d1,d2)
		  ehold<-eigen(Vtrue,symmetric=TRUE)
      L<-ehold$vectors[,nassocs:1]    #A or L
      G<-sort(ehold$values)      #D
      iG<-diag(G^-.5)    #iD

      # Construct Y = TT*A so that
			# E{(Y-b)*(Y-b)'} = E{(TT*A-b)*(TT*A-b)'}
			#				  = T*V*T'
			#				  = I

      TT<-iG%*%t(L)
		  Y<-TT%*%Atrue
		  Z<-TT%*%U

		  res<-(Y-Z%*%atrue)	# residuals in orthogonalized space
		  invT<-qr.solve(TT)

      bootlist=NULL
      for (i in 1:nreps)
      {
        randindex<-sample(1:nassocs,replace=TRUE)	# vector of random indices
        #randindex=1:nassocs					# retain order
        YY<-Z%*%atrue+res[randindex]	# create new values of Y with random residuals
        A<-invT%*%YY	# back-transformed data
        pstart<-dtrue+c(0,.1)
        estRand<-optim(pstart,pegls,control=list(maxit=maxit))
        MSEFullrand<-estRand$value
        d1rand<-abs(estRand$par[1])
  	    d2rand<-abs(estRand$par[2])

        # Calculate EGLS coef w estimated ds
        V1<-(d1rand^tau1)*(1-d1rand^(2*initV1))/(1-d1rand^2)
        V2<-(d2rand^tau2)*(1-d2rand^(2*initV2))/(1-d2rand^2)
        V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
        V2<-V2/det(V2)^(1/nspp2)
        V<-kronecker(V2,V1)
        invV<-qr.solve(V)
        arand<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))   #NOTE: Ives in his Matlab code uses a Left matrix division symbol (\)

        bootlist<-rbind(bootlist,c(d1rand, d2rand, t(arand)))
      }
      nr<-dim(bootlist)[1]
      nc<-dim(bootlist)[2]

      #Calculate bootstrapped CIs
      alpha<-0.05  # alpha is always 0.05, but could change here
      conf<-NULL
      for(j in 1:nc)
      {
        bootconf<-quantile(bootlist[,j],probs = c(alpha/2, 1-alpha/2))
        conf<-rbind(conf,c(bootconf[1],bootconf[2]))
      }
      signal.strength<-data.frame(cbind(conf[1:2,1],dtrue,conf[1:2,2]))
      rownames(signal.strength)<-c("d1","d2")
      colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")

      #organize output
      CI.boot<-conf
      rownames(CI.boot)<-c("d1","d2","intercept",colnames(U)[-1])
      colnames(CI.boot)<-c("booted lower CI 95%","booted upper CI 95%")
      colnames(bootlist)<-c("d1","d2","intercept",colnames(U)[-1])
      output<-list(MSE=MSEs,signal.strength=signal.strength,assocs=assocs,covars1=covars1,covars2=covars2,
                    coefficients=data.frame(coefs),CI.boot=CI.boot,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=bootlist,phylocovs=phylocovs)
      class(output)<-"pblm"
      return(output)

    } else {
    ########
    # If bootstrapping not performed

    conf<-matrix(NA,2,2)
    signal.strength<-data.frame(cbind(conf[1,],c(d1,d2),conf[2,]))
    rownames(signal.strength)<-c("d1","d2")
    colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")
    output<-list(MSE=MSEs,signal.strength=signal.strength,assocs=assocs,covars1=covars1,covars2=covars2,coefficients=data.frame(coefs),CI.boot=NULL,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=NULL,phylocovs=phylocovs)
    class(output)<-"pblm"
    return(output)
    }
  }
}


###################################################################################################################################################
####### Takes a PBLM object and makes predictions for species not in phylogeny

pblmpredict.traits.novel<-function(x,tree1.w.novel=NULL,tree2.w.novel=NULL, covars1.novel=NULL,covars2.novel=NULL,calcboot=FALSE)
{
  if (!identical(class(x),"pblm")) stop("x must be of class pblm")
  if(is.null(x$phylocovs$V1) | is.null(x$phylocovs$V2)) stop("a pblm fit with phylogenies must be supplied")

  sppnames1.orig<-rownames(x$phylocovs$V1)
  sppnames2.orig<-rownames(x$phylocovs$V2)
  nspp1<-length(sppnames1.orig)
  nspp2<-length(sppnames2.orig)
  initV1.<-x$phylocovs$V1.orig
  initV2.<-x$phylocovs$V2.orig
  assocs.<-x$assocs
  covars1.orig<-as.matrix(x$covars1)
  covars2.orig<-as.matrix(x$covars2)
        	#tau1.init<-matrix(diag(initV1.),nspp1,nspp1) + matrix(diag(initV1.),nspp1,nspp1)-2*initV1.
 	        #tau2.init<-matrix(diag(initV2.),nspp2,nspp2) + matrix(diag(initV2.),nspp2,nspp2)-2*initV2.
  mu<-x$coefficients[grep("full",rownames(x$coefficients)),2]

  novel.spp1<-setdiff(tree1.w.novel$tip,sppnames1.orig)
  novel.spp2<-setdiff(tree2.w.novel$tip,sppnames2.orig)
    n.novels1<-length(novel.spp1)
    n.novels2<-length(novel.spp2)
    if(n.novels2>0 && n.novels1>0){stop("Sorry! Function does not yet work when species are missing from both phylogenies. Instead run one at a time.")}

      pairnames.orig <- rownames(x$variates)

  if(n.novels1>0)
  {
    assocs<-matrix(99,n.novels1,nspp2)
    rownames(assocs)<-novel.spp1
    colnames(assocs)<-sppnames2.orig
  } else {
    assocs<-matrix(99,n.novels2,nspp1)
    rownames(assocs)<-novel.spp2
    colnames(assocs)<-sppnames1.orig
    assocs<-t(assocs)
  }

    A<-as.matrix(as.vector(as.matrix(assocs)))

  sppnames1<-rownames(assocs)
  sppnames2<-colnames(assocs)
  #make names of species pairs

  if(n.novels1>0){
  pairnames=NULL  # make a vector of pairwise comparison names
  for (o in 1:(nspp2))
  {
    for (u in 1:n.novels1)
    {
      pairnames<-c(pairnames,paste(sppnames2[o],sppnames1[u],sep="-"))
    }
  }
  } else {

  pairnames=NULL  # make a vector of pairwise comparison names
  for (o in 1:(nspp1))
  {
    for (u in 1:n.novels2)
    {
      pairnames<-c(pairnames,paste(sppnames2[u],sppnames1[o],sep="-"))
    }
  }
  }
#######################################################################################################################################
  #Edit and Clean the novel phylogenies
  if(is.null(tree1.w.novel))          #tree1 is the phylogeny for the rows
  {
    V1.novel<-initV1.
  } else {
    #Clean Tree
    if(is(tree1.w.novel)[1]=="phylo")
    {
      if(is.null(tree1.w.novel$edge.length)){tree1.w.novel<-compute.brlen(tree1.w.novel, 1)}  #If phylo has no given branch lengths
      V1.novel<-vcv.phylo(tree1.w.novel,corr=TRUE)
      #V1.novel<-V1.novel[rownames(assocs),rownames(assocs)]
    } else {
      V1.novel<-tree1.w.novel#[rownames(assocs),rownames(assocs)]
    }
   }

  if(is.null(tree2.w.novel))         #tree2 is the phylogeny for the columns
  {
    V2.novel<-initV2.
  } else {
    #Clean Tree
    if(is(tree2.w.novel)[1]=="phylo")
    {
      if(is.null(tree2.w.novel$edge.length)){tree2.w.novel<-compute.brlen(tree2.w.novel, 1)}  #If phylo has no given branch lengths
      V2.novel<-vcv.phylo(tree2.w.novel,corr=TRUE)
      #V1.novel<-V1.novel[rownames(assocs),rownames(assocs)]
    } else {
      V2.novel<-tree2.w.novel#[rownames(assocs),rownames(assocs)]
    }
   }

   #scale V with novels with the calculated d values
    d1<-abs(x$sig[1,2])
    d2<-abs(x$sig[2,2])
    nspp1<-dim(V1.novel)[1]
    nspp2<-dim(V2.novel)[1]

    V1.novel<-V1.novel/det(V1.novel)^(1/nspp1)   # scale covariance matrices (this reduces numerical problems caused by
  	V2.novel<-V2.novel/det(V2.novel)^(1/nspp2)   # determinants going to infinity or zero)

  	# tau = tau_i + tau_j where tau_i equals the node to tip distance
  	tau1<-matrix(diag(V1.novel),nspp1,nspp1) + matrix(diag(V1.novel),nspp1,nspp1)-2*V1.novel
  	tau2<-matrix(diag(V2.novel),nspp2,nspp2) + matrix(diag(V2.novel),nspp2,nspp2)-2*V2.novel

    V1<-(d1^tau1)*(1-d1^(2*V1.novel))/(1-d1^2)
    V2<-(d2^tau2)*(1-d2^(2*V2.novel))/(1-d2^2)
    V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
    V2<-V2/det(V2)^(1/nspp2)
    V<-kronecker(V2,V1,make.dimnames=T)            #the values in this V will be different than the V in x since the novel species is included
    nams<-rownames(V)
    nams<-sub(":","-",nams)
    rownames(V)<-nams
    colnames(V)<-nams
    novel.assocs<-setdiff(rownames(V),pairnames.orig)
    orig.assocs<-pairnames.orig
    Vyy<-V[orig.assocs,orig.assocs]

#######################################################################################################################################
    #Edit and Clean the novel covariates
      #Clean Covariates
  #If the covariate applies to both sets, then it should be in the matrix of the longer set
  covnames<-NULL
  C1covs<-NULL
  covars1<-covars1.novel
  if(is.null(covars1.novel))
  {
    C1<-as.matrix(x$covars1)
  } else {
    if(is.null(dim(covars1)))
    {
      C1<-matrix(covars1,n.novels1,nspp2,byrow=FALSE)
      if(is.factor(covars1))
      {
        C1<-as.matrix(as.vector(C1))
        C1covs<-cbind(C1covs,C1)
        C1<-as.matrix(model.matrix(~as.factor(C1)-1)[,-1])
        colnames(C1)<-paste(rep("covar1",length(levels(covars1))-1),levels(covars1)[-1],sep="-")
      } else {
        C1<-as.matrix(as.vector(C1))
        C1covs<-cbind(C1covs,C1)
      }
      covnames<-c(covnames,"covar1")
    } else {
      C1<-NULL
      for(i in 1:dim(covars1)[2])
      {
        C1hold<-matrix(covars1[,i],n.novels1,nspp2,byrow=FALSE)
        if(is.factor(covars1[,i]))
        {
          C1hold<-as.matrix(as.vector(C1hold))
          C1covs<-cbind(C1covs,C1hold)
          C1hold<-as.matrix(model.matrix(~as.factor(C1hold)-1)[,-1])
          colnames(C1hold)<-paste(rep(colnames(covars1)[i],length(levels(covars1[,i]))-1),levels(covars1[,i])[-1],sep="-")
          C1<-cbind(C1,C1hold)
        } else {
          C1hold<-as.matrix(as.vector(C1hold))
          C1covs<-cbind(C1covs,C1hold)
          colnames(C1hold)<-colnames(covars1)[i]
          C1<-cbind(C1,C1hold)
        }
        covnames<-c(covnames,colnames(covars1)[i])
      }
    }
  }

  covnames<-NULL
  C2covs<-NULL
  covars2<-covars2.novel
  if(is.null(covars2.novel))
  {
    C2<-as.matrix(x$covars2)
  } else {
    if(is.null(dim(covars2)))
    {
      C1<-matrix(covars2,n.novels2,nspp1,byrow=FALSE)
      if(is.factor(covars2))
      {
        C2<-as.matrix(as.vector(C2))
        C2covs<-cbind(C2covs,C2)
        C2<-as.matrix(model.matrix(~as.factor(C2)-1)[,-1])
        colnames(C2)<-paste(rep("covar2",length(levels(covars2))-1),levels(covars2)[-1],sep="-")
      } else {
        C2<-as.matrix(as.vector(C2))
        C2covs<-cbind(C2covs,C2)
      }
      covnames<-c(covnames,"covar2")
    } else {
      C2<-NULL
      for(i in 1:dim(covars2)[2])
      {
        C2hold<-matrix(covars2[,i],n.novels2,nspp1,byrow=FALSE)
        if(is.factor(covars2[,i]))
        {
          C2hold<-as.matrix(as.vector(C2hold))
          C2covs<-cbind(C2covs,C2hold)
          C2hold<-as.matrix(model.matrix(~as.factor(C2hold)-1)[,-1])
          colnames(C2hold)<-paste(rep(colnames(covars2)[i],length(levels(covars2[,i]))-1),levels(covars2[,i])[-1],sep="-")
          C2<-cbind(C2,C2hold)
        } else {
          C2hold<-as.matrix(as.vector(C2hold))
          C2covs<-cbind(C2covs,C2hold)
          colnames(C2hold)<-colnames(covars2)[i]
          C2<-cbind(C2,C2hold)
        }
        covnames<-c(covnames,colnames(covars2)[i])
      }
    }
  }

  # Make U, the combined matrix of covariates
  U<-NULL
  U<-cbind(1,C1,C2)

#  if(!is.null(C1))
#  {
#      U<-cbind(U,C1)
#  }
#  if(!is.null(C2))
#  {
#      U<-cbind(U,C2)
#  }
#
#
##  U<-NULL
##  if(is.null(C1) & is.null(C2))
##  {
##    U<-rep(1,length(A))    #need to edit
##  } else {
##
##    if(is.null(C1))
##    {
##      U<-rep(1,length(A))  #need to edit
##    } else {
##      U<-cbind(rep(1,length(A)),C1)
##    }
##
##    if(is.null(C2))
##    {
##    } else {
##      U<-cbind(U,C2)
##    }
##  }

U.novel<-U
  #Set up loop for the first (rows) dimension
  #if(n.novels1>0){
  #cors.1 <- NULL
  #obs.novels1 <- NULL
  #preds.1 <- NULL
  #nams.1 <- NULL
  #se <- NULL
  #other.assocs<-1:length(pairnames)

    #for(i in 1:n.novels1)
    #{

      #novel.assocs<-grep(novel.spp1[i],rownames(V))
      #other.assocs<-other.assocs[-1*novel.assocs]

      Vxy<-V[novel.assocs,orig.assocs]
      invVyy<-qr.solve(Vyy)

      #A amd U matrix
          A.U<-x$variates[rownames(Vyy),]   #variate data for the other assocs and traits with the "novel" removed
          A <-A.U[,1]
          b0<-rep(1,length(orig.assocs))
          U <- as.matrix(cbind(b0,A.U[,-1]))            #might need to remove b0
          #U <- as.matrix(x$fit.variates[-1 * novel.assocs,c(-1,-grep(rownames(x$fit.variates)[1],colnames(x$fit.variates)):-dim(x$fit.variates)[2])])

          #covariate data for novel associations
          #b0<-rep(1,length(novel.assocs))
          #U.novel <- as.matrix(cbind(b0,x$variates[novel.assocs,-1]))

          #mu <- solve((t(U) %*% invVyy %*% U), (t(U) %*% invVyy %*% A))
          mu. <- solve((t(U) %*% invVyy %*% U), (t(U) %*% invVyy %*% A))

          dxgy <- Vxy %*% invVyy %*% (A - (U%*%mu))
          estmu <- U.novel%*%mu + dxgy
          #Vxgy<-V[novel.assocs,novel.assocs]-Vxy%*%invVyy%*%t(Vxy)
          #cors.1<-c(cors.1,cor(estmu,x$variates[novel.assocs,1]))
          #preds.1 <- c(preds.1, estmu)
          #se<-c(se,data.frame(diag(Vxgy)))
          #obs.novels1 <- c(obs.novels1,x$variates[novel.assocs,1])
          #nams.1<-c(nams.1,rownames(estmu))

    if(is.null(covars2.novel))
    {
      assoc.pred<-matrix(estmu,n.novels1,nspp2)
      rownames(assoc.pred)<-novel.spp1
      colnames(assoc.pred)<-sppnames2.orig
    }
    if(is.null(covars1.novel))
    {
      assoc.pred<-matrix(estmu,n.novels2,nspp1)
      rownames(assoc.pred)<-novel.spp2
      colnames(assoc.pred)<-sppnames1.orig
    }

    estmu.point<-estmu
    assoc.pred.point<-assoc.pred
      V1.novel.orig<-V1.novel
      V2.novel.orig<-V2.novel

    if(calcboot)
    {
      nboots<-dim(x$bootvalues)[1]
      mboots<-dim(x$bootvalues)[2]
      #nspp1<-dim(V1.novel.orig)[1]
      #nspp2<-dim(V2.novel.orig)[1]
      estmu.boot <- NULL
    if(n.novels1>0){
    assoc.pred.boot <- array(dim=c(n.novels1,nspp2,nboots))
    rownames(assoc.pred.boot)<-novel.spp1
    colnames(assoc.pred.boot)<-sppnames2.orig
    }
    if(n.novels2>0){
    assoc.pred.boot <- array(dim=c(n.novels2,nspp1,nboots))
    rownames(assoc.pred.boot)<-novel.spp2
    colnames(assoc.pred.boot)<-sppnames1.orig
    }

      for(i in 1:nboots) #loop across the boot values
      {
        mu<-x$bootvalues[i,3:mboots]
        d1<-abs(x$bootvalues[i,1])
        d2<-abs(x$bootvalues[i,2])

        #scale V with novels with the calculated d values

         V1.novel<-V1.novel.orig/det(V1.novel.orig)^(1/nspp1)   # scale covariance matrices (this reduces numerical problems caused by
  	     V2.novel<-V2.novel.orig/det(V2.novel.orig)^(1/nspp2)   # determinants going to infinity or zero)

        # tau = tau_i + tau_j where tau_i equals the node to tip distance
        tau1<-matrix(diag(V1.novel),nspp1,nspp1) + matrix(diag(V1.novel),nspp1,nspp1)-2*V1.novel
  	    tau2<-matrix(diag(V2.novel),nspp2,nspp2) + matrix(diag(V2.novel),nspp2,nspp2)-2*V2.novel

        V1<-(d1^tau1)*(1-d1^(2*V1.novel))/(1-d1^2)
        V2<-(d2^tau2)*(1-d2^(2*V2.novel))/(1-d2^2)
        V1<-V1/det(V1)^(1/nspp1)   # model of coevolution
        V2<-V2/det(V2)^(1/nspp2)
        V<-kronecker(V2,V1,make.dimnames=T)            #the values in this V will be different than the V in x since the novel species is included
        nams<-rownames(V)
        nams<-sub(":","-",nams)
        rownames(V)<-nams
        colnames(V)<-nams
        novel.assocs<-setdiff(rownames(V),pairnames.orig)
        orig.assocs<-pairnames.orig
        Vyy<-V[orig.assocs,orig.assocs]
        Vxy<-V[novel.assocs,orig.assocs]
        invVyy<-qr.solve(Vyy)

      #A amd U matrix
          #A.U<-x$variates[rownames(Vyy),]   #variate data for the other assocs and traits with the "novel" removed
          #A <-A.U[,1]
          #b0<-rep(1,length(orig.assocs))
          #U <- as.matrix(cbind(b0,A.U[,-1]))            #might need to remove b0

          #mu. <- solve((t(U) %*% invVyy %*% U), (t(U) %*% invVyy %*% A))

          dxgy <- Vxy %*% invVyy %*% (A - (U%*%mu))
          estmu. <- U.novel%*%mu + dxgy
          estmu.boot <-cbind(estmu.boot,estmu.)
          #Vxgy<-V[novel.assocs,novel.assocs]-Vxy%*%invVyy%*%t(Vxy)
          #cors.1<-c(cors.1,cor(estmu,x$variates[novel.assocs,1]))
          #preds.1 <- c(preds.1, estmu)
          #se<-c(se,data.frame(diag(Vxgy)))
          #obs.novels1 <- c(obs.novels1,x$variates[novel.assocs,1])
          #nams.1<-c(nams.1,rownames(estmu))

#    assoc.pred.<-matrix(estmu.,n.novels1,nspp2)

    assoc.pred.boot[,,i] <- estmu.
  } #end loop across the boot values

      lob<-list(estmu.point=estmu.point,assoc.pred.point=assoc.pred.point,estmu.boot=estmu.boot,assoc.pred.boot=assoc.pred.boot)
      class(lob)<-"pblm.predict"
      return(lob)

  } else {#end if statement for boot

    lob<-list(estmu.point=estmu.point,assoc.pred.point=assoc.pred.point,estmu.boot=NULL,assoc.pred.boot=NULL)
    class(lob)<-"pblm.predict"
    return(lob)

  }
}


#######################################################################################################################
################ Some Utility functions for pblm objects


####################### Extract Observed values from a pblm
get.obs.pblm<-function(x,sp=NULL)
{

  if(is.null(sp))
  {
  ll<-x$variates$associations
  names(ll)<-rownames(x$variates)
  return(data.frame(ll))
  } else {
  ll.<-NULL
    for(lk in 1:length(sp)){
      ll<-x$variates$associations[grep(sp[lk],rownames(x$variates))]
      names(ll)<-rownames(x$variates)[grep(sp[lk],rownames(x$variates))]
      ll.<-rbind(ll.,ll)
  }
  return(data.frame(t(ll.)))

  }
}

get.predicted.orig<-function(k.orig,sp=NULL)
{
  if(class(k.orig)=="pblm.predict")
  {

    if(is.null(sp))
    {

      return(data.frame(rows=k.orig$preds.1,cols=k.orig$preds.2 ))

    } else {

      jj<-data.frame(rows=k.orig$preds.1,cols=k.orig$preds.2 )
      ll<-NULL
      for(lk in 1:length(sp))
      {
        ll<-rbind(ll,jj[grep(sp[lk],rownames(jj)),])
      }
      return(data.frame(ll))

    }
  }

}


get.predicted.novel<-function(k,rows=TRUE)
{
    bon<-function(uu)
    {
      data.frame(t(uu$assoc.pred.point))
    }

    nams<-function(uu,rows=rows)
    {
      if(rows)
      {
        paste(colnames(uu$assoc.pred.point),rownames(uu$assoc.pred.point),sep="-")
      } else {
        paste(rownames(uu$assoc.pred.point),colnames(uu$assoc.pred.point),sep="-")
      }
    }

  bon(k)->pp
  rownames(pp)<-nams(k,rows=rows)
  return(data.frame(pp))
}


get.predicted.list<-function(tt,rows=TRUE)
{
  bon<-function(uu)
  {
    (uu$assoc.pred.point)
  }
  nams<-function(uu,rows=rows)
  {
    if(rows)
    {
      paste(colnames(uu$assoc.pred.point),rownames(uu$assoc.pred.point),sep="-")
    } else {
      paste(rownames(uu$assoc.pred.point),colnames(uu$assoc.pred.point),sep="-")
    }
  }
  unlist(lapply(tt,bon))->pp
  names(pp)<-unlist(lapply(tt,nams,rows))
  return(data.frame(pp))
}

#Summary pblm

summary.pblm<-function(x)
{
if(class(x)=="pblm")
{

  ans<-list()
  ans$MSE<-x$MSE
  ans$signal<-x$signal.strength
  ans$coef<-x$coefficients[grep("full",rownames(x$coefficients)),]
  ans
}
}
