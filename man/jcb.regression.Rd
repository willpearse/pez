\name{jcb.regression}
\alias{jcb.regression}
\title{Plot a jcb.plot (sensu Cavender-Bares et al. 2004 Figure 6)}
\usage{
  jcb.regression(data,
    eco.rnd = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
    eco.method = c("quantile", "lm", "mantel"),
    eco.permute = 1000, eco.swap = 1000,
    evo.method = c("lambda", "delta", "kappa"),
    evo.permute = 1000, ...)
}
\arguments{
  \item{data}{a comparative community ecology object on
  which to run the analyses}

  \item{eco.rnd}{null distribution with which to compare
  your community data - one of: taxa.labels, richness,
  frequency, sample.pool, phylogeny.pool, independentswap,
  trialswap (as implemented in 'picante')}

  \item{eco.method}{regression to perform - one of: lm,
  quantile, mantel}

  \item{eco.permute}{the number of null permutations to
  perform}

  \item{eco.swap}{the number of independent swaps to
  perform (if doing so)}

  \item{evo.method}{how to measure phylogenetic inertia -
  one of: lambda, delta, kappa, D.c}

  \item{evo.permute}{number of permutations of phylogenetic
  inertia (D.c only)}
}
\description{
  \code{jcb.regression} calculates traits' phylogenetic
  inertia, and then plots this as a function of trait
  similarity among communities
}
\details{
  This is extremely unchcked, so beware!
}
\examples{
\dontrun{
data(phylocom)
data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
jcb.regression(data, eco.permute=10)
plot(jcb.regression(data, permute=10, method="lm))
plot(jcb.regression(data, permute=10, method="lm", altogether=FALSE))
}
}
\author{
  Will Pearse
}

