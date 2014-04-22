\name{eco.phy.regression}
\alias{eco.phy.regression}
\title{Run an eco.phy.regression regression (Cavender-Bares et al. 2004)}
\usage{
eco.phy.regression(data, randomisation = c("taxa.labels", "richness",
  "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  permute = 0, method = c("quantile", "lm", "mantel"), ...)
}
\arguments{
  \item{data}{a comparative community ecology object on
  which to run the regression}

  \item{randomisation}{what kind of null distributions to
  compare your data with - one of: taxa.labels, richness,
  frequency, sample.pool, phylogeny.pool, independentswap,
  trialswap (as implemented in 'picante')}

  \item{permute}{the number of null permutations to
  perform}
}
\description{
\code{eco.phy.regression} runs an eco.phy.regression
regression, in the sense of Cavender-Bares et al. (2004),
on a comparative.comm object.
}
\details{
This is extremely unchcked, so beware!
}
\examples{
\dontrun{
data(phylocom)
data <- comparative.comm(phylocom$phy, phylocom$sample)
eco.phy.regression(data, permute=10)
}
}
\author{
Will Pearse, Jeannine Cavender-Bares
}

