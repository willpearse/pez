\name{eco.trait.regression}
\alias{eco.trait.regression}
\title{Run an eco.trait.regression regression (Cavender-Bares et al. 2004)}
\usage{
  eco.trait.regression(data,
    randomisation = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
    permute = 0, method = c("quantile", "lm", "mantel"),
    swap.iter = 1000, altogether = TRUE, ...)
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

  \item{altogether}{use distance matrix based on all traits
  (default TRUE), or perform separate regressions for each
  trait}
}
\description{
  \code{eco.trait.regression} runs a regression of
  ecological community coexistence versus species
  trait-distances
}
\details{
  This is extremely unchcked, so beware!
}
\note{
  Like the eco.trait and eco.env methods, this is a
  data-hungry method. Warnings will be generated if any of
  the methods cannot be fitted properly (the examples below
  give toy examples of this). In such cases the summary and
  plot methods of these functions may generate errors; use
  'traceback()' to examine where these are coming from, and
  consider whether you want to be working with the data
  generating these errors. I am loathe to hide these errors
  or gloss over them, because they represent the reality of
  your data!
}
\examples{
\dontrun{
data(phylocom)
data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
eco.trait.regression(data, permute=10)
plot(eco.trait.regression(data, permute=10, method="lm))
plot(eco.trait.regression(data, permute=10, method="lm", altogether=FALSE))
}
}
\author{
  Will Pearse, Jeannine Cavender-Bares
}

