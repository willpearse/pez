\name{dissimilarity}
\alias{dissimilarity}
\title{Calculate dissimilarity metrics across communities}
\usage{
  dissimilarity(data,
    metric = c("all", "unifrac", "pcd", "phylosor"),
    pa = TRUE, permute = 1000)
}
\arguments{
  \item{data}{a comparative community ecology object}

  \item{metric}{specify (a) particular metric(s) to
  calculate (unifrac, pcd, phylosor), or the default 'all'}

  \item{pa}{If TRUE (default), all metrics are calculated
  across presence-absence matrices, and species abundances
  are ignored}

  \item{permute}{Number of permutations for metric
  (currently only for PCD)}
}
\value{
  cc.dissimilarity object (a named list with the output
  from each metric)
}
\description{
  \code{dissimilarity} calculates dissimilarity metrics in
  comparative.comm communities
}
\details{
  This calculates UniFrac, PCD (and its subcomponents), and
  the phylogenetic Sorenson's index, all defined as
  dissimilarity metrics in Pearse et al.
}
\note{
  This function uses XXX's version of the PCD function, not
  that included in picante
}
\examples{
\dontrun{
data(phylocom)
data <- comparative.comm(phylocom$phy, phylocom$sample)
dissimilarity(data)
dispersion(data, "pcd")
}
}
\author{
  Matt Helmus, Will Pearse
}

