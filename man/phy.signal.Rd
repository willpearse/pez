\name{phy.signal}
\alias{phy.signal}
\title{Calculate a lambda, delta, or kappa for a comparative community ecology dataset}
\usage{
phy.signal(data, traits = TRUE, method = c("lambda", "delta", "kappa"))
}
\arguments{
  \item{data}{a comparative community ecology object on
  which you want to calculate phylogenetic signal}

  \item{traits}{if TRUE (default) calculate phylogenetic
  signal of species' traits; otherwise calculate it of
  community abundance}
}
\value{
Named numeric vector, where each element is a trait or
community.
}
\description{
\code{classic.phy.signal} calculates Pagel's lambda, delta,
or kappa for all community abundances or trait values.
Essentially a wrapper around
}
\details{
This may not necessarily be a good thing to do! You may not
have trait or community data where these measures make
sense!
}
\examples{
\dontrun{
data(phylocom)
data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
classic.phy.signal(data)
classic.phy.signal(data, traits=FALSE, method="kappa")
}
}
\author{
Will Pearse, Jeannine Cavender-Bares
}
\references{
R. P. Freckleton, P. H. Harvey, and M. Pagel. Phylogenetic
analysis and comparative data: A test and review of
evidence. American Naturalist, 160:712-726, 2002.
}

