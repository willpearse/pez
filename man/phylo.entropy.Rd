\name{phylo.entropy}
\alias{phylo.entropy}
\title{\code{phylo.entropy} calculate phylogenetic entropy across a comparative.comm object}
\usage{
phylo.entropy(data)
}
\arguments{
  \item{data}{comparative.comm object}
}
\value{
named vector of entropy values for each community in the
input data
}
\description{
\code{phylo.entropy} calculate phylogenetic entropy across
a comparative.comm object
}
\details{
Something and a description of what it does with a
reference!
}
\examples{
\dontrun{
data(pez)
data <- comparative.comm(phy.only$phy, phy.only$sample)
phylo.entropy(data)
}
}

