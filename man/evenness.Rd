\name{evenness}
\alias{evenness}
\title{Calculate evenness metrics across communities}
\usage{
evenness(data, metric = c("all", "rao", "taxon", "entropy", "cadotte", "pst",
  "lambda", "delta", "kappa"))
}
\arguments{
  \item{data}{a comparative community ecology object}

  \item{metric}{specify (a) particular metric(s) to
  calculate (all, rao, delta, entropy, cadotte, pst), or
  the default 'all'}
}
\value{
cc.evenness object (a named list with the output from each
metric)
}
\description{
\code{evenness} calculates evenness metrics in
comparative.comm communities
}
\details{
This calculates Rao's quadratic entropy, taxonomic
diversity, phylogenetic entropy, Cadotte's abundance
measures, and Pst, all defined as evenness metrics in
Pearse et al.
}
\examples{
\dontrun{
data(phylocom)
data <- comparative.comm(phylocom$phy, phylocom$sample)
evenness(data)
evenness(data, "rao")
}
}
\author{
Matt Helmus, Will Pearse
}

