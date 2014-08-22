\name{D.c}
\alias{D.c}
\title{Calculate D(c) across communities}
\usage{
D.c(data, permute = 1000, traits = TRUE)
}
\arguments{
  \item{data}{a comparative community ecology object}

  \item{permute}{the number of null permutations to
  perform}

  \item{traits}{if TRUE (default) calculate phylogenetic
  signal of species' traits; otherwise calculate it of
  community abundance}
}
\value{
Matrix with columns for each trait or community, and rows
for Dc value, and P(Dc=0) and P(Dc=1).
}
\description{
\code{D.c} runs the extension of D (Dc) across a
comparative community object
}
\details{
This is heavily based on David Orme's implementation of
'phylo.d', which is in turn heavily based on Susanne
Fritz's implementation of 'phylo.d'. There is no 'D'
function, because Dc calculated across binary data should
be equivalent to D. Please don't use this without
contacting Will Pearse first - this isn't published, and
I'd like to be given some credit for it! Dc=0 is the
Brownian expectation, Dc=1 is the random expectation.
}
\examples{
\dontrun{
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits)
D.c(data, permute=1000)
D.c(data, traits=FALSE)
}
}
\author{
Will Pearse
}

