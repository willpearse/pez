\name{dispersion}
\alias{dispersion}
\title{Calculate dispersion phylogenetic biodiversity metrics across communities}
\usage{
dispersion(data, metric = c("all", "sesmpd", "sesmntd", "sespd", "innd", "d"),
  permute = 1000, ...)
}
\arguments{
  \item{data}{a \code{comparative.comm} object}

  \item{permute}{the number of null permutations to
  perform}

  \item{metric}{specify particular metrics to calculate,
  default is \code{all}}

  \item{...}{additional arguments for modifying the
  behaviour of various metrics (see details)}
}
\value{
a \code{phy.structure} list object of metric values
}
\description{
\code{dispersion} calculates phylogenetic biodiversity
metrics
}
\details{
Calculates various metrics of phylogenetic biodiversity
that are categorized as \emph{dispersion} metrics by Pearse
\emph{et al.} (2014) This calculates $SES_{MPD}$,
$SES_{MNTD}$, $SES_{PD}$, INND, and D. All these are
defined as dispersion metrics in Pearse et al., Dc=0 is the
Brownian expectation, Dc=1 is the random expectation.
}
\examples{
\dontrun{
data(phylocom)
data <- comparative.comm(phylocom$phy, phylocom$sample)
dispersion(data)
dispersion(data, 100, "sesmpd")
}
}
\author{
M.R. Helmus, Will Pearse
}
\references{
Pearse W.D., Purvis A., Cavender-Bares J. & Helmus M.R.
(2014). Metrics and Models of Community Phylogenetics. In:
Modern Phylogenetic Comparative Methods and Their
Application in Evolutionary Biology. Springer Berlin
Heidelberg, pp. 451-464.

\code{sesmpd,sesmntd} Webb C.O. (2000). Exploring the
phylogenetic structure of ecological communities: An
example for rain forest trees. American Naturalist, 156,
145-155.

\code{sespd} Webb C.O., Ackerly D.D. & Kembel S.W. (2008).
Phylocom: software for the analysis of phylogenetic
community structure and trait evolution. Bioinformatics
Applications Note, 24, 2098-2100.

\code{innd} Ness J.H., Rollinson E.J. & Whitney K.D.
(2011). Phylogenetic distance can predict susceptibility to
attack by natural enemies. Oikos, 120, 1327-1334.

\code{d} Fritz S.A. & Purvis A. (2010). Selectivity in
Mammalian Extinction Risk and Threat Types: a New Measure
of Phylogenetic Signal Strength in Binary Traits.
Conservation Biology, 24, 1042-1051.
}

