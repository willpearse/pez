\name{plot.fpcomSims}
\alias{plot.fpcomSims}
\title{Plot fpcomSims objects}
\usage{
\method{plot}{fpcomSims}(x, y, cex.add = 0.8, plottype = c("distance",
  "traitgram", "gradient", "ordination"), ...)
}
\arguments{
  \item{x}{An \code{fpcomSims} object.}

  \item{y}{Not used at the moment.}

  \item{cex.add}{The cex graphics parameter for additional
  graphical.  elements not produced by the
  \code{\link{plot}} command.}

  \item{plottype}{The type of plot to produce -- see
  details.}

  \item{...}{Additional parameters to be passed to
  \code{\link{plot}}.}
}
\value{
No return value, called for its side-effect of producing a
plot.
}
\description{
This plot method can produce four different graphical
summaries of the output of the \code{\link{fpcomSims}}
function.
}
\details{
\describe{

\item{\code{plottype} equals \code{"distance"}}{A
scatterplot of the phylogenetic versus functional distances
between the species pairs are produced. Species numbers are
used to label the points. The \code{\link{cophenetic}}
function is used to compute the phylogenetic distances and
the \code{\link{dist}} function is used to compute the
functional Euclidean distances. Both distances are
standardised such that the maximum distance is one.}

\item{\code{plottype} equals \code{"traitgram"}}{A
\code{\link{traitgram}} is produced, which combines both
the phylogeny and the observed trait.}

\item{\code{plottype} equals \code{"gradient"}}{A
scatterplot of the probabilities of occurrence versus the
observed gradient is produced.  The species are identified
by their numbers. Note that these probabilities of
occurrence depend only on the two gradients and the two
traits, and therefore are not subject to variation among
sites with identical gradient values. Such variation is
expressed in the \code{comm} element of \code{fpcomSims}
objects, which gives not probability of occurrence but
occurrence itself. Note that if the \code{p} argument to
\code{fpcomSims} is one, then only the observed trait and
gradient determine probability of occurrence and the plot
consists of perfect sigmoid curves.  But if \code{p} is
zero then only the unknown gradient and trait determine
probability of occurrence and the plot is just noise. In
this latter case we would expect phylogenetic distance to
provide more important information about communities.}

\item{\code{plottype} equals \code{"ordination"}}{A special
type of ordination of the species is produced. This
ordination is based on the probabilities of occurrence and
not on occurrence itself, and so it is able to fully
explain all variation on two axes (if \code{p} is not
either zero of one). Therefore, such an ordination is not
possible in practice with real data but it is useful in
this case for fully representing distances between species
in species distribution space. Technically it is a singular
value decomposition of the logit transformed probabilities
of occurrence. Note well the percentages of variation
explained by the two axes, and that this ordination is
gauranteed to have 100 percent of the variation explained
by the first two axes.} }
}

