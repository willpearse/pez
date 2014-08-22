\name{cc.barplot}
\alias{cc.barplot}
\title{Plots a barplot of a trait in a comparative community object}
\usage{
cc.barplot(data, traits, scaling = FALSE, ranging = TRUE, yranging = NULL,
  tipsontop = TRUE, ceti = 1, csub = 0, f.phylog = 1/(1 + ncol(values)),
  clabel.leaves = 0, grid.col = gray(0.85), connect.col = "black",
  vert.line.col = grey(0.5), vert.line.lty = 1,
  bar.col = heat.colors(dim(values)[1], alpha = 0.8), spacing = 2.5,
  border.col = "black", draw.connecting.lines = FALSE, draw.grid = FALSE,
  draw.vert.line = TRUE, phylog.lwd = 1, ...)
}
\arguments{
  \item{data}{a comparative community ecology object}

  \item{trait}{name of trait you wish to plot (only one!)}
}
\value{
Named numeric vector, where each element is a trait or
community.
}
\description{
\code{cc.barplot} plots a barplot of a trait in a
comparative community object
}
\details{
This may not necessarily be a good thing to do!
}
\examples{
\dontrun{
data(phylocom)
data <- comparative.comm(phylocom$phy, phylocom$sample, traits=phylocom$traits)
cc,barplotc(data, c("traitA", "traitB"))
}
}
\author{
Matt Helmus, Will Pearse
}

