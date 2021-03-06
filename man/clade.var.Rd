% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clade.var.R
\name{clade.var}
\alias{clade.var}
\title{Calculation clade-level beta-dispersion}
\usage{
clade.var(data, alpha = 0.05, clade = NULL)
}
\arguments{
\item{data}{\code{\link{comparative.comm}} object}

\item{alpha}{Critical values (DEFAULT: 5%) for construction of
confidence intervals. See notes for details about this.}

\item{clade}{If given (by default not: NULL) then values for only
these clades will be returned. Should be a numeric vector of
any length; negative values are interpreted as excluding
certain clades. Rarely, if ever, needed.}
}
\value{
a \code{data.frame} with each clade's observed
    \code{variance}, Faith's \code{pd} (calculated by summing the
    branch lengths of each clade), the number of species in the
    clade, the expected variance (from the analytical expectation;
    see 'notes'), and then the upper and lower confidence interval
    bounds for expectation of that variance (see also 'notes';
    calculated at a critical value of \code{alpha}).
}
\description{
As described in Pearse et al. (in press), this function calculations
phylogenetic beta-dispersion for each clade within a
phylogeny. Beta-dispersion is the variation of presences \emph{and
absences} across \emph{multiple} sites, and by focusing on clades'
dispersions patterns it is possible to unpick processes occurring
across a system. This method is based around species' variances in
While this code reports an analytical expectation (and so p-value)
for
}
\note{
As described in Pearse et al. (in press), there are two ways
    of assessing the statistical significance of this approach: an
    analytical expectation and through a null permutation
    test. This code returns the analytical expectation; to perform
    a bootstrap, simply run this function many times across data
    you have permuted in some way. The assumptions implicit in the
    analytical test (e.g., iid and Normally distributed
    abundances/occurrences) will rarely, if ever, be met, and so I
    advise against using it as anything other than a general guide.
}
\examples{
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits)
clade.var(data)
}
\references{
Pearse WD, Legendre P, Peres-Neto, PR, & Davies TJ (in
    press). The interaction of phylogeny and community structure:
    Linking the community composition and trait evolution of
    clades. Global Ecology & Biogeography (available online at
    bioRxiv: 404111).
}
\author{
Will Pearse
}
