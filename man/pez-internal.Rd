\name{pez-internal}
\alias{.Random.seed}
\alias{.eco.evo.clade.regression}
\alias{.eco.evo.regression}
\alias{.eco.null}
\alias{.eco.phy.regression}
\alias{.plot.regression}
\alias{.prepare.regression.output}
\alias{.summary.regression}
\alias{summary.phy.structure}
\alias{print.phy.structure}
\alias{coef.phy.structure}
\title{Internal pez Functions}
\alias{.removeErrors}
\description{
 Internal pez functions
}
\details{


  These are not to be called by the user, with the exception of
  \code{print.phy.structure}, \code{summary.phy.structure}, and
  \code{coef.phy.structure}. These three functions are extremely simple,
  and are simply wrappers to get information about these objects
  extracted. \code{coef.phy.structure} is the same as calling
  \code{object$coefs}, as is explained in the documentation for
  \code{\link{shape}}, \code{\link{evenness}}, \code{\link{dispersion}},
  and \code{\link{dissimilarity}}.
}
\keyword{internal}