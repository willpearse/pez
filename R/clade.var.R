#' Calculation clade-level beta-dispersion
#'
#' As described in Pearse et al. (in press), this function calculations
#' phylogenetic beta-dispersion for each clade within a
#' phylogeny. Beta-dispersion is the variation of presences \emph{and
#' absences} across \emph{multiple} sites, and by focusing on clades'
#' dispersions patterns it is possible to unpick processes occurring
#' across a system. This method is based around species' variances in
#' While this code reports an analytical expectation (and so p-value)
#' for
#' 
#' @param data \code{\link{comparative.comm}} object
#' @param alpha Critical values (DEFAULT: 5%) for construction of
#'     confidence intervals. See notes for details about this.
#' @param clade If given (by default not: NULL) then values for only
#'     these clades will be returned. Should be a numeric vector of
#'     any length; negative values are interpreted as excluding
#'     certain clades. Rarely, if ever, needed.
#' @note As described in Pearse et al. (in press), there are two ways
#'     of assessing the statistical significance of this approach: an
#'     analytical expectation and through a null permutation
#'     test. This code returns the analytical expectation; to perform
#'     a bootstrap, simply run this function many times across data
#'     you have permuted in some way. The assumptions implicit in the
#'     analytical test (e.g., iid and Normally distributed
#'     abundances/occurrences) will rarely, if ever, be met, and so I
#'     advise against using it as anything other than a general guide.
#' @return a \code{data.frame} with each clade's observed
#'     \code{variance}, Faith's \code{pd} (calculated by summing the
#'     branch lengths of each clade), the number of species in the
#'     clade, the expected variance (from the analytical expectation;
#'     see 'notes'), and then the upper and lower confidence interval
#'     bounds for expectation of that variance (see also 'notes';
#'     calculated at a critical value of \code{alpha}).
#' @author Will Pearse
#' @references Pearse WD, Legendre P, Peres-Neto, PR, & Davies TJ (in
#'     press). The interaction of phylogeny and community structure:
#'     Linking the community composition and trait evolution of
#'     clades. Global Ecology & Biogeography (available online at
#'     bioRxiv: 404111).
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' clade.var(data)
#' @importFrom stats qchisq var
#' @importFrom caper clade.matrix
#' @export
clade.var <- function(data, alpha=0.05, clade=NULL){
    #Setup
    ci <- function(x, alpha=0.05){
        upper <- ((length(x)-1)*var(x)) / qchisq(alpha/2, length(x)-1)
        lower <- ((length(x)-1)*var(x)) / qchisq(1-(alpha/2), length(x)-1)
        return(cbind(upper, lower))
    }
    if(!inherits(data, "comparative.comm")) stop("Must use pez::comparative.comm object")
    clade.mat <- clade.matrix(data$phy)$clade.matrix

    variance <- pd <- n.spp <- numeric(nrow(clade.mat))
    cis <- matrix(NA, nrow=nrow(clade.mat), ncol=2)
    for(i in seq_len(nrow(clade.mat))){
        spp <- unname(clade.mat[i,]==1)
        n.spp[i] <- sum(spp)
        variance[i] <- var(rowSums(as.matrix(comm(data)[,spp])))
        cis[i,] <- ci(rowSums(as.matrix(comm(data)[,spp])), alpha)
        if(n.spp[i] == 1) pd[i] <- 0 else pd[i] <- sum(drop_tip(phy(data),setdiff(species(data),species(data)[spp]))$edge.length)/2
    }
    expect.var <- apply(clade.mat, 1, function(x) sum(variance[which(unname(x)==1)]))
    output <- data.frame(cbind(variance, pd, n.spp, expect.var, cis))
    names(output) <- c("variance", "pd", "n.spp", "expect.var", "upper.ci", "lower.ci")
    if(!is.null(clade))
        return(output[clade,])
    return(output)
}
