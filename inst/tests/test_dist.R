require(testthat)
require(pez)
require(picante)
data(phylocom)

test_that("PA comm dist", {
    cc <- comparative.comm(phylocom$phy,
                           ifelse(phylocom$sample, 1, 0),
                           phylocom$traits, warn=FALSE)
    distByHand <- 1-as.dist(crossprod(cc$comm)/(dim(cc)[1] - crossprod(1-cc$comm)))
    expect_that(comm.dist(cc), equals(distByHand))
})

## TODO:  non-PA dist matrix / overlap warning

test_that("func dist", {
    cc <- comparative.comm(phylocom$phy, phylocom$sample, warn=FALSE)
    expect_error(traits.dist(cc))
    ## TODO:  non-error case
})

test_that("func phylo dist", {
    cc <- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, warn=FALSE)
    aVec <- seq(0, 1, by = 0.2)
    pVec <- 0:4
    PDist <- phylo.dist(cc)
    FDist <- traits.dist(cc)
    for(a in aVec) {
        for(p in pVec) {
            distByHand <- (a*(PDist/max(PDist))^p +
                           (1-a)*(FDist/max(FDist))^p)^(1/p)
            expect_that(funct.phylo.dist(cc, a, p), equals(distByHand))
        }
    }
    PDist[] <- PDist/max(PDist)
    FDist[] <- FDist/max(FDist)
    expect_that(funct.phylo.dist(cc, 1, 2), is_equivalent_to(PDist))
    expect_that(funct.phylo.dist(cc, 0, 2), is_equivalent_to(FDist))
})

test_that("traitgram", {
    cc <- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, warn=FALSE)
    ccNophy <- within(cc, phy <- NULL)
    ccNotrait <- within(cc, data <- NULL)
    expect_that(traitgram.cc(ccNophy), throws_error())
    expect_that(traitgram.cc(ccNotrait), throws_error())
    traitgram.cc(cc)
    ccPCA <- within(cc, data <- princompOne(data))
})
