require(testthat)
require(pez)
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

## TODO:  phylo and func-phylo dists
