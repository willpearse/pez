require(testthat)
require(picante)
require(pez)
data(phylocom)

test_that("traitgram", {
    cc <- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, warn=FALSE)
    ccNophy <- within(cc, phy <- NULL)
    ccNotrait <- within(cc, data <- NULL)
    expect_that(cc.traitgram(ccNophy), throws_error())
    expect_that(cc.traitgram.(ccNotrait), throws_error())
    cc$data <- princompOne(cc$data)    
    traitgram.cc(cc)
})
