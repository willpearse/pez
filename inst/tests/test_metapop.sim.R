require(testthat)
require(pez)
context("Metapopulation and phylogeny simulation")

#Tests
# - no explicit test of edge2phylo because these test it
test_that("sim.meta.comm", expect_equal(names(sim.meta.comm()), c("comm","environment")))
test_that("sim.meta.phy.comm", expect_equal(class(sim.meta.phy.comm()), c("comparative.comm", "comparative.data")))
