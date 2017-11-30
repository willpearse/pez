require(testthat)
require(pez)
require(caper)
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits, warn=FALSE)

context("phy.signal")

test_that("lambda", {
    set.seed(123)
    lambda <- phy.signal(data, method="lambda")
    expect_equal(names(lambda), names(data$data))
    expect_equivalent(round(lambda,3), c(1, 0.294))
})

test_that("delta", {
    set.seed(123)
    delta <- phy.signal(data, method="delta")
    expect_equal(names(delta), names(data$data))
    expect_equivalent(round(delta,3), c(0.387, 3))
})

test_that("kappa", {
    set.seed(123)
    kappa <- phy.signal(data, method="kappa")
    expect_equal(names(kappa), names(data$data))
    expect_equivalent(round(kappa,3), c(0,0.061))
})
