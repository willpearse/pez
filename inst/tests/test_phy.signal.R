require(testthat)
require(pez)
require(caper)
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits, warn=FALSE, vcv=FALSE)

context("phy.signal")

test_that("lambda", {
    set.seed(123)
    lambda <- phy.signal(data, method="lambda")
    expect_that(names(lambda), equals(names(data$data)))
    expect_that(lambda, is_equivalent_to(c(1, 0.294137803455436)))
})

test_that("delta", {
    set.seed(123)
    delta <- phy.signal(data, method="delta")
    expect_that(names(delta), equals(names(data$data)))
    expect_that(delta, is_equivalent_to(c(0.386503004143207, 3)))
})

test_that("kappa", {
    set.seed(123)
    kappa <- phy.signal(data, method="kappa")
    expect_that(names(kappa), equals(names(data$data)))
    expect_that(kappa, is_equivalent_to(c(1e-06, 0.0608586780172654)))
})
