#I'm not thoroughly testing the output of the quantile, mantel, etc. methods, because those are maintained by other people
# - though perhaps a few worked examples? Might help to check my own code...
#Setup
require(pez)
require(testthat)
require(picante)
data(phylocom)
data <- comparative.comm(phylocom$phylo, phylocom$sample, traits=phylocom$traits, warn=FALSE)

context("fingerprint.regression")

test_that("quantile", {
  set.seed(123)
  expect_warning(basic.quantile <<- fingerprint.regression(data=data, eco.permute=100))
  expect_is(basic.quantile, "fingerprint.regression")
  expect_equal(basic.quantile$eco.method, "quantile")
  expect_equal(basic.quantile$evo.method, "lambda")
  expect_equal(names(basic.quantile), c("evo", "eco", "evo.method", "eco.method"))
  expect_equivalent(basic.quantile$evo, c(1, 1, 1e-06, 1))
  set.seed(123)
  expect_warning(test <<- eco.trait.regression(data, "taxa.labels", 100, "quantile", altogether=FALSE))
  expect_equal(basic.quantile$eco$raw, test)
  set.seed(123)
  expect_equal(basic.quantile$evo, phy.signal(data, method="lambda"))
})

test_that("lm", {
  set.seed(123)
  expect_warning(basic.lm <<- fingerprint.regression(data, eco.permute=100, eco.method="lm", eco.rnd="richness", evo.method="delta"))
  expect_is(basic.lm, "fingerprint.regression")
  expect_equal(basic.lm$eco.method, "lm")
  expect_equal(basic.lm$evo.method, "delta")
  expect_equal(names(basic.lm), c("evo", "eco", "evo.method", "eco.method"))
  set.seed(123)
  expect_warning(test <<- eco.trait.regression(data, "richness", 100, "lm", altogether=FALSE))
  expect_equal(basic.lm$eco$raw, test)
  set.seed(123)
  expect_equal(basic.lm$evo, phy.signal(data, method="delta"))
})

test_that("mantel", {
  set.seed(123)
  basic.mantel <<- fingerprint.regression(data, eco.permute=10, eco.method="mantel", eco.rnd="trialswap", evo.method="kappa")
  expect_is(basic.mantel, "fingerprint.regression")
  expect_equal(basic.mantel$eco.method, "mantel")
  expect_equal(basic.mantel$evo.method, "kappa")
  expect_equal(names(basic.mantel), c("evo", "eco", "evo.method", "eco.method"))
  set.seed(123)
  test <<- eco.trait.regression(data, "trialswap", 10, "mantel", altogether=FALSE)
  expect_equal(basic.mantel$eco$raw, test)
  set.seed(123)
  expect_equal(basic.mantel$evo, phy.signal(data, method="kappa"))
})
