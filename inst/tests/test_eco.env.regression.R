#I'm not thoroughly testing the output of the quantile, mantel, etc. methods, because those are maintained by other people
#This needs proper examples, because at the moment the second column of environmental data are a nonsense (doesn't affect tests, but you will need decent examples!)
#Setup
require(testthat)
require(picante)
data(phylocom)
dummy.env <- data.frame(letters=rep(letters[1:3],2), row.names=rownames(phylocom$sample), stringsAsFactors=TRUE)
data <- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, dummy.env, warn=FALSE)

context("eco.env.regression")

test_that("quantile", {
  set.seed(123)
  basic.quantile <<- eco.env.regression(data)
  set.seed(123)
  expect_equal(basic.quantile$method, "quantile")
  expect_equal(basic.quantile$method, eco.env.regression(data, method="quantile")$method)
  expect_equal(names(basic.quantile), c("observed", "randomisations", "obs.slope", "rnd.slopes", "method", "permute", "randomisation", "type", "altogether", "data"))
  expect_equal(basic.quantile$permute, 0)
  expect_equal(basic.quantile$altogether, TRUE)
  expect_equal(basic.quantile$method, "quantile")
  expect_equal(basic.quantile$randomisations, list())
  expect_equal(basic.quantile$obs.slope, setNames(0.5830952,"as.numeric(trait.mat)"))
  expect_equal(basic.quantile$data, data)
  expect_equal(basic.quantile$type, "eco.env.regression")
  expect_equal(class(basic.quantile), "eco.xxx.regression")
})

test_that("mantel", {
  set.seed(123)
  basic.mantel <<- eco.env.regression(data, method="mantel")
  expect_equal(names(basic.mantel), c("observed", "randomisations", "obs.slope", "rnd.slopes", "method", "permute", "randomisation", "type", "altogether", "data"))
  expect_equal(basic.mantel$permute, basic.quantile$permute)
  expect_equal(basic.mantel$altogether, basic.quantile$altogether)
  expect_equal(basic.mantel$method, "mantel")
  expect_equal(basic.mantel$randomisations, basic.quantile$randomisations)
  expect_equal(basic.mantel$obs.slope, 0.5382758)
  expect_equal(basic.mantel$data, data)
  expect_equal(basic.mantel$type, "eco.env.regression")
  expect_equal(class(basic.mantel), "eco.xxx.regression")
  
  set.seed(123)
  complex.mantel <<- eco.env.regression(data, method="mantel", randomisation="frequency", permute=10, altogether=FALSE)
  expect_equal(names(complex.mantel), c("", "altogether", "type", "data", "permute", "method"))
  expect_equal(complex.mantel$permute, 10)
  expect_equal(complex.mantel$altogether, FALSE)
  expect_equal(complex.mantel$method, "mantel")
  expect_equal(complex.mantel$data, data)
  expect_equal(names(complex.mantel[[1]]), c("observed", "randomisations", "obs.slope", "rnd.slopes", "method", "permute", "randomisation", "type", "altogether"))
  expect_equal(complex.mantel[[1]]$observed$statistic, 0.5382758)
  expect_equal(round(complex.mantel[[1]]$randomisations[[9]]$statistic,4), 0.0058)
  expect_equal(complex.mantel$type, "eco.env.regression.list")
  expect_equal(class(complex.mantel), "eco.xxx.regression.list")
})

test_that("lm", {
  set.seed(123)
  basic.lm <<- eco.env.regression(data, method="lm")
  set.seed(123)
  expect_equal(names(basic.lm), names(basic.quantile))
  expect_equal(basic.lm$permute, basic.quantile$permute)
  expect_equal(basic.lm$altogether, basic.quantile$altogether)
  expect_equal(basic.lm$method, "lm")
  expect_equal(basic.lm$randomisations, basic.lm$randomisations)
  expect_equal(round(basic.lm$obs.slope,4), setNames(0.6194,"as.numeric(trait.mat)"))
  expect_equal(basic.lm$data, data)
  expect_equal(basic.lm$type, "eco.env.regression")
  expect_equal(class(basic.lm), "eco.xxx.regression")
  
  set.seed(123)
  complex.lm <- eco.env.regression(data, method="lm", randomisation="independentswap", permute=10, altogether=FALSE)
  expect_equal(round(coef(complex.lm[[1]]$observed),4), setNames(c(0.0144, 0.6194),c("(Intercept)", "as.numeric(trait.mat)")))
  expect_equal(names(complex.lm), names(complex.mantel))
  expect_equal(complex.lm$permute, complex.mantel$permute)
  expect_equal(complex.lm$altogether, complex.mantel$altogether)
  expect_equal(complex.lm$method, "lm")
  expect_equal(complex.lm$data, data)
  expect_equal(names(complex.lm[[1]]), c("observed", "randomisations", "obs.slope", "rnd.slopes", "method", "permute", "randomisation", "type", "altogether"))
  expect_equal(round(coef(complex.lm[[1]]$randomisations[[5]]),4), setNames(c(0.1516,0.1137),c("(Intercept)","as.numeric(trait.mat)")))
  expect_equal(class(complex.lm), "eco.xxx.regression.list")
})
