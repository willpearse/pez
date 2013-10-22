require(testthat)
require(picante)
data(phylocom)

context("comparative community object")

test_that("comm - phy", {
  simple <<- comparative.comm(phylocom$phy, phylocom$sample, warn=FALSE)
  expect_that(simple, equals(simple[]))
  expect_that(simple, equals(simple[,]))
  expect_that(ncol(simple$comm), equals(ncol(phylocom$sample)))
  expect_that(nrow(simple$comm), equals(nrow(phylocom$sample)))
  expect_that(simple$traits, equals(NULL))
  expect_that(simple$env, equals(NULL))
  expect_that(simple$comm[3,], equals(phylocom$sample[3,]))
  expect_that(simple[1:3,]$comm[3,], equals(simple$comm[3,]))
  expect_that(simple[,1:3]$comm[1,], equals(simple$comm[1,1:3]))
  
  t <- phylocom$sample
  colnames(t)[1] <- "deliberately wrong"
  dropping <- comparative.comm(phylocom$phy, t, warn=FALSE)
  expect_that(ncol(dropping$comm), equals(ncol(simple$comm)-1))
  expect_that(nrow(dropping$comm), equals(nrow(simple$comm)))
})

test_that("comm - phy - traits", {
  traits <<- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, warn=FALSE)
  expect_that(ncol(traits$comm), equals(ncol(simple$comm)))
  expect_that(nrow(traits$comm), equals(nrow(simple$comm)))
  expect_that(traits$comm[3,], equals(simple$comm[3,]))
  expect_that(traits$traits[3,], equals(phylocom$traits[3,]))
  expect_that(traits$env, equals(NULL))
  expect_that(traits[1:3,]$comm[3,], equals(simple[1:3,]$comm[3,]))
  expect_that(traits[,1:3]$comm[1,], equals(simple[,1:3]$comm[1,]))
  expect_that(traits[1:3,]$traits[3,], equals(traits$traits[3,]))
  expect_that(traits[,1:2]$traits[1,], equals(traits$traits[1,]))
})

test_that("comm - phy - traits - env", {
  dummy.env <- data.frame(letters=letters[1:6], repeated=rep("A", 6), row.names=rownames(phylocom$sample))
  env <<- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, dummy.env, warn=FALSE)
  expect_that(ncol(env$comm), equals(ncol(simple$comm)))
  expect_that(nrow(env$comm), equals(nrow(simple$comm)))
  expect_that(env$comm[3,], equals(simple$comm[3,]))
  expect_that(env$traits[3,], equals(phylocom$traits[3,]))  
  expect_that(env[1:3,]$comm[3,], equals(simple[1:3,]$comm[3,]))
  expect_that(env[,1:3]$comm[1,], equals(simple[,1:3]$comm[1,]))
  expect_that(env[1:3,]$traits[3,], equals(traits$traits[3,]))
  expect_that(env[,1:2]$traits[1,], equals(traits$traits[1,]))
  expect_that(env[1:3,]$env[,1], equals(dummy.env[1:3,1]))
  expect_that(env[,1:3]$env[,1], equals(dummy.env[,1]))
})