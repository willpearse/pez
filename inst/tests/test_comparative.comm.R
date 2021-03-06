require(testthat)
require(picante)
require(pez)
require(caper)
data(phylocom)

context("comparative community object")

test_that("comm - phy", {
  simple <<- comparative.comm(phylocom$phy, phylocom$sample, warn=FALSE)
  expect_equal(simple, simple[])
  expect_equal(simple, simple[,])
  expect_equal(ncol(simple$comm), ncol(phylocom$sample))
  expect_equal(nrow(simple$comm), nrow(phylocom$sample))
  expect_equal(simple$data, NULL)
  expect_equal(simple$env, NULL)
  expect_equal(simple$comm[3,], phylocom$sample[3,simple$phy$tip.label])
  expect_equal(simple[1:3,]$comm[3,], simple$comm[3,])
  expect_equal(simple[,1:3]$comm[1,], simple$comm[1,1:3])
  expect_equal(sort(colnames(simple$comm)), sort(simple$phy$tip.label))
  expect_equal(simple[1:2,], simple[c("clump1","clump2a")])
  expect_equal(simple[1:2,], simple[c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)])
  expect_equal(simple[,1:10], simple[,species(simple)[1:10]])
  expect_equal(simple[,1:10], simple[,c(rep(TRUE,10),rep(FALSE,50))])

  t <- simple
  expect_equal(species(t), colnames(simple$comm))
  expect_equal(species(t), simple$phy$tip.label)
  species(t) <- letters[1:25]
  expect_equal(species(t), letters[1:25])
  species(t)[1:3] <- c("w", "t", "f")
  expect_equal(species(t), c("w","t","f",letters[4:25]))
  expect_equal(sites(t), rownames(simple$comm))
  sites(t) <- letters[1:6]
  expect_equal(sites(t), letters[1:6])
  sites(t)[1:3] <- c("w", "t", "f")
  expect_equal(sites(t), c("w","t","f",letters[4:6]))
  expect_true(is.null(env.names(simple)))
  expect_true(is.null(trait.names(simple)))
  
  t <- phylocom$sample
  colnames(t)[1] <- "deliberately wrong"
  dropping <- comparative.comm(phylocom$phy, t, warn=FALSE)
  expect_equal(ncol(dropping$comm), ncol(simple$comm)-1)
  expect_equal(nrow(dropping$comm), nrow(simple$comm))
  expect_equal(names(simple), c("phy","comm","data","env","dropped","names"))
  expect_equal(rownames(dropping$vcv), colnames(dropping$vcv))
  expect_true(is.null(simple$data))
  expect_true(is.null(simple$env))
})

test_that("comm - phy - traits", {
  trts <<- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, warn=FALSE)
  expect_equal(ncol(trts$comm), ncol(simple$comm))
  expect_equal(nrow(trts$comm), nrow(simple$comm))
  expect_equal(trts$comm[3,], simple$comm[3,])
  expect_equal(trts$data[3,], phylocom$traits[3,])
  expect_equal(trts$env, NULL)
  expect_equal(trts[1:3,]$comm[3,], simple[1:3,]$comm[3,])
  expect_equal(trts[,1:3]$comm[1,], simple[,1:3]$comm[1,])
  expect_equal(trts[1:3,]$data[3,], trts$data[3,])
  expect_equal(trts[,1:2]$data[1,], trts$data[1,])
  expect_equal(sort(colnames(trts$comm)), sort(trts$phy$tip.label))
  expect_equal(sort(rownames(trts$data)), sort(trts$phy$tip.label))
  expect_equal(names(trts), names(simple))
  expect_true(is.null(trts$env))
  expect_equal(trts[1:2,], trts[c("clump1","clump2a")])
  expect_equal(trts[1:2,], trts[c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)])
  expect_equal(trts[,1:10], trts[,species(trts)[1:10]])
  expect_equal(trts[,1:10], trts[,c(rep(TRUE,10),rep(FALSE,50))])
})

test_that("comm - phy - traits - env", {
  dummy.env <- data.frame(letters=letters[1:6], repeated=rep("A", 6), row.names=rownames(phylocom$sample))
  env <- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, dummy.env, warn=FALSE)
  expect_equal(ncol(env$comm), ncol(simple$comm))
  expect_equal(nrow(env$comm), nrow(simple$comm))
  expect_equal(env$comm[3,], simple$comm[3,])
  expect_equal(env$data[3,], phylocom$traits[3,])
  expect_equal(env[1:3,]$comm[3,], simple[1:3,]$comm[3,])
  expect_equal(env[,1:3]$comm[1,], simple[,1:3]$comm[1,])
  expect_equal(env[1:3,]$data[3,], trts$data[3,])
  expect_equal(env[,1:2]$data[1,], trts$data[1,])
  expect_equal(env[1:3,]$env[,1], dummy.env[1:3,1])
  expect_equal(env[,1:3]$env[,1], dummy.env[,1])
  expect_equal(sort(colnames(env$comm)), sort(env$phy$tip.label))
  expect_equal(rownames(env$comm), rownames(env$env))
  expect_equal(sort(rownames(env$data)), sort(env$phy$tip.label))
  expect_equal(names(env), names(simple))
  expect_equal(env[1:2,], env[c("clump1","clump2a")])
  expect_equal(env[1:2,], env[c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)])
  expect_equal(env[,1:10], env[,species(env)[1:10]])
  expect_equal(env[,1:10], env[,c(rep(TRUE,10),rep(FALSE,50))])
  
  t <- env
  expect_equal(species(t), colnames(env$comm))
  expect_equal(species(t), env$phy$tip.label)
  expect_equal(sites(t), rownames(env$env))
  expect_equal(species(t), rownames(env$data))
  species(t) <- letters[1:25]
  expect_equal(species(t), letters[1:25])
  species(t)[1:3] <- c("wa", "ta", "fa")
  expect_equal(species(t), c("wa","ta","fa",letters[4:25]))
  expect_equal(sites(t), rownames(env$comm))
  expect_equal(sites(t), rownames(env$env))
  sites(t) <- letters[1:6]
  expect_equal(sites(t), letters[1:6])
  sites(t)[1:3] <- c("w", "t", "fa")
  expect_equal(sites(t), c("w","t","fa",letters[4:6]))
  expect_equal(env.names(env), colnames(env$env))
  expect_equal(trait.names(env), colnames(env$data))

  traits(t)$test <- letters[nrow(traits(t))]
  expect_equal(trait.names(t), c(trait.names(env),"test"))
  env(t)$test <- letters[nrow(comm(t))]
  expect_equal(names(env(t)), c("letters","repeated","test"))
  expect_equal(comm(t), t$comm)
  tt <- t
  expect_warning(comm(t) <- comm(t)[-2,])
  expect_equal(comm(t), tt$comm[-2,])
  expect_equal(phy(t), t$phy)
  expect_equal(phy(t), tree(t))
 })
