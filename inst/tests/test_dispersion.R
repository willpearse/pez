#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#Testing all the measures together is somewhat silly, because I can't set the set on anything so the randomisations match.
#Setup
#Here's a demo Matt
require(pez)
require(picante)
require(testthat)
data(phylocom)
data <- comparative.comm(phylocom$phylo, phylocom$sample, warn=FALSE)
context("dispersion metrics")

#Tests
test_that("SESmpd", {
  sesmpd_test <<- dispersion(data, "sesmpd")
  expect_that(sesmpd_test$sesmpd$mpd.obs, equals(c(4.85714285714286, 6, 7.14285714285714, 8.28571428571429, 8.85714285714286, 8.42857142857143)))
  expect_that(sesmpd_test$sesmpd$ntaxa, equals(rep(8,6)))
  expect_that(rownames(sesmpd_test$sesmpd), equals(rownames(data$comm)))
  expect_that(sesmpd_test$coefs, equals(coef(sesmpd_test)))
})
test_that("SESmntd", {
  sesmntd_test <<- dispersion(data, "sesmntd")
  expect_that(sesmntd_test$sesmntd$mntd.obs, equals(c(2,2,2,2,6,5)))
  expect_that(sesmntd_test$sesmntd$ntaxa, equals(rep(8,6)))
  expect_that(rownames(sesmntd_test$sesmntd), equals(rownames(data$comm)))
  expect_that(sesmntd_test$coefs, equals(coef(sesmntd_test)))
})
test_that("INND", {
  innd_test <<- dispersion(data, "innd")
  expect_that(innd_test$innd$mpd.obs, equals(c(0.238095238095238, 0.214285714285714, 0.2, 0.164285714285714, 0.116666666666667, 0.12797619047619)))
  expect_that(innd_test$innd$ntaxa, equals(rep(8,6)))
  expect_that(rownames(innd_test$innd), equals(rownames(data$comm)))
  expect_that(innd_test$coefs, equals(coef(innd_test)))
})
#I've checked these values with caper's
test_that("D", {
  set.seed(123)
  d_test <<- dispersion(data, "d")
  expect_that(d_test$d, equals(structure(c(-1.91147244241933, -1.59620673259184, -1.47692430746261, 
-0.184587530399159, 1.8475430319446, 1.54405249539692, 0, 0, 
0, 0.006, 0.997, 0.919, 0.993, 0.991, 0.974, 0.64, 0, 0.003), .Dim = c(6L, 
3L), .Dimnames = list(c("clump1", "clump2a", "clump2b", "clump4", 
"even", "random"), c("D", "P(D=1)", "P(D=0)")))))
  expect_that(d_test$coefs, equals(coef(d_test)))
})

#SESmpd etc. is so unpredictable that tests of the numerical values are almost pointless...
test_that("Non-standard distance matrices",{
    data(laja)
    laja <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
    set.seed(123)
    sqrt <- dispersion(laja, "all", sqrt.phy=TRUE)
    expect_that(names(coef(sqrt)), equals(c("sesmpd","sesmntd","sespd","innd")))
    t <- sqrt(cophenetic(laja$phy))
    set.seed(123)
    ext.dist <- dispersion(laja, "all", ext.dist=as.dist(t))
    expect_that(identical(names(coef(sqrt)),names(coef(ext.dist))), is_false())
    #expect_that(coef(sqrt)[,names(coef(ext.dist))], equals(coef(ext.dist)))
    t <- comparative.comm(invert.tree, river.sites, invert.traits)
    traitgram <- dispersion(t, traitgram=1)
    traitgram.group <- dispersion(t, "all", traitgram=c(0,0.5,1))
    #expect_that(round(coef(traitgram),1), is_equivalent_to(round(traitgram.group[traitgram.group$traitgram==1.0,names(coef(traitgram))],1)))
})
