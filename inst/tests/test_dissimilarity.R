#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#PCD will have to be re-done when you fix the output... remove the NAs and remove the arrays!
#PCD isn't deterministic, so this is hard to test...
#Some aspect of as.dist.default in the 'call' slot of a phylosor and unifrac object is messing up equality testing
#Setup
require(testthat)
require(picante)
require(pez)
data(phylocom)
data <- comparative.comm(phylocom$phylo, phylocom$sample, warn=FALSE)
context("dissimilarity metrics")

#Tests
test_that("UniFrac", {
  unifrac_test <<- dissimilarity(data, "unifrac")
  expect_that(unifrac_test$unifrac, is_equivalent_to(structure(c(0.625, 0.64, 0.8125, 0.789473684210526, 0.771428571428571, 0.653846153846154, 0.607142857142857, 0.763157894736842, 0.742857142857143, 0.571428571428571, 0.736842105263158, 0.714285714285714, 0.470588235294118, 0.675675675675676, 0.538461538461538), Labels = c("clump1", "clump2a", "clump2b", "clump4", "even", "random"), Size = 6L, class = "dist", Diag = FALSE, Upper = FALSE)))
  expect_that(unifrac_test$pcd, equals(NULL))
  expect_that(unifrac_test$phylosor, equals(NULL))
  expect_that(unifrac_test$comdist, equals(NULL))
})

#pcd isn't always the same answer, so these tests are poor...
test_that("PCD", {
  set.seed(123)
  pcd_test <<- dissimilarity(data, "pcd")
  expect_that(pcd_test$unifrac, equals(NULL))
  set.seed(123)
  t <- pcd(data$comm, data$phy)
  expect_that(pcd_test$pcd$PCD, equals(t$PCD, tolerance=0.02))
  expect_that(pcd_test$pcd$PCDc, equals(t$PCDc, tolerance=0.02))
  expect_that(pcd_test$pcd$PCDp, equals(t$PCDp, tolerance=0.02))
  expect_that(pcd_test$phylosor, equals(NULL))
  expect_that(pcd_test$comdist, equals(NULL))
})

test_that("PhyloSor", {
  phylosor_test <<- dissimilarity(data, "phylosor")
  expect_that(phylosor_test$unifrac, equals(NULL))
  expect_that(phylosor_test$pcd, equals(NULL))
  t <- as.dist(1 - as.matrix(phylosor(data$comm, data$phy)))
  #...the call is different between the two...
  expect_that(phylosor_test$phylosor, is_equivalent_to(t))
  expect_that(phylosor_test$comdist, equals(NULL))
})

test_that("comdist", {
  comdist_test <<- dissimilarity(data, "comdist")
  expect_that(comdist_test$unifrac, equals(NULL))
  expect_that(comdist_test$pcd, equals(NULL))
  expect_that(comdist_test$phylosor, equals(NULL))
  expect_that(comdist_test$comdist, equals(comdist(data$comm, cophenetic(data$phy), abundance=TRUE)))
  t <- dissimilarity(data, "comdist", abundance=FALSE)
  expect_that(identical(t, comdist_test), is_false())
})

test_that("Non-standard distance matrices",{
    data(laja)
    laja <<- comparative.comm(invert.tree, river.sites, warn=FALSE)
    sqrt <- dissimilarity(laja, "all", sqrt.phy=TRUE)
    expect_that(names(sqrt), equals(c("unifrac","pcd","phylosor","comdist", "type")))
    t <- sqrt(cophenetic(laja$phy))
    ext.dist <- dissimilarity(laja, "all", ext.dist=as.dist(t))
    expect_that(ext.dist$comdist, equals(sqrt$comdist))
    t <- comparative.comm(invert.tree, river.sites, invert.traits)
    traitgram <- dissimilarity(t, traitgram=1)
})

test_that("Each measure is the same as calculated together", {
  set.seed(123)
  all_test <- dissimilarity(data)
  expect_that(all_test$unifrac, equals(unifrac_test$unifrac))
  #expect_that(all_test$pcd, equals(pcd_test$pcd)) - impossible, see above
  expect_that(all_test$phylosor, equals(phylosor_test$phylosor))
  expect_that(all_test$comdist, equals(comdist_test$comdist))
})
