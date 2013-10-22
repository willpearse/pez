#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#Testing all the measures together is somewhat silly, because I can't set the set on anything so the randomisations match.
#Setup
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
})
test_that("SESmntd", {
  sesmntd_test <<- dispersion(data, "sesmntd")
  expect_that(sesmntd_test$sesmntd$mntd.obs, equals(c(2,2,2,2,6,5)))
  expect_that(sesmntd_test$sesmntd$ntaxa, equals(rep(8,6)))
  expect_that(rownames(sesmntd_test$sesmntd), equals(rownames(data$comm)))
})
test_that("INND", {
  innd_test <<- dispersion(data, "innd")
  expect_that(innd_test$innd$mpd.obs, equals(c(0.238095238095238, 0.214285714285714, 0.2, 0.164285714285714, 0.116666666666667, 0.12797619047619)))
  expect_that(innd_test$innd$ntaxa, equals(rep(8,6)))
  expect_that(rownames(innd_test$innd), equals(rownames(data$comm)))
})
test_that("D", {
  set.seed(123)
  d_test <<- dispersion(data, "d")
  expect_that(d_test$d, is_equivalent_to(matrix(c(0.245391114795942, 0.0553744718496124, 0.296182917482368, 0.393634696092272, 1.57306330183479, 0.598846457288309, 0.05, 0.023, 0.053, 0.077, 0.915, 0.169, 0.355, 0.476, 0.349, 0.262, 0.001, 0.178), nrow=6, ncol=3)))
  expect_that(rownames(d_test$d), equals(rownames(data$comm)))
  expect_that(colnames(d_test$d), equals(c("D.c","P(D.c=1)","P(D.c=0)")))
})
#test_that("Each measure is the same as calculated together", {
#  set.seed(123)
#  all_test <- dispersion(data)
#  expect_that(sesmpd_test$sesmpd, equals(all_test$sesmpd))
#  expect_that(sesmntd_test$sesmntd, equals(all_test$sesmntd))
#  expect_that(innd_test$innd, equals(all_test$innd))
#  expect_that(d_test$d, equals(all_test$d))
#})