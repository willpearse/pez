#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#Rather strange to test entropy when you know it's wrong!...
#Setup
require(testthat)
data(phylocom)
data <- comparative.comm(phylocom$phylo, phylocom$sample, warn=FALSE)
context("phylo.entropy")

#Tests
test_that("It works...", {
  output <- phylo.entropy(data)
  expect_that(names(output), equals(rownames(data$comm)))
  expect_that(output, is_equivalent_to(c(4.15888308335967, 4.62549821485909, 5.2620123831539, 6.64830674427379, 8.31776616671934, 7.22716285080812)))
})