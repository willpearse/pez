#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#PCD will have to be re-done when you fix the output... remove the NAs and remove the arrays!
#PCD isn't deterministic, so this is hard to test...
#Some aspect of as.dist.default in the 'call' slot of a phylosor and unifrac object is messing up equality testing
#Setup
require(testthat)
data(phylocom)
data <- comparative.comm(phylocom$phylo, phylocom$sample, warn=FALSE)
context("dissimilarity metrics")

#Tests
test_that("UniFrac", {
  unifrac_test <<- dissimilarity(data, "unifrac")
  expect_that(unifrac_test$unifrac, is_equivalent_to(structure(c(0.625, 0.64, 0.8125, 0.789473684210526, 0.771428571428571, 0.653846153846154, 0.607142857142857, 0.763157894736842, 0.742857142857143, 0.571428571428571, 0.736842105263158, 0.714285714285714, 0.470588235294118, 0.675675675675676, 0.538461538461538), Labels = c("clump1", "clump2a", "clump2b", "clump4", "even", "random"), Size = 6L, class = "dist", Diag = FALSE, Upper = FALSE)))
  expect_that(unifrac_test$pcd, equals(NULL))
  expect_that(unifrac_test$phylosor, equals(NULL))
})

test_that("PCD", {
  set.seed(123)
  pcd_test <<- dissimilarity(data, "pcd")
  expect_that(pcd_test$unifrac, equals(NULL))
  expect_that(pcd_test$pcd$PCD, equals(structure(c(NA, NA, NA, NA, NA, NA, 1.19859053329401, NA, NA, NA, NA, NA, 1.12480880494661, 1.09885935603036, NA, NA, NA, NA, 1.60654145062316, 0.983618916257841, 0.922881824502088, NA, NA, NA, 1.49785292594533, 1.35894296255503, 1.26336032683872, 0.692355146728863, NA, NA, 1.49578630826059, 1.34558030091584, 1.23041022411978, 1.19074481768663, 0.918018156083656, NA), .Dim = c(6L, 6L), .Dimnames = list(c("clump1", "clump2a", "clump2b", "clump4", "even", "random"), c("clump1", "clump2a", "clump2b", "clump4", "even", "random")))))
  expect_that(pcd_test$pcd$PCDc, equals(structure(c(NA, NA, NA, NA, NA, NA, 0.735294117647059, NA, NA, NA, NA, NA, 0.735294117647059, 0.735294117647059, NA, NA, NA, NA, 1.10294117647059, 0.735294117647059, 0.735294117647059, NA, NA, NA, 1.10294117647059, 1.10294117647059, 1.10294117647059, 0.735294117647059, NA, NA, 1.10294117647059, 1.10294117647059, 1.10294117647059, 1.10294117647059, 1.10294117647059, NA), .Dim = c(6L, 6L), .Dimnames = list(c("clump1", "clump2a", "clump2b", "clump4", "even", "random"), c("clump1", "clump2a", "clump2b", "clump4", "even", "random")))))
  expect_that(pcd_test$pcd$PCDp, equals(structure(c(NA, NA, NA, NA, NA, NA, 1.63008312527985, NA, NA, NA, NA, NA, 1.52973997472738, 1.49444872420129, NA, NA, NA, NA, 1.45659758189833, 1.33772172611066, 1.25511928132284, NA, NA, NA, 1.35805331952377, 1.23210828604989, 1.14544669633377, 0.941602999551253, NA, NA, 1.35617958615627, 1.2199928061637, 1.11557193653527, 1.07960863470255, 0.832336461515848, NA), .Dim = c(6L, 6L), .Dimnames = list(c("clump1", "clump2a", "clump2b", "clump4", "even", "random"), c("clump1", "clump2a", "clump2b", "clump4", "even", "random")))))
  expect_that(pcd_test$pcd$PSVpool, equals(0.823333333333333))
  expect_that(pcd_test$pcd$PSVmncd, equals(structure(c(0.7487, 0.678963551587302, 0.623292673987885, 0.576011079532113, 0.526413087521158, 0.465702378097296, 0.43409677346687, 0.398520403901477), .Dim = 8L)))
  expect_that(pcd_test$phylosor, equals(NULL))
})

test_that("PhyloSor", {
  phylosor_test <<- dissimilarity(data, "phylosor")
  expect_that(phylosor_test$unifrac, equals(NULL))
  expect_that(phylosor_test$pcd, equals(NULL))
  expect_that(phylosor_test$phylosor, is_equivalent_to(structure(c(0.545454545454545, 0.529411764705882, 0.315789473684211, 0.347826086956522, 0.372093023255814, 0.514285714285714, 0.564102564102564, 0.382978723404255, 0.409090909090909, 0.6, 0.416666666666667, 0.444444444444444, 0.692307692307692, 0.489795918367347, 0.631578947368421), Labels = c("clump1", "clump2a", "clump2b", "clump4", "even", "random"), Size = 6L, class = "dist", Diag = FALSE, Upper = FALSE)))
})


test_that("Each measure is the same as calculated together", {
  set.seed(123)
  all_test <- dissimilarity(data)
  expect_that(all_test$unifrac, equals(unifrac_test$unifrac))
  expect_that(all_test$pcd, equals(pcd_test$pcd))
  expect_that(all_test$phylosor, equals(phylosor_test$phylosor))
})