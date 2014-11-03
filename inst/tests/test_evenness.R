#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#Rather strange to test entropy when you know it's wrong!...
#Setup
require(testthat)
require(picante)
require(pez)
require(caper)
data(phylocom)
data <- comparative.comm(phylocom$phylo, phylocom$sample, warn=FALSE)
context("evenness metrics")

#Tests
test_that("Rao", {
  rao_test <<- evenness(data, "rao")
  expect_that(rao_test$rao, is_equivalent_to(c(2.125, 2.47222222222222, 2.91666666666667, 3.47222222222222, 3.875, 3.5546875)))
  expect_that(names(rao_test$rao), equals(rownames(data$comm)))
    for(i in seq_along(rao_test))
      if(!names(rao_test)[i] %in% c("type","coefs","rao"))
          expect_that(rao_test[[i]], equals(NULL))
  expect_that(rao_test$rao, is_equivalent_to(rao_test$coefs$rao))
})

test_that("taxon", {
  set.seed(123)
  taxon_test <<- evenness(data, "taxon")
  expect_that(taxon_test$taxon$Species, is_equivalent_to(c(8, 8, 8, 8, 8, 8)))
  expect_that(names(taxon_test$taxon$Species), equals(rownames(data$comm)))
  expect_that(taxon_test$taxon, equals(taxondive(data$comm, cophenetic(data$phy))))
  expect_that(taxon_test$coefs$taxon.Delta, is_equivalent_to(taxon_test$taxon$D[1:6]))
  expect_that(taxon_test$coefs$taxon.DeltaStar, is_equivalent_to(taxon_test$taxon$Dstar[1:6]))
  expect_that(taxon_test$coefs$taxon.LambdaPlus, is_equivalent_to(taxon_test$taxon$Lambda[1:6]))
  expect_that(taxon_test$coefs$taxon.DeltaPlus, is_equivalent_to(taxon_test$taxon$Dplus[1:6]))
  expect_that(taxon_test$coefs$taxon.S.DeltaPlus, is_equivalent_to(taxon_test$taxon$SDplus[1:6]))
  for(i in seq_along(taxon_test))
      if(!names(taxon_test)[i] %in% c("type","coefs","taxon"))
          expect_that(taxon_test[[i]], equals(NULL))
})

test_that("Entropy", {
  entropy_test <<- evenness(data, "entropy")
  expect_that(entropy_test$entropy, is_equivalent_to(c(2.07944154167984, 4.62549821485909, 5.2620123831539, 
6.64830674427379, 8.31776616671934, 7.22716285080812)))
  expect_that(names(entropy_test$entropy), equals(rownames(data$comm)))
  expect_that(entropy_test$entropy, is_equivalent_to(entropy_test$coefs$entropy))
  for(i in seq_along(entropy_test))
      if(!names(entropy_test)[i] %in% c("type","coefs","entropy"))
          expect_that(entropy_test[[i]], equals(NULL))
})

test_that("PAE", {
  pae <<- evenness(data, "pae")
  for(i in seq_along(pae))
      if(!names(pae)[i] %in% c("type","coefs","pae"))
          expect_that(pae[[i]], equals(NULL))
  expect_that(pae$pae, is_equivalent_to(c(-0.3125, 0, 0.0555555555555556, 0.227272727272727, 0.366666666666667, 0.666666666666667)))
})

test_that("IAC", {
  iac <<- evenness(data, "iac")
  for(i in seq_along(iac))
      if(!names(iac)[i] %in% c("type","coefs","iac"))
          expect_that(iac[[i]], equals(NULL))
  expect_that(iac$iac, is_equivalent_to(c(0.5, 0.75, 0.75, 0.6875, 0.416666666666667, 0.916666666666667)))
})

test_that("Haed", {
  haed <<- evenness(data, "haed")
  for(i in seq_along(haed))
      if(!names(haed)[i] %in% c("type","coefs","Haed"))
          expect_that(haed[[i]], equals(NULL))
  expect_that(haed$Haed, is_equivalent_to(c(c(2.07944154167984, 2.43261844647745, 2.42601513195981, 2.42601513195981, 2.07944154167984, 2.67765404685526))))
})

test_that("Eaed", {
  eaed <<- evenness(data, "eaed")
  for(i in seq_along(eaed))
      if(!names(eaed)[i] %in% c("type","coefs","Eaed"))
          expect_that(eaed[[i]], equals(NULL))
    expect_that(eaed$Eaed, is_equivalent_to(c(1, 0.978957679027899, 0.976300309778954, 0.976300309778954, 1, 0.965759553653586)))

})
  
test_that("Delta", {
  set.seed(123)
  delta_test <<- evenness(data, "delta")
  expect_that(delta_test$delta$values, equals(c(NA, 1e-06, 1e-06, 1e-06, NA, 3)))
  expect_that(delta_test$delta$models[[1]], equals(NA))
  expect_that(coef(delta_test$delta$models[[2]]), is_equivalent_to(1.5))
  expect_that(delta_test$delta$values, is_equivalent_to(delta_test$coefs$delta))
  for(i in seq_along(delta_test))
      if(!names(delta_test)[i] %in% c("type","coefs","delta"))
          expect_that(delta_test[[i]], equals(NULL))
})

test_that("Kappa", {
  set.seed(123)
  kappa_test <<- evenness(data, "kappa")
  expect_that(kappa_test$kappa$values, equals(c(NA, 3, 3, 3, NA, 1.65281212634147)))
  expect_that(kappa_test$kappa$models[[1]], equals(NA))
  expect_that(coef(kappa_test$kappa$models[[2]]), is_equivalent_to(1.5))
  expect_that(kappa_test$kappa$values, is_equivalent_to(kappa_test$coefs$kappa))
  for(i in seq_along(kappa_test))
      if(!names(kappa_test)[i] %in% c("type","coefs","kappa"))
          expect_that(kappa_test[[i]], equals(NULL))
})

test_that("Lambda", {
  set.seed(123)
  lambda_test <<- evenness(data, "lambda")
  expect_that(lambda_test$lambda$values, equals(c(NA, 1, 1, 1, NA, 1e-06)))
  expect_that(lambda_test$lambda$models[[1]], equals(NA))
  expect_that(coef(lambda_test$lambda$models[[2]]), is_equivalent_to(1.5))
  expect_that(lambda_test$lambda$values, is_equivalent_to(lambda_test$coefs$lambda))
  for(i in seq_along(lambda_test))
      if(!names(lambda_test)[i] %in% c("type","coefs","lambda"))
          expect_that(lambda_test[[i]], equals(NULL))
})

test_that("MPD", {
  mpd_test <<- evenness(data, "mpd")
  expect_that(mpd_test$mpd, is_equivalent_to(mpd_test$coefs$mpd))
  expect_that(mpd_test$mpd, equals(c(c(4.25, 4.94444444444444, 5.83333333333333, 6.94444444444444, 7.75, 7.109375))))
  for(i in seq_along(mpd_test))
      if(!names(mpd_test)[i] %in% c("type","coefs","mpd"))
          expect_that(mpd_test[[i]], equals(NULL))
})

test_that("MNTD", {
  mntd_test <<- evenness(data, "mntd")
  expect_that(mntd_test$mntd, equals(mntd(data$comm, cophenetic(data$phy), abundance.weighted=TRUE)))
  for(i in seq_along(mntd_test))
      if(!names(mntd_test)[i] %in% c("type","coefs","mntd"))
          expect_that(mntd_test[[i]], equals(NULL))
})

test_that("dist FD", {
    data(laja)
    laja <<- comparative.comm(invert.tree, river.sites, warn=FALSE)
    #dist.fd <<- evenness(laja, "dist.fd")
    # - tested below; same problems as in shape (not enough info etc.)
    #for(i in seq_along(dist.fd))
    #    if(!names(dist.fd)[i] %in% c("type","coefs"))
    #        expect_that(dist.fd[[i]], equals(NULL))
})

test_that("Non-standard distance matrices",{
    sqrt <- evenness(laja, "all", sqrt.phy=TRUE)
    expect_that(names(coef(sqrt)), equals(c("rao","taxon.Delta","taxon.DeltaStar","taxon.LambdaPlus","taxon.DeltaPlus","taxon.S.DeltaPlus","entropy","pae","iac","Haed","Eaed","mpd","mntd","FRic","FEve","FDiv","FDis","RaoQ")))
    t <- sqrt(cophenetic(laja$phy))
    ext.dist <- evenness(laja, "all", ext.dist=as.dist(t))
    expect_that(identical(names(coef(sqrt)),names(coef(ext.dist))), is_false())
    expect_that(coef(sqrt)[,names(coef(ext.dist))], equals(coef(ext.dist)))
    t <- comparative.comm(invert.tree, river.sites, invert.traits)
    traitgram <- evenness(t, traitgram=1)
    traitgram.group <- evenness(t, "all", traitgram=c(0,0.5,1))
    expect_that(coef(traitgram), is_equivalent_to(traitgram.group[traitgram.group$traitgram==1.0,names(coef(traitgram))]))
})
