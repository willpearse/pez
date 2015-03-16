#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#Setup
require(testthat)
require(pez)
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits, warn=FALSE)
context("Metric wrappers")

#Shape
test_that("Shape",{
    expect_that(shape(data), is_equivalent_to(cbind(.psv(data,abundance.weighted=FALSE),.psr(data,abundance.weighted=FALSE),.mpd(data,abundance.weighted=FALSE),.mntd(data,abundance.weighted=FALSE),.taxon(data,abundance.weighted=FALSE),.eigen.sum(data,abundance.weighted=FALSE),.eed(data,abundance.weighted=FALSE),.hed(data,abundance.weighted=FALSE))))
    t <- .sqrt.phy(data)
    expect_that(shape(data,sqrt.phy=TRUE), is_equivalent_to(cbind(.mpd(t,abundance.weighted=FALSE),.mntd(t,abundance.weighted=FALSE),.taxon(t,abundance.weighted=FALSE),.eigen.sum(t,abundance.weighted=FALSE))))
    expect_that(round(shape(data,traitgram=0.5),2), is_equivalent_to(as.data.frame(matrix(c(0.54, 0.63, 0.5, 0.56, 0.52, 0.59, 0.51, 0.45, 0.53, 0.48, 0.52, 0.23, 0.39, 0.28, 0.3, 0.32, 0.27, 0.23, 0.27, 0.3, 0.26, 0.32, 0.54, 0.63, 0.5, 0.56, 0.52, 0.59, 0.51, 0.45, 0.53, 0.48, 0.52, 0.54, 0.63, 0.5, 0.56, 0.52, 0.59, 0.51, 0.45, 0.53, 0.48, 0.52, 0.03, 0.02, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.03, 0.02, 0.03, 0.54, 0.63, 0.5, 0.56, 0.52, 0.59, 0.51, 0.45, 0.53, 0.48, 0.52, 15.55, 5.07, 11.97, 10.06, 6.18, 12.43, 14.19, 8.08, 8.43, 9.67, 7.86, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01),ncol=8))))
    expect_that(shape(data,traitgram=c(0,1))[,1:8], is_equivalent_to(rbind(shape(data,traitgram=0),shape(data,traitgram=1))))
    expect_that(shape(data,ext.dist=as.dist(cophenetic(data$phy))), is_equivalent_to(cbind(.mpd(data),.mntd(data),.taxon(data),.eigen.sum(data))))
})

#Evenness
test_that("Evenness",{
    expect_that(evenness(data), is_equivalent_to(cbind(.rao(data,abundance.weighted=TRUE),.phylo.entropy(data,abundance.weighted=TRUE),.pae(data,abundance.weighted=TRUE),.iac(data,abundance.weighted=TRUE),.haed(data,abundance.weighted=TRUE),.eaed(data,abundance.weighted=TRUE),.mpd(data,abundance.weighted=TRUE),.mntd(data,abundance.weighted=TRUE),.taxon(data,abundance.weighted=TRUE),.pse(data,abundance.weighted=TRUE))))
    t <- .sqrt.phy(data)
    expect_that(evenness(data,sqrt.phy=TRUE), is_equivalent_to(cbind(.rao(t,abundance.weighted=TRUE),.phylo.entropy(t,abundance.weighted=TRUE),.pae(t,abundance.weighted=TRUE),.iac(t,abundance.weighted=TRUE),.haed(t,abundance.weighted=TRUE),.eaed(t,abundance.weighted=TRUE),.mpd(t,abundance.weighted=TRUE),.mntd(t,abundance.weighted=TRUE),.taxon(t,abundance.weighted=TRUE))))
    expect_that(round(evenness(data,traitgram=0.5),2), is_equivalent_to(as.data.frame(matrix(c(37.95, 12.58, 20.3, 10.01, 10.74, 3.86, 22.19, 1.26, 459.76, 21.97, 26.3, 0.74, 0.63, 0.82, 0.67, 0.78, 0.84, 0.72, 0.85, 0.6, 0.7, 0.68, 0.21, 0.06, 0.34, 0.13, 0.31, 0.51, 0.28, 0.31, 0.02, 0.25, 0.14, 0.23, 0.24, 0.2, 0.24, 0.29, 0.3, 0.2, 0.22, 0.23, 0.18, 0.23, 0.21, 0.06, 0.34, 0.13, 0.31, 0.51, 0.28, 0.32, 0.02, 0.25, 0.14, 0.51, 0.53, 0.43, 0.39, 0.46, 0.65, 0.4, 0.45, 0.47, 0.43, 0.45, 0.03, 0.02, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.03, 0.02, 0.03, 0.54, 0.63, 0.5, 0.56, 0.52, 0.59, 0.51, 0.45, 0.53, 0.48, 0.52, 15.55, 5.07, 11.97, 10.06, 6.18, 12.43, 14.19, 8.08, 8.43, 9.67, 7.86),ncol=9))))
    expect_that(evenness(data,traitgram=c(0,1))[,1:9], is_equivalent_to(rbind(evenness(data,traitgram=0),evenness(data,traitgram=1))))
    expect_that(evenness(data,ext.dist=as.dist(cophenetic(data$phy))), is_equivalent_to(cbind(.mpd(data,abundance.weighted=TRUE),.mntd(data,abundance.weighted=TRUE),.taxon(data,abundance.weighted=TRUE))))
})

#Dispersion
# - do when you've fixed INND etc.

#Dissimilarity
# - pointless doing anything smart here given that none of these (even with set.seed) are deterministic
test_that("Dissimilarity",{
    expect_that(names(dissimilarity(data,permute=10)), equals(c("unifrac","pcd","phylosor","comdist")))
    expect_that(names(dissimilarity(data, abundance.weighted=TRUE,permute=10)), equals("comdist"))
    expect_that(names(dissimilarity(data,traitgram=0.5,permute=10)), equals("comdist"))
    expect_that(names(dissimilarity(data,ext.dist=as.dist(cophenetic(data$phy)),permute=10)), equals("comdist"))            
})
