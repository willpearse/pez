#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#Setup
if(FALSE){
require(testthat)
require(pez)
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits, warn=FALSE)
context("Metric wrappers")

#Shape
test_that("Shape",{
    expect_warning(expect_equivalent(pez.shape(data), cbind(.pd(data,abundance.weighted=FALSE),.psv(data,abundance.weighted=FALSE),.psr(data,abundance.weighted=FALSE),.mpd(data,abundance.weighted=FALSE),.mntd(data,abundance.weighted=FALSE),.vpd(data,abundance.weighted=FALSE),.vntd(data,abundance.weighted=FALSE),.mipd(data,abundance.weighted=FALSE),.innd(data,abundance.weighted=FALSE),.taxon(data,abundance.weighted=FALSE),.eigen.sum(data,abundance.weighted=FALSE),.eed(data,abundance.weighted=FALSE),.hed(data,abundance.weighted=FALSE),.scheiner(data,abundance.weighted=FALSE))))
    t <- .sqrt.phy(data)
    expect_warning(expect_equivalent(pez.shape(t,sqrt.phy=TRUE), cbind(.pd(t,abundance.weighted=FALSE),.mpd(t,abundance.weighted=FALSE),.mntd(t,abundance.weighted=FALSE),.vpd(data,abundance.weighted=FALSE),.vntd(data,abundance.weighted=FALSE),.mipd(t,abundance.weighted=FALSE),.innd(t,abundance.weighted=FALSE),.taxon(t,abundance.weighted=FALSE),.eigen.sum(t,abundance.weighted=FALSE),.scheiner(t,abundance.weighted=FALSE))))
    expect_true({pez.shape(data,traitgram=0.5);TRUE})
    expect_warning(expect_equivalent(pez.shape(data,traitgram=c(0,1))[,1:10], rbind(pez.shape(data,traitgram=0),pez.shape(data,traitgram=1))))
    expect_warning(expect_equivalent(pez.shape(data,ext.dist=as.dist(cophenetic(data$phy))), cbind(.mpd(data),.mntd(data),.mipd(data),.innd(data),.taxon(data),.eigen.sum(data))))
})

#Evenness
test_that("Evenness",{
    expect_equivalent(pez.evenness(data), cbind(.rao(data,abundance.weighted=TRUE),.phylo.entropy(data,abundance.weighted=TRUE),.pae(data,abundance.weighted=TRUE),.iac(data,abundance.weighted=TRUE),.haed(data,abundance.weighted=TRUE),.eaed(data,abundance.weighted=TRUE),.mpd(data,abundance.weighted=TRUE),.mntd(data,abundance.weighted=TRUE),.mipd(data,abundance.weighted=TRUE),.innd(data,abundance.weighted=TRUE),.taxon(data,abundance.weighted=TRUE),.pse(data,abundance.weighted=TRUE),.scheiner(data,abundance.weighted=TRUE)))
    t <- .sqrt.phy(data)
    expect_equivalent(pez.evenness(data,sqrt.phy=TRUE), cbind(.rao(t,abundance.weighted=TRUE),.phylo.entropy(t,abundance.weighted=TRUE),.pae(t,abundance.weighted=TRUE),.iac(t,abundance.weighted=TRUE),.haed(t,abundance.weighted=TRUE),.eaed(t,abundance.weighted=TRUE),.mpd(t,abundance.weighted=TRUE),.mntd(t,abundance.weighted=TRUE),.mipd(t,abundance.weighted=TRUE),.innd(t,abundance.weighted=TRUE),.taxon(t,abundance.weighted=TRUE),.scheiner(t,abundance.weighted=TRUE)))
    expect_true({pez.evenness(data,traitgram=0.5);TRUE})
    expect_equivalent(pez.evenness(data,traitgram=c(0,1))[,1:11], rbind(pez.evenness(data,traitgram=0),pez.evenness(data,traitgram=1)))
    expect_equivalent(pez.evenness(data,ext.dist=as.dist(cophenetic(data$phy))), cbind(.mpd(data,abundance.weighted=TRUE),.mntd(data,abundance.weighted=TRUE),.mipd(data,abundance.weighted=TRUE),.innd(data,abundance.weighted=TRUE),.taxon(data,abundance.weighted=TRUE)))
})

#Dispersion
# - there are no meaningful tests of this, because I can't fix the value of SESmpd etc.
test_that("Dispersion",{
    expect_true({pez.dispersion(data, permute=10);TRUE})
    expect_true({pez.dispersion(data, abundance=TRUE, permute=10);TRUE})
    expect_true({pez.dispersion(data, traitgram=0.5, permute=10);TRUE})
    expect_true({pez.dispersion(data, traitgram=c(0,1), permute=10);TRUE})
    expect_true({pez.dispersion(data, sqrt.phy=TRUE, permute=10);TRUE})
})

#Dissimilarity
# - pointless doing anything smart here given that none of these (even with set.seed) are deterministic
test_that("Dissimilarity",{
    expect_equal(names(pez.dissimilarity(data,permute=10)), c("unifrac","pcd","phylosor","comdist"))
    expect_equal(names(pez.dissimilarity(data, abundance.weighted=TRUE,permute=10)), "comdist")
    expect_equal(names(pez.dissimilarity(data,traitgram=0.5,permute=10)), "comdist")
    expect_equal(names(pez.dissimilarity(data,ext.dist=as.dist(cophenetic(data$phy)),permute=10)), "comdist")
})

#Endemisim
test_that("Endemism",{
    expect_equivalent(pez.endemism(data), data.frame(PE=.pe(data),BED=.bed(data)))
})

#Generic wrappers
test_that("generic.metrics", {
    expect_equal(generic.metrics(data, c(.mpd, .pse, .rao)), cbind(.mpd(data), .pse(data), .rao(data)))
    expect_equal(generic.metrics(data, c(.mpd, .pse, .rao), abundance.weighted=TRUE), cbind(.mpd(data, abundance.weighted=TRUE), .pse(data, abundance.weighted=TRUE), .rao(data, abundance.weighted=TRUE)))
})
test_that("generic.null", {
    null <- generic.null(data, c(.mpd, .pse))
    expect_equal(dim(null), c(11,2,5))
    expect_equal(dimnames(null), list(NULL, NULL, c("observed","null.mean","SE","SES","rank")))
})
}
