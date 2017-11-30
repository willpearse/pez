require(testthat)
require(picante)
require(pez)
data(phylocom)
data(laja)
test_that("traitgram", {
    cc <- comparative.comm(phylocom$phy, phylocom$sample, phylocom$traits, warn=FALSE)
    ccNophy <- within(cc, phy <- NULL)
    ccNotrait <- within(cc, data <- NULL)
    expect_error(cc.traitgram(ccNophy))
    expect_error(cc.traitgram.(ccNotrait))
    cc$data <- princompOne(cc$data)    
    traitgram.cc(cc)
})


test_that("ConDivSim", {
    pez <- comparative.comm(invert.tree, river.sites, invert.traits)
    pez <- pez[,colSums(pez$comm) > 0]
    pez$comm[pez$comm>1] <- 1
    expect_equal(colnames(ConDivSim(pez)), c("SpRich", "ExpMeanMPD", "ExpQ0.005MPD", "ExpQ0.025MPD", "ExpQ0.05MPD", "ExpQ0.25MPD", "ExpQ0.75MPD", "ExpQ0.95MPD", "ExpQ0.975MPD", "ExpQ0.995MPD", "ExpSdMPD"))
    expect_equal(colnames(ConDivSim(pez, "phy")), c("SpRich", "ExpMeanMPD", "ExpQ0.005MPD", "ExpQ0.025MPD", "ExpQ0.05MPD", "ExpQ0.25MPD", "ExpQ0.75MPD", "ExpQ0.95MPD", "ExpQ0.975MPD", "ExpQ0.995MPD", "ExpSdMPD"))
    expect_equal(colnames(ConDivSim(pez, as.dist(cophenetic(pez$phy)))), c("SpRich", "ExpMeanMPD", "ExpQ0.005MPD", "ExpQ0.025MPD", "ExpQ0.05MPD", "ExpQ0.25MPD", "ExpQ0.75MPD", "ExpQ0.95MPD", "ExpQ0.975MPD", "ExpQ0.995MPD", "ExpSdMPD"))
})

