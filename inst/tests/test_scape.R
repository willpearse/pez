require(testthat)
require(pez)

context("scape")
expect_equal(names(eco.scape(rcoal(64))), c("cc", "Y", "Yab", "gradient1", "gradient2", "X.joint", "X1", "X2", "nichewd", "K", "environ", "sppXs", "V.phylo", "V.phylo.rho", "V.center", "bsp1", "bspp2"))
expect_equal(names(scape(rcoal(64))), c("cc", "Y", "X.joint", "X1", "X2", "sppXs", "V.phylo", "V.phylo.rho", "V.center", "V.range", "V.repulse", "bspp1", "bspp2", "u", "wd"))
expect_equal(names(scape(rcoal(64), same.range=FALSE)), c("cc", "Y", "X.joint", "X1", "X2", "sppXs", "V.phylo", "V.phylo.rho", "V.center", "V.range", "V.repulse", "bspp1", "bspp2", "u", "wd"))
            
