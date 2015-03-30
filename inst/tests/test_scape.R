require(testthat)
require(pez)

context("scape")
expect_that(names(eco.scape(rcoal(64))), equals(c("cc", "gradient1", "gradient2", "X.joint", "X1", "X2", "nichewd", "K", "environ", "sppXs", "V.phylo", "V.phylo.rho", "V.center", "bsp1", "bspp2")))
expect_that(names(scape(rcoal(64))), equals(c("cc", "X.joint", "X1", "X2", "sppXs", "V.phylo", "V.phylo.rho", "V.center", "V.range", "V.repulse", "bspp1", "bspp2", "u", "wd")))
            
