#This may need checking on other systems, because I'm using exact value checks with no consideration of precision...
#Setup
require(testthat)
require(pez)
require(picante)
data(phylocom)
data <- comparative.comm(phylocom$phylo, phylocom$sample, warn=FALSE)
context("Generic metrics")

test_that("PSV", expect_that(.psv(data), equals(c(0.485714285714286, 0.6, 0.714285714285714, 0.828571428571429, 0.885714285714286, 0.842857142857143))))
test_that("PSR", expect_that(.psr(data), equals(c(3.88571428571429, 4.8, 5.71428571428571, 6.62857142857143, 7.08571428571429, 6.74285714285714))))
test_that("MPD", expect_that(.mpd(data), equals(c(4.85714285714286, 6, 7.14285714285714, 8.28571428571429, 8.85714285714286, 8.42857142857143))))
test_that("MNTD", expect_that(.mntd(data), equals(c(2,2,2,2,6,5))))
test_that("Taxon", expect_that(.taxon(data, abundance.weighted=TRUE), is_equivalent_to(as.data.frame(matrix(c(4.85714285714286, 5.39393939393939, 6.36363636363636, 7.57575757575758, 8.85714285714286, 7.58333333333333, 4.85714285714286, 5.74193548387097, 6.7741935483871, 8.06451612903226, 8.85714285714286, 8.42592592592593, 2.12244897959184, 5.71428571428572, 11.265306122449, 7.3469387755102, 2.12244897959182, 3.81632653061224, 4.85714285714286, 6, 7.14285714285714, 8.28571428571429, 8.85714285714286, 8.42857142857143, 38.8571428571429, 48, 57.1428571428571, 66.2857142857143, 70.8571428571429, 67.4285714285714), ncol=5)))))
test_that("Taxon", expect_that(.taxon(data, abundance.weighted=TRUE), is_equivalent_to(as.data.frame(matrix(c(4.85714285714286, 5.39393939393939, 6.36363636363636, 7.57575757575758, 8.85714285714286, 7.58333333333333, 4.85714285714286, 5.74193548387097, 6.7741935483871, 8.06451612903226, 8.85714285714286, 8.42592592592593, 2.12244897959184, 5.71428571428572, 11.265306122449, 7.3469387755102, 2.12244897959182, 3.81632653061224, 4.85714285714286, 6, 7.14285714285714, 8.28571428571429, 8.85714285714286, 8.42857142857143, 38.8571428571429, 48, 57.1428571428571, 66.2857142857143, 70.8571428571429, 67.4285714285714), ncol=5)))))
test_that("Taxon", expect_that(.taxon(data), is_equivalent_to(as.data.frame(matrix(c(4.85714285714286, 6, 7.14285714285714, 8.28571428571429, 8.85714285714286, 8.42857142857143, 4.85714285714286, 6, 7.14285714285714, 8.28571428571429, 8.85714285714286, 8.42857142857143, 2.12244897959184, 5.71428571428572, 11.265306122449, 7.3469387755102, 2.12244897959182, 3.81632653061224, 4.85714285714286, 6, 7.14285714285714, 8.28571428571429, 8.85714285714286, 8.42857142857143, 38.8571428571429, 48, 57.1428571428571, 66.2857142857143, 70.8571428571429, 67.4285714285714), ncol=5)))))
test_that("PD", expect_that(.pd(data, abundance.weighted=FALSE), is_equivalent_to(cbind(c(16, 17, 18, 22, 30, 27), c(-4.88235294117647, -4.82352941176471, -3.82352941176471, 0.176470588235294, 9.11764705882352, 4.23529411764706)))))
test_that("PD", expect_that(.pd(data, abundance.weighted=TRUE), is_equivalent_to(cbind(c(16, 17, 18, 22, 30, 27), c(-5.66666666666666, -4.66666666666667, -3.66666666666667, 0.333333333333332, 8.33333333333333, 5.33333333333333)))))
test_that("Colless", expect_that(.colless(data), is_equivalent_to(c(0, 0, 0, 0, 0, 5))))
test_that("Gamma", expect_that(.gamma(data), is_equivalent_to(c(-1.41421356237309, -0.707106781186547, -0.157134840263678, -0.385694607919935, -2.9227080289044, -2.19988776369148))))
test_that("Eigenvectors", {
    expect_that(.eigen.sum(data), is_equivalent_to(c(1.85989806057967e-32, 9.80257655161354e-05, 0.0580509448862977, 0.0456171567606935, 0.0434666292330911, 0.0476659578856296)))
  expect_that(.eigen.sum(data, which.eigen=2), is_equivalent_to(c(3.30159419037812e-32, 0.082086748214039, 0.0137252988148037, 0.0415764004365078, 0.0384905400587625, 0.0452932152543931)))
})
test_that("EED", expect_that(.eed(data), equals(c(3.202317, 3.104838, 3.013283, 2.697896, 2.241960, 2.391608), tolerance=0.00001)))
test_that("HED", expect_that(.hed(data), equals(setNames(c(6.659032, 6.456330, 6.265945, 5.610116, 4.662026, 4.973210),rownames(data$comm)), tolerance=0.00001)))
test_that("Rao", expect_that(.rao(data), is_equivalent_to(c(2.125, 2.47222222222222, 2.91666666666667, 3.47222222222222, 3.875, 3.5546875))))
test_that("Entropy", expect_that(.phylo.entropy(data), is_equivalent_to(c(2.07944154167984, 4.62549821485909, 5.2620123831539, 6.64830674427379, 8.31776616671934, 7.22716285080812))))
test_that("PAE", expect_that(.pae(data), is_equivalent_to(c(-0.3125, 0, 0.0555555555555556, 0.227272727272727, 0.366666666666667, 0.666666666666667))))
test_that("IAC", expect_that(.iac(data), is_equivalent_to(c(0.5, 0.75, 0.75, 0.6875, 0.416666666666667, 0.916666666666667))))
test_that("Haed", expect_that(.haed(data), is_equivalent_to(c(c(2.07944154167984, 2.43261844647745, 2.42601513195981, 2.42601513195981, 2.07944154167984, 2.67765404685526)))))
test_that("Eaed", expect_that(.eaed(data), is_equivalent_to(c(1, 0.978957679027899, 0.976300309778954, 0.976300309778954, 1, 0.965759553653586))))
set.seed(123); test_that("Delta", expect_that(.delta(data), equals(c(NA, 1e-06, 1e-06, 1e-06, NA, 3))))
set.seed(123); test_that("Kappa", expect_that(.kappa(data), equals(c(NA, 3, 3, 3, NA, 1.65281212634147))))
set.seed(123); test_that("Lambda", expect_that(.lambda(data), is_equivalent_to(c(NA, 1, 1, 1, NA, 1e-06))))


#Dispersion
test_that("SESmpd", expect_that(.ses.mpd(data)[,2], equals(c(4.85714285714286, 6, 7.14285714285714, 8.28571428571429, 8.85714285714286, 8.42857142857143))))
test_that("SESmntd", expect_that(.ses.mntd(data)[,2], equals(c(2,2,2,2,6,5))))
set.seed(123);test_that("D", expect_that(.d(data), equals(structure(c(-1.91147244241933, -1.59620673259184, -1.47692430746261, -0.184587530399159, 1.8475430319446, 1.54405249539692, 0, 0, 
0, 0.006, 0.997, 0.919, 0.993, 0.991, 0.974, 0.64, 0, 0.003), .Dim = c(6L, 
3L), .Dimnames = list(c("clump1", "clump2a", "clump2b", "clump4", 
"even", "random"), c("D", "P(D=1)", "P(D=0)"))))))

#Dissimilarity
test_that("UniFrac", expect_that(.unifrac(data), is_equivalent_to(structure(c(0.625, 0.64, 0.8125, 0.789473684210526, 0.771428571428571, 0.653846153846154, 0.607142857142857, 0.763157894736842, 0.742857142857143, 0.571428571428571, 0.736842105263158, 0.714285714285714, 0.470588235294118, 0.675675675675676, 0.538461538461538), Labels = c("clump1", "clump2a", "clump2b", "clump4", "even", "random"), Size = 6L, class = "dist", Diag = FALSE, Upper = FALSE))))
#...PCD is not reproducible, so no test :-(
#...comdist has no good test either :-(









if(FALSE){          
          test_that("Non-standard distance matrices",{
    sqrt <- shape(laja, "all", sqrt.phy=TRUE)
    expect_that(names(sqrt), equals(c("psv","psr","mpd","mntd","pd","pd.ivs","colless","gamma","taxon","eigen.sum","Eed","Hed","dist.fd","type","coefs")))
    t <- sqrt(cophenetic(laja$phy))
    ext.dist <- shape(laja, "all", ext.dist=as.dist(t))
    expect_that(identical(names(coef(sqrt)),names(coef(ext.dist))), is_false())
    expect_that(coef(sqrt)[,names(coef(ext.dist))], equals(coef(ext.dist)))
    t <- comparative.comm(invert.tree, river.sites, invert.traits)
    traitgram <- shape(t, traitgram=1)
    traitgram.group <- shape(t, "all", traitgram=c(0,0.5,1))
    expect_that(coef(traitgram), is_equivalent_to(traitgram.group[traitgram.group$traitgram==1.0,names(coef(traitgram))]))
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

                    test_that("MPD", expect_that(.mpd(data, abdundance.weighted=TRUE), equals(c(c(4.25, 4.94444444444444, 5.83333333333333, 6.94444444444444, 7.75, 7.109375)))))
test_that("MNTD", expect_that(.mntd(data, abdundance.weighted=TRUE), equals(mntd(data$comm, cophenetic(data$phy), abundance.weighted=TRUE))))


#SESmpd etc. is so unpredictable that tests of the numerical values are almost pointless...
test_that("Non-standard distance matrices",{
    data(laja)
    laja <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
    set.seed(123)
    sqrt <- dispersion(laja, "all", sqrt.phy=TRUE)
    expect_that(names(coef(sqrt)), equals(c("sesmpd","sesmntd","sespd","innd")))
    t <- sqrt(cophenetic(laja$phy))
    set.seed(123)
    ext.dist <- dispersion(laja, "all", ext.dist=as.dist(t))
    expect_that(identical(names(coef(sqrt)),names(coef(ext.dist))), is_false())
    #expect_that(coef(sqrt)[,names(coef(ext.dist))], equals(coef(ext.dist)))
    t <- comparative.comm(invert.tree, river.sites, invert.traits)
    traitgram <- dispersion(t, traitgram=1)
    traitgram.group <- dispersion(t, "all", traitgram=c(0,0.5,1))
    #expect_that(round(coef(traitgram),1), is_equivalent_to(round(traitgram.group[traitgram.group$traitgram==1.0,names(coef(traitgram))],1)))
})

test_that("INND", expect_that(innd_test$innd$mpd.obs, equals(c(0.238095238095238, 0.214285714285714, 0.2, 0.164285714285714, 0.116666666666667, 0.12797619047619)))
          test_that("PhyloSor", expect_that(.phylosor(data), is_equivalent_to(as.numeric(c(0.454545454545455, 0.470588235294118, 0.684210526315789, 0.652173913043478, 0.627906976744186, 0.485714285714286, 0.435897435897436, 0.617021276595745, 0.590909090909091, 0.4, 0.583333333333333, 0.555555555555556, 0.307692307692308, 0.510204081632653, 0.368421052631579), Labels = c("clump1", "clump2a", "clump2b", "clump4", "even", "random"), Size = 6L, call = as.dist.default(m = 1 - as.matrix(output)), class = "dist", Diag = FALSE, Upper = FALSE))))

          
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

      }
