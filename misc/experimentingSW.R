library(pez)
library(ape)
set.seed(1)

n <- 8 # number of sites
m <- 5 # number of species

env <- data.frame(x = rnorm(n))
coord <- data.frame(lat = rnorm(n), long = rnorm(n))
rownames(env) <- rownames(coord) <- letters[1:n]

traits <- data.frame(z1 = rnorm(m), z2 = rnorm(m))
rownames(traits) <- LETTERS[1:m]
newick <- "((A,(B,E)),(D,C));"
tree <- chronos(compute.brlen(read.tree(text = newick)))
comm <- matrix(rbinom(n*m, 1, 0.5), n, m)
dimnames(comm) <- list(rownames(env), rownames(traits))

comparative.comm(tree, comm, traits, env)
