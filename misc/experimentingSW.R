library(pez)
library(picante)
library(ape)
#library(ade4)
#library(subscript)
set.seed(1)

n <- 8 # number of sites
m <- 5 # number of species

env <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
coord <- data.frame(lat = rnorm(n), long = rnorm(n))
rownames(env) <- rownames(coord) <- letters[1:n]

traits <- data.frame(z1 = rnorm(m), z2 = rnorm(m))
rownames(traits) <- LETTERS[1:m]
newick <- "((A,(B,E)),(D,C));"
## tree <- chronos(compute.brlen(read.tree(text = newick)))
tree <- compute.brlen(read.tree(text = newick))
comm <- matrix(rbinom(n*m, 1, 0.8), n, m)
dimnames(comm) <- list(rownames(env), rownames(traits))

cc <- comparative.comm(tree, comm, traits, env)
cc.shape <- shape(cc, removeErrors = FALSE)
cc.shape <- shape(cc, removeErrors = TRUE)


psd(comm, tree)
psd(comm, vcv(tree))

psv(comm, tree)
psv(comm, vcv(tree))

shape(cc, "psv")


plot(cc)
cc.barplot(cc, "z1")

cc[1:2, 1:2]

## cc.barplot(cc, )
