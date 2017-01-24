pez.richness <- function(data, abundance.weighted=FALSE, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, which.eigen=1, quick=TRUE, q=0.0001)
{
    ## Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
        stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
    if(!is.null(traitgram)){
        if(length(traitgram) > 1){
            output <- vector("list", length(traitgram))
            for(i in seq_along(output))
                output[[i]] <- cbind(Recall(data, sqrt.phy, traitgram=traitgram[i], which.eigen=which.eigen, traitgram.p=traitgram.p, q=q), traitgram[i], sites(data))
            output <- do.call(rbind, output)
            names(output)[ncol(output)-1] <- "traitgram"
            names(output)[ncol(output)] <- "site"
            rownames(output) <- NULL
            return(output)
        } else {
            dist <- as.matrix(funct.phylo.dist(data, traitgram, traitgram.p))
            traitgram <- TRUE
        }      
    } else traitgram <- FALSE
    
    if(!is.null(ext.dist)){
        dist <- .check.ext.dist(ext.dist, species(data), ncol(data$comm))
        ext.dist <- TRUE
    } else ext.dist <- FALSE
    
    ## Setup
    if(!abundance.weighted)
        data$comm[data$comm > 0] <- 1
    coefs <- data.frame(row.names=rownames(data$comm))
    if(sqrt.phy)
        data <- .sqrt.phy(data)
    if(traitgram==FALSE & ext.dist==FALSE)
        dist <- cophenetic(data$phy)
    
    ## Filter metrics according to suitability and calculate
    functions <- setNames(c(.pd,.psr), c("pd","psr"))
    if(sqrt.phy == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "colless", "gamma", "eed", "hed")]
    if(traitgram == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "pd", "colless", "gamma", "eed", "hed", "scheiner")]
    if(ext.dist == TRUE)
        functions <- functions[!names(functions) %in% c("pd", "psv", "psr", "colless", "gamma", "eed", "hed", "scheiner")]
    if(!is.binary.tree(data$phy) & "colless" %in% names(functions))
        warning("Cannot compute Colless' index with non-binary tree")
    output <- lapply(functions, function(x) try(x(data, dist=dist, abundance.weighted=abundance.weighted, which.eigen=which.eigen), silent=TRUE))
    
    ## Clean up output and return
    output <- Filter(function(x) !inherits(x, "try-error"), output)
    output <- do.call(cbind, output)
    return(as.data.frame(output))
}

pez.divergence <- function(data, abundance.weighted=FALSE, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, which.eigen=1, quick=TRUE, q=0.0001)
{
    ## Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
        stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
    if(!is.null(traitgram)){
        if(length(traitgram) > 1){
            output <- vector("list", length(traitgram))
            for(i in seq_along(output))
                output[[i]] <- cbind(Recall(data, sqrt.phy, traitgram=traitgram[i], which.eigen=which.eigen, traitgram.p=traitgram.p, q=q), traitgram[i], sites(data))
            output <- do.call(rbind, output)
            names(output)[ncol(output)-1] <- "traitgram"
            names(output)[ncol(output)] <- "site"
            rownames(output) <- NULL
            return(output)
        } else {
            dist <- as.matrix(funct.phylo.dist(data, traitgram, traitgram.p))
            traitgram <- TRUE
        }      
    } else traitgram <- FALSE
    
    if(!is.null(ext.dist)){
        dist <- .check.ext.dist(ext.dist, species(data), ncol(data$comm))
        ext.dist <- TRUE
    } else ext.dist <- FALSE
    
    ## Setup
    if(!abundance.weighted)
        data$comm[data$comm > 0] <- 1
    coefs <- data.frame(row.names=rownames(data$comm))
    if(sqrt.phy)
        data <- .sqrt.phy(data)
    if(traitgram==FALSE & ext.dist==FALSE)
        dist <- cophenetic(data$phy)
    
    ## Filter metrics according to suitability and calculate
    functions <- setNames(c(.scheiner,.mpd,.mntd,.psv), c(".scheiner","mpd","mntd","psv"))
    if(abundance.weighted)
        functions <- append(functions, setNames(.pse,"pse"))
    if(sqrt.phy == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "colless", "gamma", "eed", "hed")]
    if(traitgram == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "pd", "colless", "gamma", "eed", "hed", "scheiner")]
    if(ext.dist == TRUE)
        functions <- functions[!names(functions) %in% c("pd", "psv", "psr", "colless", "gamma", "eed", "hed", "scheiner")]
    if(!is.binary.tree(data$phy) & "colless" %in% names(functions))
        warning("Cannot compute Colless' index with non-binary tree")
    output <- lapply(functions, function(x) try(x(data, dist=dist, abundance.weighted=abundance.weighted, which.eigen=which.eigen), silent=TRUE))
    
    ## Clean up output and return
    output <- Filter(function(x) !inherits(x, "try-error"), output)
    output <- do.call(cbind, output)
    return(as.data.frame(output))
}

pez.regularity <- function(data, abundance.weighted=FALSE, sqrt.phy=FALSE, traitgram=NULL, traitgram.p=2, ext.dist=NULL, which.eigen=1, quick=TRUE, q=0.0001)
{
    ## Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    if(sum(c(!is.null(traitgram), sqrt.phy, !is.null(ext.dist))) > 1)
        stop("Confusion now hath made its masterpiece!\nYou have specified more than one thing to do with a distance matrix.")
    if(!is.null(traitgram)){
        if(length(traitgram) > 1){
            output <- vector("list", length(traitgram))
            for(i in seq_along(output))
                output[[i]] <- cbind(Recall(data, sqrt.phy, traitgram=traitgram[i], which.eigen=which.eigen, traitgram.p=traitgram.p, q=q), traitgram[i], sites(data))
            output <- do.call(rbind, output)
            names(output)[ncol(output)-1] <- "traitgram"
            names(output)[ncol(output)] <- "site"
            rownames(output) <- NULL
            return(output)
        } else {
            dist <- as.matrix(funct.phylo.dist(data, traitgram, traitgram.p))
            traitgram <- TRUE
        }      
    } else traitgram <- FALSE
    
    if(!is.null(ext.dist)){
        dist <- .check.ext.dist(ext.dist, species(data), ncol(data$comm))
        ext.dist <- TRUE
    } else ext.dist <- FALSE
    
    ## Setup
    if(!abundance.weighted)
        data$comm[data$comm > 0] <- 1
    coefs <- data.frame(row.names=rownames(data$comm))
    if(sqrt.phy)
        data <- .sqrt.phy(data)
    if(traitgram==FALSE & ext.dist==FALSE)
        dist <- cophenetic(data$phy)
    
    ## Filter metrics according to suitability and calculate
    functions <- setNames(c(.colless,.gamma,.iac,.vpd,.taxon,.vntd,.eed,.hed,.haed), c(".colless","gamma","iac","vpd","taxon","vntd","eed","hed","haed"))    
    if(sqrt.phy == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "colless", "gamma", "eed", "hed")]
    if(traitgram == TRUE)
        functions <- functions[!names(functions) %in% c("psv", "psr", "pd", "colless", "gamma", "eed", "hed", "scheiner")]
    if(ext.dist == TRUE)
        functions <- functions[!names(functions) %in% c("pd", "psv", "psr", "colless", "gamma", "eed", "hed", "scheiner")]
    if(!is.binary.tree(data$phy) & "colless" %in% names(functions))
        warning("Cannot compute Colless' index with non-binary tree")
    output <- lapply(functions, function(x) try(x(data, dist=dist, abundance.weighted=abundance.weighted, which.eigen=which.eigen), silent=TRUE))
    
    ## Clean up output and return
    output <- Filter(function(x) !inherits(x, "try-error"), output)
    output <- do.call(cbind, output)
    return(as.data.frame(output))
}
