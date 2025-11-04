checkGesecaArgs <- function(E, pathways){
    # Error if pathways is not a list
    if (!is.list(pathways)) {
        stop("pathways should be a list with each element containing names of the stats argument")
    }

    # Error if E matrix doesn't have rownames
    if (is.null(rownames(E))) {
        stop("E rows should be named")
    }

    # Error if stats names are NA
    if (any(is.na(rownames(E)))) {
        stop("NAs in rownames(E) are not allowed")
    }

    # Error for duplicate gene names
    if (any(duplicated(rownames(E)))) {
        stop("Duplicate rownames(E) are not allowed")
    }

    # Error if E has non-finite values
    if (any(!is.finite(E))){
        stop("Not all E values are finite numbers")
    }
}

gesecaPreparePathways <- function(E, pathways, minSize, maxSize){
    res <- preparePathways(pathways, universe = rownames(E), minSize, maxSize)

    return(res)
}


# E should be row-centered
calcGesecaScores <- function(indxs, E){
    return(sum(colSums(E[indxs, , drop=FALSE])**2))
    # return(sum(abs(colSums(E[indxs, ]))))
}

# E should be row-centered
getGesecaLeadingEdge <- function(indxs, E, extend=FALSE, sort=TRUE) {
    profile <- colSums(E[indxs, , drop=FALSE])
    profile <- profile / sqrt(sum(profile**2))
    Ey <- E[indxs, , drop=FALSE]
    if (extend) {
        Ey <- E
    }
    weights <- (Ey %*% profile)[,1]
    weights <- weights / sqrt(apply(Ey**2, 1, sum))
    if (sort) {
        weights <- base::sort(weights, decreasing = TRUE)
    }
    return(names(weights[weights >= 0.5]))
}

# E should be row-centered
getGesecaLeadingEdge2 <- function(genes, E, extend=FALSE, sort=TRUE) {
    genes <- intersect(genes, rownames(E))
    profile <- colSums(E[genes,  , drop=FALSE])
    profile <- profile / sqrt(sum(profile**2))

    weights <- (E %*% profile)[,1]
    weights <- sort(weights, decreasing=TRUE)
    res <- calcGseaStat(stats = weights,
                 selectedStats = fmatch(genes, names(weights)),
                 returnLeadingEdge = TRUE, scoreType = "pos")$leadingEdge
    res <- names(weights)[res]
    return(res)
}
