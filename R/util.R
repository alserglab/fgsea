preparePathways <- function(pathways, universe, minSize, maxSize){
    universeArg <- deparse(match.call()$universe)
    # Error if pathways is not a list
    if (!is.list(pathways)) {
        stop("pathways should be a list with each element containing gene identifiers")
    }

    # Error if stats is not named
    if (is.null(universe)) {
        .stopf("%s should not be null", universeArg)
    }

    # Error if stats names are NA
    if (any(is.na(universe))) {
        .stopf("NAs in %s are not allowed", universeArg)
    }

    # Error if stats names are empty string
    if (any(universe == "")) {
        .stopf("Empty strings are not allowed in %s", universeArg)
    }

    # Error for duplicate gene names
    if (any(duplicated(universe))) {
        .stopf("Duplicate values in %s not allowed", universeArg)
    }

    minSize <- max(minSize, 1)
    maxSize <- min(maxSize, length(universe)-1)

    pathwaysFiltered <- lapply(pathways, function(p) { unique(na.omit(fmatch(p, universe))) })
    pathwaysSizes <- lengths(pathwaysFiltered)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)

    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    list(filtered=pathwaysFiltered,
         sizes=pathwaysSizes)
}

.stopf <- function(fmt, ...) {
    stop(sprintf(fmt, ...))
}
