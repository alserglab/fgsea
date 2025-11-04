#' Runs simple variant for performing gene sets co-regulation analysis
#'
#' This function is based on the rude Monte Carlo sampling approach
#' and P-value calculation accuracy is limited to `1 / nperm` value.
#' @param pathways List of gene sets to check.
#' @param E expression matrix, rows corresponds to genes, columns corresponds to samples.
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param center a logical value indicating whether the gene expression should be centered to have zero mean before the analysis takes place.
#' The default is TRUE. The value is passed to \link[base]{scale}.
#' @param scale a logical value indicating whether the gene expression should be scaled to have unit variance before the analysis takes place.
#' The default is FALSE. The value is passed to \link[base]{scale}.
#' @param nperm Number of permutations to do. Minimal possible nominal p-value is about 1/nperm
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param BPPARAM Parallelization parameter used in bplapply.
#'
#' @import BiocParallel
#' @import fastmatch
#' @import data.table
#' @return A table with GESECA results. Each row corresponds to a tested pathway. The columns are the following
#' \itemize{
#' \item pathway -- name of the pathway as in `names(pathways)`;
#' \item pctVar -- percent of explained variance along gene set;
#' \item pval -- P-value that corresponds to the gene set score;
#' \item padj -- a BH-adjusted p-value;
#' \item size -- size of the pathway after removing genes not present in `rownames(E)`.
#' }
#'
#' @examples
#' data("exampleExpressionMatrix")
#' data("examplePathways")
#' gesecaRes <- gesecaSimple(examplePathways, exampleExpressionMatrix, minSize=15, maxSize=500)
#' @export
gesecaSimple <- function(pathways,
                         E,
                         minSize    = 1,
                         maxSize    = nrow(E) - 1,
                         center     = TRUE,
                         scale      = FALSE,
                         nperm      = 1000,
                         nproc      = 0,
                         BPPARAM    = NULL){
    if (scale && any(apply(E, 1, sd) == 0)){
        stop("Cannot rescale the constant/zero gene expression rows to unit variance")
    }
    E <- t(base::scale(t(E), center=center, scale = scale))

    checkGesecaArgs(E, pathways)
    pp <- gesecaPreparePathways(E, pathways, minSize, maxSize)
    pathwayFiltered <- pp$filtered
    pathwaySizes <- pp$sizes
    pathwayNames <- names(pathwayFiltered)

    m <- length(pathwayFiltered)
    if (m == 0) {
        return(data.table(pathway = character(),
                          pctVar  = numeric(),
                          pval    = numeric(),
                          padj    = numeric(),
                          size    = integer()))
    }

    pathwayScores <- sapply(pathwayFiltered, calcGesecaScores, E = E)

    granularity <- max(1000, ceiling(nperm / 128))
    permPerProc <- rep(granularity, floor(nperm / granularity))
    if (nperm - sum(permPerProc) > 0) {
        permPerProc <- c(permPerProc, nperm - sum(permPerProc))
    }

    seeds <- sample.int(10^9, length(permPerProc))

    BPPARAM <- setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)

    pvals <- gesecaSimpleImpl(pathwayScores, pathwaySizes,
                              pathwayNames, permPerProc,
                              seeds, m, E, BPPARAM)
    pvals[, pctVar := pathwayScores]

    totalVar <- sum(rowSums(E**2))
    pvals[, pctVar := pctVar / size / totalVar * 100]

    setnames(pvals, "nGeScore", "nMoreExtreme")

    setcolorder(pvals, c("pathway", "pctVar", "pval",
                         "padj", "nMoreExtreme", "size"))
    pvals <- pvals[order(pval)]
    return(pvals)
}



gesecaSimpleImpl <- function(pathwayScores,
                             pathwaySizes,
                             pathwayNames,
                             permPerProc,
                             seeds,
                             toKeepLength,
                             E,
                             BPPARAM){
    K <- max(pathwaySizes)
    universe <- seq_len(nrow(E))

    counts <- bplapply(seq_along(permPerProc), function(i){
        set.seed(seeds[i])
        nperm1 <- permPerProc[i]
        scores <- gesecaCumScores(nperm1, pathwaySizes, E)
        geScore <- t(apply(scores, 1, function(x){
            return(x >= pathwayScores)
        }))
        geScore <- colSums(geScore)
        data.table(pathway=seq_len(toKeepLength),
                   geScore=geScore)
    }, BPPARAM = BPPARAM)
    counts <- rbindlist(counts)

    pathway=geScore=NULL
    .="damn notes"
    pvals <- counts[, list(nGeScore = sum(geScore)), by = .(pathway)]

    pvals[, pval := (nGeScore + 1) / (sum(permPerProc) + 1)]
    pvals[, padj := p.adjust(pval, method = "BH")]

    pvals[, size := pathwaySizes[pathway]]
    pvals[, pathway := pathwayNames[pathway]]

    return(pvals)
}


gesecaCumScores <- function(nperm, pathwaySizes, E){
    K <- max(pathwaySizes)
    universe <- seq_len(nrow(E))
    uPathwaySizes <- unique(pathwaySizes)

    scores <- lapply(seq_len(nperm), function(i){
        randIndxs <- sample.int(nrow(E), K)

        Esubset <- E[randIndxs, , drop=FALSE]

        Ecum <- do.call(cbind, apply(Esubset, 2, cumsum, simplify = FALSE))
        Esq <- (Ecum ^ 2)

        res <- unname(rowSums(Esq)[pathwaySizes])
        return(res)
    })
    scores <- do.call(rbind, scores)
    colnames(scores) <- pathwaySizes
    return(scores)
}

