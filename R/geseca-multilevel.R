#' Runs multilevel Monte-Carlo variant for performing gene sets co-regulation analysis
#'
#' This function is based on the adaptive multilevel splitting Monte Carlo approach
#' and allows to estimate arbitrarily small P-values for the task of analyzing
#' variance along a set of genes.
#' @param pathways List of gene sets to check.
#' @param E expression matrix, rows corresponds to genes, columns corresponds to samples.
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param center a logical value indicating whether the gene expression should be centered to have zero mean before the analysis takes place.
#' The default is TRUE. The value is passed to \link[base]{scale}.
#' @param scale a logical value indicating whether the gene expression should be scaled to have
#' unit variance before the analysis takes place.
#' The default is FALSE. The value is passed to \link[base]{scale}.
#' @param sampleSize sample size for conditional sampling.
#' @param eps This parameter sets the boundary for calculating P-values.
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param BPPARAM Parallelization parameter used in bplapply.
#' @param nPermSimple Number of permutations in the simple geseca implementation
#' for preliminary estimation of P-values.
#'
#' @import BiocParallel
#' @import fastmatch
#' @import data.table
#' @return A table with GESECA results. Each row corresponds to a tested pathway.
#' The columns are the following
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
#' gr <- geseca(examplePathways, exampleExpressionMatrix, minSize=15, maxSize=500)
#' @export
geseca <- function(pathways,
                   E,
                   minSize     = 1,
                   maxSize     = nrow(E) - 1,
                   center      = TRUE,
                   scale       = FALSE,
                   sampleSize  = 101,
                   eps         = 1e-50,
                   nproc       = 0,
                   BPPARAM     = NULL,
                   nPermSimple = 1000)
{
    if (scale && any(apply(E, 1, sd) == 0)){
        stop("Cannot rescale a constant/zero gene expression rows to unit variance")
    }
    E <- t(base::scale(t(E), center=center, scale = scale))

    checkGesecaArgs(E, pathways)
    pp <- gesecaPreparePathways(E, pathways, minSize, maxSize)

    pathwayFiltered <- pp$filtered
    pathwaySizes <- pp$sizes
    pathwayNames <- names(pathwayFiltered)

    totalVar <- sum(rowSums(E**2))

    m <- length(pathwayFiltered)
    if (m == 0) {
        return(data.table(pathway = character(),
                          pctVar  = numeric(),
                          pval    = numeric(),
                          padj    = numeric(),
                          log2err = numeric(),
                          size    = integer()))
    }

    # Throw a warning message if sample size is less than 3
    if (sampleSize < 3){
        warning("sampleSize is too small, so sampleSize = 3 is set.")
        sampleSize <- max(3, sampleSize)
    }
    if (sampleSize %% 2 == 0){
        sampleSize <-  sampleSize + 1
    }

    eps <- max(0, min(1, eps))

    pathwayScores <- sapply(pathwayFiltered, calcGesecaScores, E = E)

    granularity <- max(1000, ceiling(nPermSimple / 128))
    permPerProc <- rep(granularity, floor(nPermSimple / granularity))
    if (nPermSimple - sum(permPerProc) > 0) {
        permPerProc <- c(permPerProc, nPermSimple - sum(permPerProc))
    }

    seeds <- sample.int(10^9, length(permPerProc))
    BPPARAM <- setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)


    grSimple <- gesecaSimpleImpl(pathwayScores, pathwaySizes,
                                 pathwayNames, permPerProc,
                                 seeds, m, E, BPPARAM)
    grSimple[,  gsScore := pathwayScores]

    roughEstimator <- log2((grSimple$nGeScore + 1) / (nPermSimple + 1))
    simpleError <- getSimpleError(roughEstimator, grSimple$nGeScore, nPermSimple)
    multilevelError <- sapply((grSimple$nGeScore + 1) / (nPermSimple + 1),
                              multilevelError, sampleSize)

    # todo: check this if
    if (all(multilevelError >= simpleError)){
        grSimple[, log2err := 1/log(2) * sqrt(trigamma(nGeScore + 1) -
                                                  trigamma((nPermSimple + 1)))]

        grSimple[, pctVar := gsScore / size / totalVar * 100]

        setorder(grSimple, pathway)
        grSimple <- grSimple[, .(pathway, pctVar, pval, padj,
                                 log2err,size)]
        grSimple <- grSimple[]
        return(grSimple)
    }

    dtGrSimple <- grSimple[multilevelError >= simpleError]
    dtGrSimple[, log2err := 1 / log(2) * sqrt(trigamma(nGeScore + 1) -
                                                  trigamma(nPermSimple + 1))]


    dtGrMultilevel <- grSimple[multilevelError < simpleError]
    mPathwaysList <- split(dtGrMultilevel, by = "size")

    # In most cases, this gives a speed increase with parallel launches.
    indxs <- sample(1:length(mPathwaysList))
    mPathwaysList <- mPathwaysList[indxs]





    # seed=sample.int(1e9, size=1)
    seeds <- sample.int(10^9, length(mPathwaysList))
    pvals <- bplapply(seq_along(mPathwaysList), function(i){
        x <- mPathwaysList[[i]]
        # scaledScore <- x[, pctVar] # this is pctVar
        scores <- x[["gsScore"]]
        size <- unique(x[["size"]])
        return(gesecaCpp(E           = E,
                         inpScores   = scores,
                         genesetSize = size,
                         sampleSize  = sampleSize,
                         seed        = seeds[i],
                         eps         = eps))
    }, BPPARAM = BPPARAM)

    pvals <- rbindlist(unlist(pvals, recursive = FALSE))

    result <- rbindlist(mPathwaysList)

    result[, pval := pvals$pval]
    result[, log2err := pvals$log2err]
    result[pval < eps, c("pval", "log2err") := list(eps, NA)]

    result[, padj := p.adjust(pval, method = "BH")]

    result <- rbindlist(list(result, dtGrSimple), use.names = TRUE)
    result[, pctVar := gsScore / size / totalVar * 100]
    result[, nGeScore := NULL]


    if (nrow(result[pval==eps & is.na(log2err)])){
        warning("For some pathways, in reality P-values are less than ",
                paste(eps),
                ". You can set the `eps` argument to zero for better estimation.")
    }
    result <- result[, .(pathway, pctVar, pval, padj,
                         log2err, size)]
    result <- result[order(pval)]
    return(result)
}

# This function finds error for rough P-value estimators.
# Function is based on the Clopper-Pearson interval.
getSimpleError <- function(roughEstimator, x, n, alpha = 0.025){
    leftBorder <- log2(qbeta(alpha,
                             shape1 = x,
                             shape2 = n - x + 1))
    rightBorder <- log2(qbeta(1 - alpha,
                              shape1 = x + 1,
                              shape2 = n - x))
    simpleError <- 0.5 * pmax(roughEstimator - leftBorder, rightBorder - roughEstimator)
    return(simpleError)
}

#' Collapse list of enriched pathways to independent ones (GESECA version, highly experimental).
#'
#' @param gesecaRes Table with results of running geseca(), should be filtered
#'                 by p-value, for example by selecting ones with padj < 0.01.
#' @param pathways List of pathways, should contain all the pathways present in
#'                 `gesecaRes`.
#' @param E expression matrix, the same as in `geseca()`.
#' @param center a logical value indicating whether the gene expression should be centered to have zero mean before the analysis takes place.
#' The default is TRUE. The value is passed to \link[base]{scale}.
#' @param scale a logical value indicating whether the gene expression should be scaled to have
#' unit variance before the analysis takes place.
#' The default is FALSE. The value is passed to \link[base]{scale}.
#' @param eps eps prameter for internal gesecaMultilevel runs. Default: min(c(1e-50, gesecaRes$pval))
#' @param checkDepth how much pathways to check against
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param BPPARAM Parallelization parameter used in bplapply.
collapsePathwaysGeseca <- function(gesecaRes,
                             pathways,
                             E,
                             center      = TRUE,
                             scale       = FALSE,
                             eps=min(c(1e-50, gesecaRes$pval)),
                             checkDepth=10,
                             nproc       = 0,
                             BPPARAM     = NULL) {
    universe <- rownames(E)
    E <- t(base::scale(t(E), center=center, scale = scale))

    pathways <- pathways[gesecaRes$pathway]
    pathways <- lapply(pathways, intersect, universe)

    pvalCondMat <- matrix(nrow=length(pathways),
                          ncol=length(pathways),
                          dimnames = list(names(pathways), names(pathways)))

    BPPARAM <- setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)

    message("Calculating pairwise simple p-values...")

    pvalCondRows <- BiocParallel::bplapply(seq_along(pathways), function(i) {

        p2n <- names(pathways)[i]
        p2 <- pathways[[p2n]]

        u1 <- setdiff(universe, p2)
        u2 <- p2

        pval1 <- setNames(rep(1, length(pathways)), names(pathways))
        pval2 <- setNames(rep(1, length(pathways)), names(pathways))

        gesecaRes1 <- gesecaSimple(pathways = pathways,
                                   E=E[u1, , drop=FALSE], center=FALSE, scale=FALSE,
                                   nperm = 100, BPPARAM = SerialParam())
        gesecaRes2 <- gesecaSimple(pathways = pathways,
                                   E=E[u2, , drop=FALSE], center=FALSE, scale=FALSE,
                                   nperm = 100, BPPARAM = SerialParam())

        pval1[gesecaRes1$pathway] <- gesecaRes1$pval
        pval2[gesecaRes2$pathway] <- gesecaRes2$pval

        pval1[pval1 < 1/100] <- 0
        pval2[pval2 < 1/100] <- 0

        # warnings about zero p-values
        suppressWarnings(pvalCond <- apply(cbind(pval1, pval2), 1, aggregation::fisher))

        pvalCond
    }, BPPARAM = BPPARAM)

    pvalCondMat <- do.call(rbind, pvalCondRows)
    rownames(pvalCondMat) <- names(pathways)
    diag(pvalCondMat) <- NA
    pvalCondMax <- apply(pvalCondMat, 2, max, na.rm=TRUE)
    parentPathway <- setNames(names(pathways)[apply(pvalCondMat, 2, which.max)],
                              names(pathways))


    checkDepth <- 10
    message("Calculating pairwise multilevel p-values...")
    condTable <- rbindlist(bplapply(which(pvalCondMax == 0), function(i) {

        p1n <- names(pathways)[i]
        p1 <- pathways[[i]]

        pathwaysToCheck <- names(pathways)[which(pvalCondMat[, p1n] == 0)]

        distances <- sapply(pathways[pathwaysToCheck], function(p2) { length(intersect(p1, p2))**2/
                length(p1) / length(p2) } )

        # distances <- distances[-i]
        distances <- sort(distances, decreasing = TRUE)
        pathwaysToCheck <- names(head(distances, checkDepth))

        maxPvalCond <- 0
        parentPathway <- NULL

        for (p2n in pathwaysToCheck) {
            p2 <- pathways[[p2n]]

            u1 <- setdiff(universe, p2)
            u2 <- p2

            chisq <- qchisq(maxPvalCond, df=2*2, lower.tail = FALSE)
            pvalBound <- exp(-chisq / 2)
            # if at least one of gesecaRes1$pval, gesecaRes2$pval is below pvalBound,
            # then pvalCond is less than maxPvalCond
            eps1 <- max(pvalBound/2, 1e-50)

            suppressWarnings({ # warnings about reaching eps
                gesecaRes1 <- geseca(pathways = list(p=p1),
                                     E=E[u1, , drop=FALSE], center=FALSE, scale=FALSE,
                                     eps=eps1, BPPARAM = SerialParam())
                gesecaRes2 <- geseca(pathways = list(p=p1),
                                     E=E[u2, , drop=FALSE], center=FALSE, scale=FALSE,
                                     eps=eps1, BPPARAM = SerialParam())
            })

            pvals <- c(min(gesecaRes1$pval, 1, na.rm=TRUE),
                       min(gesecaRes2$pval, 1, na.rm=TRUE))
            pvalCond <- aggregation::fisher(pvals)
            if (length(pvals) == 0) {
                pvalCond <- 1
            }

            if (pvalCond > maxPvalCond) {
                # if (match(p2n, pathwaysToCheck) >= 5) {
                #     message(match(p2n, pathwaysToCheck), ": ", p1n, " ##### ", p2n)
                #     message("      ", maxPvalCond, " -> ", pvalCond)
                # }

                maxPvalCond <- pvalCond
                parentPathway <- p2n
            }
        }
        return(data.table(pathway=p1n, parentPathway=parentPathway, pvalCond=maxPvalCond))
    }, BPPARAM=BPPARAM))

    pvalCondMax[condTable$pathway] <- condTable$pvalCond
    parentPathway[condTable$pathway] <- condTable$parentPathway

    pvalCondReciprocal <- pvalCondMat[cbind(parentPathway, names(parentPathway))]
    names(pvalCondReciprocal) <- names(parentPathway)

    message("Calculating reciprocal multilevel p-values...")
    reciprocalCondTable <- rbindlist(bplapply(which(pvalCondReciprocal == 0), function(i) {
        p2n <- names(pathways)[i]
        p2 <- pathways[[i]]

        p1n <- parentPathway[p2n]
        p1 <- pathways[[p1n]]

        u1 <- setdiff(universe, p2)
        u2 <- p2

        suppressWarnings({ # warnings about reaching eps
            gesecaRes1 <- geseca(pathways = list(p=p1),
                                 E=E[u1, , drop=FALSE], center=FALSE, scale=FALSE,
                                 eps=eps, BPPARAM = SerialParam())
            gesecaRes2 <- geseca(pathways = list(p=p1),
                                 E=E[u2, , drop=FALSE], center=FALSE, scale=FALSE,
                                 eps=eps, BPPARAM = SerialParam())
        })

        pvals <- c(min(gesecaRes1$pval, 1, na.rm=TRUE),
                   min(gesecaRes2$pval, 1, na.rm=TRUE))

        pvalCond <- aggregation::fisher(pvals)
        if (length(pvals) == 0) {
            pvalCond <- 1
        }

        return(data.table(pathway=p1n, parentPathway=p2n, pvalCond=pvalCond))
    }, BPPARAM=BPPARAM))

    pvalCondReciprocal[reciprocalCondTable$parentPathway] <- reciprocalCondTable$pvalCond

    res <- copy(gesecaRes)

    res[, pvalCond := pvalCondMax[pathway]]
    res[, parentPathway := parentPathway[pathway]]
    res[, reciprocalPvalCond := pvalCondReciprocal[pathway]]

    res[, pScore := exp(log(pvalCond) + (log(pval)-log(pvalCond))*(log(pvalCond)/(log(pvalCond)+log(pvalCondReciprocal))))]
    # res[, pScore := sqrt(pval * pvalCond)]

    res <- res[]

    return(res)
}
