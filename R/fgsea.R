#' Wrapper to run methods for preranked gene set enrichment analysis.
#'
#' This function provide an interface to two existing functions:
#' \link[fgsea]{fgseaSimple}, \link[fgsea]{fgseaMultilevel}.
#' By default, the \link[fgsea]{fgseaMultilevel} function is used for analysis.
#' For compatibility with the previous implementation you can pass the `nperm` argument to the function.
#' @param pathways List of gene sets to check.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param gseaParam GSEA parameter value, all gene-level statis are raised to the power of `gseaParam`
#' @param ... optional arguments for functions \link[fgsea]{fgseaSimple}, \link[fgsea]{fgseaMultilevel}
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#' @export
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' fgseaRes <- fgsea(examplePathways, exampleRanks, maxSize=500)
#' # Testing only one pathway is implemented in a more efficient manner
#' fgseaRes1 <- fgsea(examplePathways[1], exampleRanks)
fgsea <- function(pathways, stats, minSize = 1, maxSize = length(stats)-1, gseaParam = 1, ...){
    arguments <- list(...)
    if ("nperm" %in% names(arguments)){
        warning("You are trying to run fgseaSimple. ",
                "It is recommended to use fgseaMultilevel. ",
                "To run fgseaMultilevel, you need to remove ",
                "the nperm argument in the fgsea function call.")
        res <- fgseaSimple(pathways = pathways, stats = stats,
                           minSize = minSize, maxSize = maxSize,
                           gseaParam = gseaParam, ...)
    }
    else{
        res <- fgseaMultilevel(pathways = pathways, stats = stats,
                               minSize = minSize, maxSize = maxSize,
                               gseaParam = gseaParam, ...)
    }
    res
}

prepareStats <- function(stats, scoreType="std", gseaParam=1) {
    # Error if stats are non-finite
    if (any(!is.finite(stats))){
        stop("Not all stats values are finite numbers")
    }

    # Warning message for ties in stats
    ties <- sum(duplicated(stats[stats != 0]))
    if (ties != 0) {
        warning("There are ties in the preranked stats (",
                paste(round(ties * 100 / length(stats), digits = 2)),
                "% of the list).\n",
                "The order of those tied genes will be arbitrary, which may produce unexpected results.")
    }

    if (all(stats > 0) & scoreType == "std"){
        warning("All values in the stats vector are greater than zero and scoreType is \"std\", ",
                "maybe you should switch to scoreType = \"pos\".")
    }

    stats <- sort(stats, decreasing=TRUE)
    stats <- abs(stats) ^ gseaParam

    # scaling weights to integer values for more stable calculations
    # sum(stats) < 2**30 avoids integer overflows
    scaleCoeff <- 2**30/sum(stats)
    if (scaleCoeff >= 1) {
        # scaling by an integer coefficient avoids information loss for integer input stats
        scaleCoeff <- floor(scaleCoeff)
    }

    # TODO: special treatment for small non-zero values?
    stats <- round(stats * scaleCoeff)
    storage.mode(stats) <- "integer"
    return(stats)
}

preparePathwaysAndStats <- function(pathways, stats, minSize, maxSize, gseaParam, scoreType){
    stats <- prepareStats(stats, scoreType, gseaParam)
    res <- preparePathways(pathways, universe=names(stats), minSize, maxSize)

    res$stats <- stats

    res
}



#' Calculates GSEA statistics for a given query gene set
#'
#' Takes \emph{O(k log k)} time, where \emph{k} is a size of `selectedSize`.
#' @param stats Named numeric vector with gene-level statistics
#'  sorted in decreasing order (order is not checked).
#' @param selectedStats Indexes of selected genes in the `stats` array.
#' @param gseaParam GSEA weight parameter (0 is unweighted, suggested value is 1).
#' @param returnAllExtremes If TRUE return not only the most extreme point, but all of them. Can be used for enrichment plot
#' @param returnLeadingEdge If TRUE return also leading edge genes.
#' @param scoreType This parameter defines the GSEA score type. Possible options are ("std", "pos", "neg")
#' @return Value of GSEA statistic if both returnAllExtremes and returnLeadingEdge are FALSE.
#' Otherwise returns list with the folowing elements:
#' \itemize{
#' \item res -- value of GSEA statistic
#' \item tops -- vector of top peak values of cumulative enrichment statistic for each gene;
#' \item bottoms -- vector of bottom peak values of cumulative enrichment statistic for each gene;
#' \item leadingGene -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#' @export
#' @examples
#' data(exampleRanks)
#' data(examplePathways)
#' ranks <- sort(exampleRanks, decreasing=TRUE)
#' es <- calcGseaStat(ranks, na.omit(match(examplePathways[[1]], names(ranks))))
calcGseaStat <- function(stats,
                         selectedStats,
                         gseaParam=1,
                         returnAllExtremes = FALSE,
                         returnLeadingEdge = FALSE,
                         scoreType         = c("std", "pos", "neg")) {
    scoreType <- match.arg(scoreType)
    S <- selectedStats
    r <- stats
    p <- gseaParam

    S <- sort(S)
    stopifnot(all(head(S, -1) < tail(S, -1)))

    m <- length(S)
    N <- length(r)
    if (m == N) {
        stop("GSEA statistic is not defined when all genes are selected")
    }
    NR <- (sum(abs(r[S])^p))
    rAdj <- abs(r[S])^p
    if (NR == 0) {
        # this is equivalent to rAdj being rep(eps, m)
        rCumSum <- seq_along(rAdj) / length(rAdj)
    } else {
        rCumSum <- cumsum(rAdj) / NR
    }


    tops <- rCumSum - (S - seq_along(S)) / (N - m)
    if (NR == 0) {
        # this is equivalent to rAdj being rep(eps, m)
        bottoms <- tops - 1 / m
    } else {
        bottoms <- tops - rAdj / NR
    }

    maxP <- max(tops)
    minP <- min(bottoms)

    switch(scoreType,
           std = geneSetStatistic <- ifelse(maxP == -minP, 0, ifelse(maxP > -minP, maxP, minP)),
           pos = geneSetStatistic <- maxP,
           neg = geneSetStatistic <- minP)

    if (!returnAllExtremes && !returnLeadingEdge) {
        return(geneSetStatistic)
    }

    res <- list(res=geneSetStatistic)
    if (returnAllExtremes) {
        res <- c(res, list(tops=tops, bottoms=bottoms))
    }
    if (returnLeadingEdge) {
        switch(scoreType,
               std = leadingEdge <- if (maxP > -minP) {
                   S[seq_along(S) <= which.max(tops)]
               } else if (maxP < -minP) {
                   rev(S[seq_along(S) >= which.min(bottoms)])
               } else {
                   NULL
               },
               pos = leadingEdge <- S[seq_along(S) <= which.max(tops)],
               neg = leadingEdge <- rev(S[seq_along(S) >= which.min(bottoms)]))

        res <- c(res, list(leadingEdge=leadingEdge))
    }
    res
}

#' Runs preranked gene set enrichment analysis.
#'
#' The function takes about \emph{O(nk^\{3/2\})} time,
#' where \emph{n} is number of permutations and \emph{k} is a maximal
#' size of the pathways. That means that setting `maxSize` parameter with a value of ~500
#' is strongly recommended.
#' @param pathways List of gene sets to check.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param scoreType This parameter defines the GSEA score type.
#' Possible options are ("std", "pos", "neg").
#' By default ("std") the enrichment score is computed as in the original GSEA.
#' The "pos" and "neg" score types are intended to be used for one-tailed tests
#' (i.e. when one is interested only in positive ("pos") or negateive ("neg") enrichment).
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param gseaParam GSEA parameter value, all gene-level statis are raised to the power of `gseaParam`
#'                  before calculation of GSEA enrichment scores.
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value;
#'  \item padj -- a BH-adjusted p-value;
#'  \item ES -- enrichment score, same as in Broad GSEA implementation;
#'  \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#'  \item nMoreExtreme` -- a number of times a random gene set had a more
#'      extreme enrichment score value;
#'  \item size -- size of the pathway after removing genes not present in `names(stats)`.
#'  \item leadingEdge -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#'
#' @export
#' @import data.table
#' @import BiocParallel
#' @import fastmatch
#' @import stats
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' fgseaRes <- fgseaSimple(examplePathways, exampleRanks, nperm=10000, maxSize=500)
#' # Testing only one pathway is implemented in a more efficient manner
#' fgseaRes1 <- fgseaSimple(examplePathways[1], exampleRanks, nperm=10000)
fgseaSimple <- function(pathways,
                        stats,
                        nperm,
                        minSize   = 1,
                        maxSize   = length(stats)-1,
                        scoreType = c("std", "pos", "neg"),
                        nproc     = 0,
                        gseaParam = 1,
                        BPPARAM   = NULL) {
    scoreType <- match.arg(scoreType)
    pp <- preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, scoreType)
    pathwaysFiltered <- pp$filtered
    pathwaysSizes <- pp$sizes
    stats <- pp$stats


    m <- length(pathwaysFiltered)

    if (m == 0) {
        return(data.table(pathway=character(),
                          pval=numeric(),
                          padj=numeric(),
                          ES=numeric(),
                          NES=numeric(),
                          nMoreExtreme=numeric(),
                          size=integer(),
                          leadingEdge=list()))
    }


    granularity <- max(1000, ceiling(nperm / 128)) # not having more than 128 threads
    permPerProc <- rep(granularity, floor(nperm / granularity))
    if (nperm - sum(permPerProc) > 0) {
        permPerProc <- c(permPerProc, nperm - sum(permPerProc))
    }

    pval=nLeZero=nGeZero=leZeroMean=geZeroMean=nLeEs=nGeEs=NULL

    seeds <- sample.int(10^9, length(permPerProc))

    BPPARAM <- setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)

    gseaStatRes <- do.call(rbind,
                lapply(pathwaysFiltered, calcGseaStat,
                       stats=stats,
                       returnLeadingEdge=TRUE,
                       scoreType=scoreType))


    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatRes[, "res"])




    pvals <- fgseaSimpleImpl(pathwayScores, pathwaysSizes,
                             pathwaysFiltered, leadingEdges,
                             permPerProc, seeds, m, stats,
                             BPPARAM, scoreType)
    if (nrow(pvals[is.na(pval)]) > 0){
        warning("There were ",
                paste(nrow(pvals[is.na(pval)])),
                " pathways for which P-values were not calculated properly due to ",
                "unbalanced gene-level statistic values")
    }

    pvals[, nLeZero := NULL]
    pvals[, nGeZero := NULL]
    pvals[, leZeroMean := NULL]
    pvals[, geZeroMean := NULL]
    pvals[, nLeEs := NULL]
    pvals[, nGeEs := NULL]

    setcolorder(pvals, c("pathway", "pval", "padj", "ES", "NES",
                         "nMoreExtreme", "size", "leadingEdge"))
    # Makes pvals object printable immediatly
    pvals <- pvals[]

    pvals
}

calcGseaStatBatch <- function(stats, selectedStats, geneRanks=seq_along(stats),
                               gseaParam=1) {
    stats <- abs(stats)^gseaParam
    calcGseaStatBatchCpp(stats, selectedStats, geneRanks)
}

#' Runs label-permuring gene set enrichment analysis.
#'
#' @param pathways List of gene sets to check.
#' @param mat Gene expression matrix. Row name should be the same as in 'pathways'
#' @param labels Numeric vector of labels for the correlation score of the same length as the number
#'               of columns in `mat`
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param gseaParam GSEA parameter value, all gene-level statis are raised to the power of `gseaParam`
#'                  before calculation of GSEA enrichment scores.
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value;
#'  \item padj -- a BH-adjusted p-value;
#'  \item ES -- enrichment score, same as in Broad GSEA implementation;
#'  \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#'  \item nMoreExtreme` -- a number of times a random gene set had a more
#'      extreme enrichment score value;
#'  \item size -- size of the pathway after removing genes not present in `names(stats)`.
#'  \item leadingEdge -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#'
#' @export
#' @examples
#' \donttest{
#' library(limma)
#' library(GEOquery)
#' es <- getGEO("GSE19429", AnnotGPL = TRUE)[[1]]
#' exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")
#' es <- es[!grepl("///", fData(es)$`Gene ID`), ]
#' es <- es[fData(es)$`Gene ID` != "", ]
#' es <- es[order(apply(exprs(es), 1, mean), decreasing=TRUE), ]
#' es <- es[!duplicated(fData(es)$`Gene ID`), ]
#' rownames(es) <- fData(es)$`Gene ID`
#'
#' pathways <- reactomePathways(rownames(es))
#' mat <- exprs(es)
#' labels <- as.numeric(as.factor(gsub(" .*", "", es$title)))
#' fgseaRes <- fgseaLabel(pathways, mat, labels, nperm = 1000, minSize = 15, maxSize = 500)
#' }
#'
#' @importFrom Matrix invPerm
fgseaLabel <- function(pathways, mat, labels, nperm,
                      minSize=1, maxSize=nrow(mat)-1,
                      nproc=0,
                      gseaParam=1,
                      BPPARAM=NULL) {

    granularity <- 100
    permPerProc <- rep(granularity, floor(nperm / granularity))
    if (nperm - sum(permPerProc) > 0) {
        permPerProc <- c(permPerProc, nperm - sum(permPerProc))
    }
    seeds <- sample.int(10^9, length(permPerProc))

    BPPARAM <- setUpBPPARAM(nproc, BPPARAM)

    tmatSc <- scale(t(mat))
    labelsSc <- scale(labels)[, 1]

    minSize <- max(minSize, 1)

    pathwaysFiltered <- lapply(pathways, function(p) { unique(na.omit(fmatch(p, rownames(mat)))) })
    pathwaysSizes <- sapply(pathwaysFiltered, length)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
    m <- length(toKeep)

    if (m == 0) {
        return(data.table(pathway=character(),
                          pval=numeric(),
                          padj=numeric(),
                          ES=numeric(),
                          NES=numeric(),
                          nMoreExtreme=numeric(),
                          size=integer(),
                          leadingEdge=list()))
    }

    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    corRanks <- var(tmatSc, labelsSc)[,1]
    ranksOrder <- order(corRanks, decreasing=TRUE)
    ranksOrderInv <- invPerm(ranksOrder)
    stats <- corRanks[ranksOrder]

    pathwaysReordered <- lapply(pathwaysFiltered, function(x) ranksOrderInv[x])

    gseaStatRes <- do.call(rbind,
                           lapply(pathwaysReordered, calcGseaStat,
                                  stats=stats,
                                  returnLeadingEdge=TRUE))


    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatRes[, "res"])

    universe <- seq_along(stats)

    counts <- bplapply(seq_along(permPerProc), function(i) {
        nperm1 <- permPerProc[i]

        labelPerms <- do.call(cbind, replicate(nperm1, sample(labelsSc), simplify = FALSE))
        randCorRanks <- var(tmatSc, labelPerms)

        randEsPs <- lapply(seq_len(nperm1), function(i) {
            randCorRanks1 <- randCorRanks[, i]
            ranksOrder <- sort.list(randCorRanks1, decreasing=TRUE)
            geneRanks <- invPerm(ranksOrder)
            stats <- randCorRanks1[ranksOrder]

            randEsP <- calcGseaStatBatch(
                stats = stats,
                selectedStats = pathwaysFiltered,
                geneRanks = geneRanks,
                gseaParam = gseaParam)
        })

        randEsPs <- do.call(cbind, randEsPs)

        leEs <- apply(sweep(randEsPs, MARGIN = 1, pathwayScores, `<=`), 1, sum)
        geEs <- apply(sweep(randEsPs, MARGIN = 1, pathwayScores, `>=`), 1, sum)
        leZero <- apply(randEsPs <= 0, 1, sum)
        geZero <- apply(randEsPs >= 0, 1, sum)
        leZeroSum <- apply(pmin(randEsPs, 0), 1, sum)
        geZeroSum <- apply(pmax(randEsPs, 0), 1, sum)

        data.table(pathway=seq_len(m),
                   leEs=leEs, geEs=geEs,
                   leZero=leZero, geZero=geZero,
                   leZeroSum=leZeroSum, geZeroSum=geZeroSum
        )
    }, BPPARAM=BPPARAM)


    counts <- rbindlist(counts)

    # Getting rid of check NOTEs
    leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
    pathway=padj=pval=ES=NES=geZeroMean=leZeroMean=NULL
    nMoreExtreme=nGeEs=nLeEs=size=NULL
    leadingEdge=NULL
    .="damn notes"


    pvals <- counts[,
                    list(pval=min((1+sum(leEs)) / (1 + sum(leZero)),
                                  (1+sum(geEs)) / (1 + sum(geZero))),
                         leZeroMean = sum(leZeroSum) / sum(leZero),
                         geZeroMean = sum(geZeroSum) / sum(geZero),
                         nLeEs=sum(leEs),
                         nGeEs=sum(geEs)
                    )
                    ,
                    by=.(pathway)]
    pvals[, padj := p.adjust(pval, method="BH")]
    pvals[, ES := pathwayScores[pathway]]
    pvals[, NES := ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean))]
    pvals[, leZeroMean := NULL]
    pvals[, geZeroMean := NULL]

    pvals[, nMoreExtreme :=  ifelse(ES > 0, nGeEs, nLeEs)]
    pvals[, nLeEs := NULL]
    pvals[, nGeEs := NULL]

    pvals[, size := pathwaysSizes[pathway]]
    pvals[, pathway := names(pathwaysFiltered)[pathway]]

    pvals[, leadingEdge := .(leadingEdges)]


    # Makes pvals object printable immediatly
    pvals <- pvals[]

    pvals
}

#' Collapse list of enriched pathways to independent ones.
#'
#' @param fgseaRes Table with results of running fgsea(), should be filtered
#'                 by p-value, for example by selecting ones with padj < 0.01.
#' @param pathways List of pathways, should contain all the pathways present in
#'                 `fgseaRes`.
#' @param stats Gene-level statistic values used for ranking, the same as
#'              in `fgsea()`.
#' @param pval.threshold Two pathways are considered dependent when p-value
#'                       of enrichment of one pathways on background of another
#'                       is greater then `pval.threshold`.
#' @param nperm Number of permutations to test for independence, should be
#'              several times greater than `1/pval.threhold`.
#'              Default value: `10/pval.threshold`.
#' @param gseaParam GSEA parameter, same as for `fgsea()`
#' @return Named list with two elments: `mainPathways` containing IDs of pathways
#'         not reducable to each other, and `parentPathways` with vector describing
#'         for all the pathways to which ones they can be reduced. For
#'         pathways from `mainPathwyas` vector `parentPathways` contains `NA` values.
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=10000, maxSize=500)
#' collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
#'                                       examplePathways, exampleRanks)
#' mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
#'                          order(-NES), pathway]
#' @export
collapsePathways <- function(fgseaRes,
                             pathways,
                             stats,
                             pval.threshold=0.05,
                             nperm=10/pval.threshold,
                             gseaParam=1) {
    universe <- names(stats)

    pathways <- pathways[fgseaRes$pathway]
    pathways <- lapply(pathways, intersect, universe)

    parentPathways <- setNames(rep(NA, length(pathways)), names(pathways))

    for (i in seq_along(pathways)) {
        p <- names(pathways)[i]
        if (!is.na(parentPathways[p])) {
            next
        }

        pathwaysToCheck <- setdiff(names(which(is.na(parentPathways))), p)
        pathwaysUp <- fgseaRes[pathway %fin% pathwaysToCheck & ES >= 0][, pathway]
        pathwaysDown <- fgseaRes[pathway %fin% pathwaysToCheck & ES < 0][, pathway]

        if (length(pathwaysToCheck) == 0) {
            break
        }

        minPval <- setNames(rep(1, length(pathwaysToCheck)), pathwaysToCheck)

        u1 <- setdiff(universe, pathways[[p]])

        fgseaResUp1 <- fgseaSimple(pathways = pathways[pathwaysUp], stats=stats[u1],
                                   nperm=nperm, maxSize=length(u1)-1, BPPARAM = SerialParam(),
                                   gseaParam=gseaParam, scoreType = "pos")
        fgseaResDown1 <- fgseaSimple(pathways = pathways[pathwaysDown], stats=stats[u1],
                                     nperm=nperm, maxSize=length(u1)-1, BPPARAM = SerialParam(),
                                     gseaParam=gseaParam, scoreType = "neg")
        fgseaRes1 <- rbindlist(list(fgseaResUp1, fgseaResDown1), use.names = TRUE)

        minPval[fgseaRes1$pathway] <- pmin(minPval[fgseaRes1$pathway], fgseaRes1$pval)

        u2 <- pathways[[p]]

        fgseaResUp2 <- fgseaSimple(pathways = pathways[pathwaysUp], stats=stats[u2],
                                   nperm=nperm, maxSize=length(u2)-1, BPPARAM = SerialParam(),
                                   gseaParam=gseaParam, scoreType = "pos")
        fgseaResDown2 <- fgseaSimple(pathways = pathways[pathwaysDown], stats=stats[u2],
                                     nperm=nperm, maxSize=length(u2)-1, BPPARAM = SerialParam(),
                                     gseaParam=gseaParam, scoreType = "neg")
        fgseaRes2 <- rbindlist(list(fgseaResUp2, fgseaResDown2), use.names = TRUE)

        minPval[fgseaRes2$pathway] <- pmin(minPval[fgseaRes2$pathway], fgseaRes2$pval)

        parentPathways[names(which(minPval > pval.threshold))] <- p
    }

    return(list(mainPathways=names(which(is.na(parentPathways))),
                parentPathways=parentPathways))
}

#' Runs preranked gene set enrichment analysis for preprocessed input data.
#'
#' @param pathwayScores Vector with enrichment scores for the `pathways`.
#' @param pathwaysSizes Vector of pathways sizes.
#' @param pathwaysFiltered Filtered pathways.
#' @param leadingEdges Leading edge genes.
#' @param permPerProc  Parallelization parameter for permutations.
#' @param seeds Seed vector
#' @param toKeepLength  Number of `pathways` that meet the condition for `minSize` and `maxSize`.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param BPPARAM Parallelization parameter used in bplapply.
#' @param scoreType This parameter defines the GSEA score type. Possible options are ("std", "pos", "neg")
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value;
#'  \item padj -- a BH-adjusted p-value;
#'  \item ES -- enrichment score, same as in Broad GSEA implementation;
#'  \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#'  \item nMoreExtreme` -- a number of times a random gene set had a more
#'      extreme enrichment score value;
#'  \item size -- size of the pathway after removing genes not present in `names(stats)`.
#'  \item leadingEdge -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
fgseaSimpleImpl <- function(pathwayScores, pathwaysSizes,
                            pathwaysFiltered, leadingEdges,
                            permPerProc, seeds,toKeepLength,
                            stats, BPPARAM, scoreType){
    K <- max(pathwaysSizes)
    universe <- seq_along(stats)

    counts <- bplapply(seq_along(permPerProc), function(i) {
        nperm1 <- permPerProc[i]
        leEs <- rep(0, toKeepLength)
        geEs <- rep(0, toKeepLength)
        leZero <- rep(0, toKeepLength)
        geZero <- rep(0, toKeepLength)
        leZeroSum <- rep(0, toKeepLength)
        geZeroSum <- rep(0, toKeepLength)
        if (toKeepLength == 1) {
            set.seed(seeds[i])
            for (j in seq_len(nperm1)) {
                randSample <- sample.int(length(universe), K)
                randEsP <- calcGseaStat(
                    stats = stats,
                    selectedStats = randSample,
                    gseaParam = 1,
                    scoreType = scoreType)
                leEs <- leEs + (randEsP <= pathwayScores)
                geEs <- geEs + (randEsP >= pathwayScores)
                leZero <- leZero + (randEsP <= 0)
                geZero <- geZero + (randEsP >= 0)
                leZeroSum <- leZeroSum + pmin(randEsP, 0)
                geZeroSum <- geZeroSum + pmax(randEsP, 0)
            }
        } else {
            aux <- calcGseaStatCumulativeBatch(
                stats = stats,
                gseaParam = 1,
                pathwayScores = pathwayScores,
                pathwaysSizes = pathwaysSizes,
                iterations = nperm1,
                seed = seeds[i],
                scoreType = scoreType)
            leEs = get("leEs", aux)
            geEs = get("geEs", aux)
            leZero = get("leZero", aux)
            geZero = get("geZero", aux)
            leZeroSum = get("leZeroSum", aux)
            geZeroSum = get("geZeroSum", aux)
        }
        data.table(pathway=seq_len(toKeepLength),
                   leEs=leEs, geEs=geEs,
                   leZero=leZero, geZero=geZero,
                   leZeroSum=leZeroSum, geZeroSum=geZeroSum
                   )
    }, BPPARAM=BPPARAM)

    counts <- rbindlist(counts)

    # Getting rid of check NOTEs
    leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
    pathway=padj=pval=ES=NES=geZeroMean=leZeroMean=NULL
    nMoreExtreme=nGeEs=nLeEs=nLeZero=nGeZero=size=NULL
    leadingEdge=NULL
    .="damn notes"



    pvals <- counts[, list(leZeroMean = sum(leZeroSum) / sum(leZero),
                           geZeroMean = sum(geZeroSum) / sum(geZero),
                           nLeZero = sum(leZero),
                           nGeZero = sum(geZero),
                           nLeEs = sum(leEs),
                           nGeEs = sum(geEs)),
                    by = .(pathway)]

    pvals[, ES := pathwayScores[pathway]]

    pvals[, NES := as.numeric(NA)]

    switch(scoreType,
           std = pvals[(ES > 0 & geZeroMean != 0) | (ES <= 0 & leZeroMean != 0),
                           NES := ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean))],
           pos = pvals[(ES >= 0 & geZeroMean != 0), NES := ES / geZeroMean],
           neg = pvals[(ES <= 0 & leZeroMean != 0), NES := ES / abs(leZeroMean)])

    pvals[, pval := as.numeric(NA)]
    pvals[!is.na(NES), pval := pmin((1+nLeEs) / (1 + nLeZero),
                        (1+nGeEs) / (1 + nGeZero))]


    pvals[, padj := as.numeric(NA)]
    pvals[!is.na(pval), padj := p.adjust(pval, method = "BH")]

    switch(scoreType,
           std = pvals[, nMoreExtreme :=  ifelse(ES > 0, nGeEs, nLeEs)],
           pos = pvals[, nMoreExtreme :=  nGeEs],
           neg = pvals[, nMoreExtreme :=  nLeEs])

    pvals[, size := pathwaysSizes[pathway]]
    pvals[, pathway := names(pathwaysFiltered)[pathway]]
    pvals[, leadingEdge := .(leadingEdges)]

    pvals
}

setUpBPPARAM <- function(nproc=0, BPPARAM=NULL){
    if (is.null(BPPARAM)) {
        if (nproc != 0) {
            if (nproc == 1) {
                result <- SerialParam(progressbar = TRUE)
            }
            if (.Platform$OS.type == "windows") {
                # windows doesn't support multicore, using snow instead
                result <- SnowParam(workers = nproc, progressbar = TRUE)
            } else {
                result <- MulticoreParam(workers = nproc, progress = TRUE)
            }
        } else {
            result <- bpparam()
        }
        return(result)
    }
    else {
        return(BPPARAM)
    }
}




