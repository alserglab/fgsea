#' Runs preranked gene set enrichment analysis.
#'
#' This feature is based on the adaptive multilevel splitting Monte Carlo approach.
#' This allows us to exceed the results of simple sampling and calculate arbitrarily small P-values.
#' @param pathways List of gene sets to check.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param sampleSize The size of a random set of genes which in turn has size = pathwaySize
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param eps This parameter sets the boundary for calculating the p value.
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
#' @param nPermSimple Number of permutations in the simple fgsea implementation
#' for preliminary estimation of P-values.
#' @param absEps deprecated, use `eps` parameter instead
#'@export
#' @import fastmatch
#' @import data.table
#' @import BiocParallel
#' @return A table with GSEA results. Each row corresponds to a tested pathway. The columns are the following
#' \itemize{
#' \item pathway -- name of the pathway as in `names(pathway)`;
#' \item pval -- an enrichment p-value;
#' \item padj -- a BH-adjusted p-value;
#' \item log2err -- the expected error for the standard deviation of the P-value logarithm.
#' \item ES -- enrichment score, same as in Broad GSEA implementation;
#' \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#' \item size -- size of the pathway after removing genes not present in `names(stats)`.
#' \item leadingEdge -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' fgseaMultilevelRes <- fgseaMultilevel(examplePathways, exampleRanks, maxSize=500)
fgseaMultilevel <- function(pathways,
                            stats,
                            sampleSize  = 101,
                            minSize     = 1,
                            maxSize     = length(stats)-1,
                            eps         = 1e-50,
                            scoreType   = c("std", "pos", "neg"),
                            nproc       = 0,
                            gseaParam   = 1,
                            BPPARAM     = NULL,
                            nPermSimple = 1000,
                            absEps      = NULL)
{
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
                          log2err=numeric(),
                          ES=numeric(),
                          NES=numeric(),
                          size=integer(),
                          leadingEdge=list()))
    }



    # Warning message for deprecated absEps parameter
    if (!is.null(absEps)){
        warning("You are using deprecated argument `absEps`. ",
                "Use `eps` argument instead. ",
                "`absEps` was assigned to `eps`.")
        eps <-  absEps
    }

    # Warning message for to small value for sampleSize
    if (sampleSize < 3){
        warning("sampleSize is too small, so sampleSize = 3 is set.")
        sampleSize <- max(3, sampleSize)
    }

    #To avoid warnings during the check
    log2err=nMoreExtreme=pathway=pval=padj=NULL
    nLeZero=nGeZero=leZeroMean=geZeroMean=nLeEs=nGeEs=isCpGeHalf=NULL
    ES=NES=size=leadingEdge=NULL
    .="damn notes"

    minSize <- max(minSize, 1)
    eps <- max(0, min(1, eps))


    if (sampleSize %% 2 == 0){
        sampleSize <-  sampleSize + 1
    }


    gseaStatRes <- do.call(rbind,
                           lapply(pathwaysFiltered,
                                  calcGseaStat,
                                  stats             = stats,
                                  returnLeadingEdge = TRUE,
                                  scoreType         = scoreType))

    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatRes[, "res"])


    seeds <- sample.int(10^9, 1)
    BPPARAM <- setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)

    simpleFgseaRes <- fgseaSimpleImpl(pathwayScores=pathwayScores, pathwaysSizes=pathwaysSizes,
                                      pathwaysFiltered=pathwaysFiltered, leadingEdges=leadingEdges,
                                      permPerProc=nPermSimple, seeds=seeds, toKeepLength=m,
                                      stats=stats, BPPARAM=SerialParam(), scoreType=scoreType)

    switch(scoreType,
           std = simpleFgseaRes[, modeFraction := ifelse(ES >= 0, nGeZero, nLeZero)],
           pos = simpleFgseaRes[, modeFraction := nGeZero],
           neg = simpleFgseaRes[, modeFraction := nLeZero])


    simpleFgseaRes[, leZeroMean := NULL]
    simpleFgseaRes[, geZeroMean := NULL]
    simpleFgseaRes[, nLeEs := NULL]
    simpleFgseaRes[, nGeEs := NULL]
    simpleFgseaRes[, nLeZero := NULL]
    simpleFgseaRes[, nGeZero := NULL]

    simpleFgseaRes[modeFraction < 10, pval := as.numeric(NA)]
    simpleFgseaRes[modeFraction < 10, padj := as.numeric(NA)]
    simpleFgseaRes[modeFraction < 10, NES := as.numeric(NA)]

    if (any(simpleFgseaRes$modeFraction < 10)){
        warning("There were ",
                paste(sum(simpleFgseaRes$modeFraction < 10)),
                " pathways for which P-values were not calculated properly due to ",
                "unbalanced (positive and negative) gene-level statistic values. ",
                "For such pathways pval, padj, NES, log2err are set to NA. ",
                "You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = ",
                paste0(format(nPermSimple * 10, scientific = FALSE), ")"))
    }

    # Storing NA fgseaSimple results in a separate data.table
    naSimpleRes <- simpleFgseaRes[is.na(pval)]
    naSimpleRes[, padj := as.numeric(NA)]
    naSimpleRes[, log2err := as.numeric(NA)]
    naSimpleRes[, modeFraction := NULL]

    simpleFgseaRes <- simpleFgseaRes[!is.na(pval)]

    leftBorder <- log2(qbeta(0.025,
                             shape1 = simpleFgseaRes$nMoreExtreme,
                             shape2=nPermSimple - simpleFgseaRes$nMoreExtreme + 1))

    rightBorder <- log2(qbeta(1 - 0.025,
                         shape1 = simpleFgseaRes$nMoreExtreme + 1,
                         shape2=nPermSimple - simpleFgseaRes$nMoreExtreme))

    crudeEstimator <- log2((simpleFgseaRes$nMoreExtreme + 1) / (nPermSimple + 1))

    simpleError <- 0.5 * pmax(crudeEstimator - leftBorder, rightBorder - crudeEstimator)


    multError <- sapply((simpleFgseaRes$nMoreExtreme + 1) / (nPermSimple + 1), multilevelError, sampleSize)


    if (all(multError >= simpleError)){
        simpleFgseaRes[, log2err := 1/log(2)*sqrt(trigamma(nMoreExtreme + 1) - trigamma((nPermSimple + 1)))]
        simpleFgseaRes[, modeFraction := NULL]
        simpleFgseaRes <- rbindlist(list(simpleFgseaRes, naSimpleRes), use.names = TRUE)

        setorder(simpleFgseaRes, pathway)
        simpleFgseaRes[, "nMoreExtreme" := NULL]
        setcolorder(simpleFgseaRes, c("pathway", "pval", "padj", "log2err",
                                      "ES", "NES", "size", "leadingEdge"))

        simpleFgseaRes <- simpleFgseaRes[]
        return(simpleFgseaRes)
    }


    dtSimpleFgsea <- simpleFgseaRes[multError >= simpleError]
    dtSimpleFgsea[, log2err := 1/log(2)*sqrt(trigamma(nMoreExtreme + 1) - trigamma(nPermSimple + 1))]
    dtSimpleFgsea[, modeFraction := NULL]



    dtMultilevel <- simpleFgseaRes[multError < simpleError]
    # Probability estimation in the denominator of the multilevel algortihm
    # (s_r(q) >= 0 in the notation of the article):
    dtMultilevel[, "denomProb" := (modeFraction + 1) / (nPermSimple + 1)]

    multilevelPathwaysList <- split(dtMultilevel, by="size")
    # In most cases, this gives a speed increase with parallel launches.
    indxs <- sample(1:length(multilevelPathwaysList))
    multilevelPathwaysList <- multilevelPathwaysList[indxs]

    seed=sample.int(1e9, size=1)

    sign <- if (scoreType %in% c("pos", "neg")) TRUE else FALSE
    cpp.res <- multilevelImpl(multilevelPathwaysList, stats,
                              sampleSize, seed, eps, sign=sign,
                              BPPARAM=BPPARAM)
    cpp.res <- rbindlist(cpp.res)


    result <- rbindlist(multilevelPathwaysList)

    # `cppMpval` - P-values that are computed in cpp code
    # `isCpGeHalf` is a flag that mathces: whether the conditional probability
    # is greater than or equal to 0.5 (see article for details)
    result[, pval := pmin(1, cpp.res$cppMPval / denomProb)]
    result[, isCpGeHalf := cpp.res$cppIsCpGeHalf]
    result[, log2err := multilevelError(pval, sampleSize = sampleSize)]
    result[isCpGeHalf == FALSE, log2err:= NA]

    if (!all(result$isCpGeHalf)){
        warning("For some of the pathways the P-values were likely overestimated. ",
                "For such pathways log2err is set to NA.")
    }

    result[, isCpGeHalf := NULL]
    result[, modeFraction := NULL]
    result[, denomProb := NULL]


    result <- rbindlist(list(result, dtSimpleFgsea, naSimpleRes), use.names = TRUE)
    result[, nMoreExtreme := NULL]

    result[pval < eps, c("pval", "log2err") := list(eps, NA)]
    result[, padj := p.adjust(pval, method = "BH")]

    if (nrow(result[pval==eps & is.na(log2err)])){
        warning("For some pathways, in reality P-values are less than ",
                paste(eps),
                ". You can set the `eps` argument to zero for better estimation.")
    }

    setcolorder(result, c("pathway", "pval", "padj", "log2err",
                          "ES", "NES", "size", "leadingEdge"))

    setorder(result, pathway)

    result <- result[]
    result
}

#'Calculates the expected error for the standard deviation of the P-value logarithm.
#'
#' @param pval P-value
#' @param sampleSize equivavlent to sampleSize in fgseaMultilevel
#'@export
#'@return The value of the expected error
#' @examples
#' expectedError <- multilevelError(pval=1e-10, sampleSize=1001)
multilevelError <- function(pval, sampleSize){
    return(sqrt(floor(-log2(pval) + 1) * (trigamma((sampleSize+1)/2) - trigamma(sampleSize+1)))/log(2))
}

#' Calculates P-values for preprocessed data.
#' @param multilevelPathwaysList List of pathways for which P-values will be calculated.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param sampleSize The size of a random set of genes which in turn has size = pathwaySize
#' @param seed `seed` parameter from `fgseaMultilevel`
#' @param eps This parameter sets the boundary for calculating the p value.
#' @param sign This option will be used in future implementations.
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @return List of P-values.
multilevelImpl <- function(multilevelPathwaysList,
                           stats,
                           sampleSize,
                           seed,
                           eps,
                           sign=FALSE,
                           BPPARAM=NULL){
    #To avoid warnings during the check
    size=ES=NULL
    res <- bplapply(multilevelPathwaysList, function(x){
        eps <- eps * min(x$denomProb)
        return(fgseaMultilevelCpp(x[, ES], stats, unique(x[, size]),
                                  sampleSize, seed, eps, sign))
    }, BPPARAM = BPPARAM)
    return(res)
}
