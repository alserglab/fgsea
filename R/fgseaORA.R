#' Simple overrepresentation analysis based on hypergeometric test
#' @param pathways List of gene sets to check.
#' @param genes Set of query genes
#' @param universe A universe from whiche `genes` were selected
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @return A table with ORA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value from hypergeometric test;
#'  \item padj -- a BH-adjusted p-value;
#'  \item foldEnrichment -- degree of enrichment relative to background;
#'  \item overlap -- size of the overlap;
#'  \item size -- size of the gene set;
#'  \item leadingEdge -- vector with overlapping genes.
#' }
#'
#' @export
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' foraRes <- fora(examplePathways, genes=tail(names(exampleRanks), 200), universe=names(exampleRanks))
fora <- function(pathways, genes, universe, minSize=1, maxSize=length(universe)-1) {
    pp <- preparePathways(pathways, universe, minSize, maxSize)
    pathwaysFiltered <- pp$filtered
    pathwaysSizes <- pp$sizes


    if (length(pathwaysFiltered) == 0){
        return(data.table(pathway=character(),
                          pval=numeric(),
                          padj=numeric(),
                          foldEnrichment=numeric(),
                          overlap=integer(),
                          size=integer(),
                          overlapGenes=list()))
    }


    if (!all(genes %in% universe)) {
        warning("Not all of the input genes belong to the universe, such genes were removed")
    }

    genesFiltered <- unique(na.omit(fmatch(genes, universe)))

    overlaps <- lapply(pathwaysFiltered, intersect, genesFiltered)

    overlapGenes <- lapply(overlaps, function(x) universe[x])

    overlapsT <- data.table(
        q=sapply(overlaps, length),
        m=sapply(pathwaysFiltered, length),
        n=length(universe)-sapply(pathwaysFiltered, length),
        k=length(genesFiltered))

    overlapsT[, es := (q/k)/(m/length(universe))]

    # q-1 because we want probability of having >=q white balls
    pathways.pvals <- with(overlapsT,
                           phyper(q-1, m, n, k, lower.tail = FALSE))

    res <- data.table(pathway=names(pathwaysFiltered),
                      pval=pathways.pvals,
                      padj=p.adjust(pathways.pvals, method="BH"),
                      foldEnrichment=overlapsT$es,
                      overlap=overlapsT$q,
                      size=overlapsT$m,
                      overlapGenes=overlapGenes)

    res <- res[order(pval),]
    res
}

#' Collapse list of enriched pathways to independent ones. Version for ORA hypergeometric test.
#'
#' @param foraRes Table with results of running fgsea(), should be filtered
#'                 by p-value, for example by selecting ones with padj < 0.01.
#' @param pathways List of pathways, should contain all the pathways present in
#'                 `fgseaRes`.
#' @param genes Set of query genes, same as in `fora()`
#' @param universe A universe from whiche `genes` were selected, same as in `fora()`
#' @param pval.threshold Two pathways are considered dependent when p-value
#'                       of enrichment of one pathways on background of another
#'                       is greater then `pval.threshold`.
#' @return Named list with two elments: `mainPathways` containing IDs of pathways
#'         not reducable to each other, and `parentPathways` with vector describing
#'         for all the pathways to which ones they can be reduced. For
#'         pathways from `mainPathwyas` vector `parentPathways` contains `NA` values.
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' foraRes <- fora(examplePathways, genes=tail(names(exampleRanks), 200), universe=names(exampleRanks))
#' collapsedPathways <- collapsePathwaysORA(foraRes[order(pval)][padj < 0.01],
#'                                              examplePathways,
#'                                              genes=tail(names(exampleRanks), 200),
#'                                              universe=names(exampleRanks))
#'
#' mainPathways <- foraRes[pathway %in% collapsedPathways$mainPathways][
#'                           order(pval), pathway]
#' @export
collapsePathwaysORA <- function(foraRes,
                             pathways,
                             genes,
                             universe,
                             pval.threshold=0.05) {

    pathways <- pathways[foraRes$pathway]
    pathways <- lapply(pathways, intersect, universe)

    parentPathways <- setNames(rep(NA, length(pathways)), names(pathways))

    for (i in seq_along(pathways)) {
        p <- names(pathways)[i]
        if (!is.na(parentPathways[p])) {
            next
        }

        pathwaysToCheck <- setdiff(names(which(is.na(parentPathways))), p)

        if (length(pathwaysToCheck) == 0) {
            break
        }

        minPval <- setNames(rep(1, length(pathwaysToCheck)), pathwaysToCheck)

        u1 <- setdiff(universe, pathways[[p]])
        foraRes1 <- fora(pathways = pathways[pathwaysToCheck],
                         genes=intersect(genes, u1),
                         universe=u1,
                         maxSize=length(u1)-1)
        minPval[foraRes1$pathway] <- pmin(minPval[foraRes1$pathway], foraRes1$pval)

        u2 <- pathways[[p]]
        foraRes2 <- fora(pathways = pathways[pathwaysToCheck],
                                genes=intersect(genes, u2),
                                universe=u2,
                                maxSize=length(u2)-1)
        minPval[foraRes2$pathway] <- pmin(minPval[foraRes2$pathway], foraRes2$pval)

        parentPathways[names(which(minPval > pval.threshold))] <- p
    }

    return(list(mainPathways=names(which(is.na(parentPathways))),
                parentPathways=parentPathways))
}
