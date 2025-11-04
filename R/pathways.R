#' Returns a list of Reactome pathways for given Entrez gene IDs
#' @param genes Entrez IDs of query genes.
#' @return A list of vectors with gene sets.
#' @export
#' @examples
#' data(exampleRanks)
#' pathways <- reactomePathways(names(exampleRanks))
reactomePathways <- function(genes) {
    stopifnot(requireNamespace("reactome.db"))
    stopifnot(requireNamespace("AnnotationDbi"))
    pathways <- na.omit(AnnotationDbi::select(reactome.db::reactome.db,
                                              keys=genes,
                                              c("PATHID"),
                                              keytype = 'ENTREZID'))

    pathways <- split(pathways$ENTREZID, pathways$PATHID)

    pathway2name <- as.data.table(AnnotationDbi::select(
        reactome.db::reactome.db,
        names(pathways),
        c("PATHNAME"),
        'PATHID'))

    # workaround for reactome.db having one-to-many PATHID to PATHNAM mapping
    PATHID=NULL
    pathway2name <- pathway2name[!duplicated(PATHID)]



    # Hack to get rid of check NOTEs
    PATHNAME=NULL

    # Remove organism prefix
    pathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]

    # workaround for reactome.db having duplicated names for pathways
    # example: Immune system (5991156 and 6096039)
    name2pathways <- split(pathway2name$PATHID, pathway2name$PATHNAME)

    pathways <- lapply(name2pathways, function(x) unique(do.call(c, pathways[x])))

    pathways <- pathways[!is.na(names(pathways))]

    pathways
}

#' Returns a list of pathways from a GMT file.
#' @param gmt.file Path to a GMT file.
#' @return A list of vectors with gene sets.
#' @examples
#' pathways <- gmtPathways(system.file(
#'      "extdata", "mouse.reactome.gmt", package="fgsea"))
#' @importFrom  utils head tail
#' @export
gmtPathways <- function(gmt.file) {
    pathwayLines <- strsplit(readLines(gmt.file), "\t")
    pathways <- lapply(pathwayLines, tail, -2)
    names(pathways) <- sapply(pathwayLines, head, 1)
    pathways
}

#' Write collection of pathways (list of vectors) to a gmt file
#' @param pathways a named list of vectors with gene ids
#' @param gmt.file name of the output file
#' @import data.table
#' @examples
#' data(examplePathways)
#' writeGmtPathways(examplePathways, tempfile("examplePathways", fileext=".gmt"))
#' @export
writeGmtPathways <- function(pathways, gmt.file) {
    dt <- data.table(pathway=names(pathways), "NA", genes=unlist(lapply(pathways, paste, collapse="\t")))
    fwrite(dt, file=gmt.file, sep="\t", col.names = FALSE, quote = FALSE)
}

#' Effeciently converts collection of pathways using AnnotationDbi::mapIds function. Parameters
#' are the sames as for mapIds except for keys, which is assumed to be a list of vectors.
#' @param x the AnnotationDb object. But in practice this will mean an object derived from an AnnotationDb object such as a OrgDb or ChipDb object.
#' @param keys a list of vectors with gene ids
#' @param column the column to search on
#' @param keytype the keytype that matches the keys used
#' @param ... other parameters passed to AnnotationDbi::mapIds
#' @seealso AnnotationDbi::mapIds
#' @examples
#' library(org.Mm.eg.db)
#' data(exampleRanks)
#' fgseaRes <- fgsea(examplePathways, exampleRanks, maxSize=500, eps=1e-4)
#' fgseaRes[, leadingEdge := mapIdsList(org.Mm.eg.db, keys=leadingEdge, column="SYMBOL", keytype="ENTREZID")]
#' @export
mapIdsList <- function(x, keys, column, keytype, ...) {
    keysFlat <- unlist(keys)
    keysUnique <- unique(keysFlat)
    ansUnique <- AnnotationDbi::mapIds(x=x, keys=keysUnique, column=column, keytype=keytype, ...)
    ansFlat <- ansUnique[keysFlat]
    ans <- split(ansFlat, rep(seq_along(keys), lengths(keys)))
    names(ans) <- names(keys)
    ans
}
