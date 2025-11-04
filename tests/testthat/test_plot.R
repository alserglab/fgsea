context("Plots")

test_that("plotGseaTable works", {
    data(examplePathways)
    data(exampleRanks)
    fgseaRes <- fgsea(examplePathways, exampleRanks,
                      minSize=15, maxSize=100, eps=0.0)
    tf <- tempfile("plot", fileext = ".png")
    topPathways <- fgseaRes[head(order(pval), n=15)][order(NES), pathway]
    p <- plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, gseaParam=0.5)
    tf <- tempfile("plot", fileext = ".png")
    ggsave(tf, plot=p)
    expect_true(TRUE) # check that didn't fail before
})

test_that("plotEnrichment works", {
    data(examplePathways)
    data(exampleRanks)
    g <- plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
                        exampleRanks)
    tf <- tempfile("plot", fileext = ".png")
    ggsave(tf, plot=g)
    expect_true(TRUE) # check that didn't fail before
})

test_that("plotEnrichment and NAs", {
    data(examplePathways)
    data(exampleRanks)
    stats <- c(exampleRanks, "bla"=NA)
    pathway <- c(examplePathways[["5991130_Programmed_Cell_Death"]], "bla")
    expect_error(g <- plotEnrichment(pathway, stats), regexp = "finite")
})

test_that("plotEnrichment and zeroes", {
    data(examplePathways)
    data(exampleRanks)
    stats <- exampleRanks
    stats[seq(1, length(stats), 2)] <- 0
    pathway <- c(examplePathways[["5991130_Programmed_Cell_Death"]])
    pd <- plotEnrichmentData(pathway, stats, gseaParam = 0)
    fr <- fgseaSimple(list(p=pathway), stats=stats, nperm=10, gseaParam = 0)
    expect_equal(pd$posES, fr$ES)
})

test_that("plotGseaTable ignores empty pathways", {
    data(examplePathways)
    data(exampleRanks)
    fgseaRes <- fgsea(examplePathways, exampleRanks,
                      minSize=15, maxSize=100, eps=0.0)

    expect_equal(length(intersect(examplePathways[477], names(exampleRanks))), 0)
    p <- plotGseaTable(examplePathways[477], exampleRanks, fgseaRes, gseaParam=0.5)
    expect_true(is(p, "gg"))
})
