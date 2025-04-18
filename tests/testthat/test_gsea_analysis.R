context("GSEA analysis")

test_that("fgseaSimple works", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    fgseaRes <- fgseaSimple(examplePathways, exampleRanks, nperm=nperm, maxSize=500)


    expect_equal(fgseaRes[23, ES], 0.5788464)
    expect_equal(fgseaRes[23, nMoreExtreme], 0)
    expect_gt(fgseaRes[1237, nMoreExtreme], 50 * nperm / 1000)

    expect_true("70385" %in% fgseaRes[grep("5991851", pathway), leadingEdge][[1]])
    expect_true(!"68549" %in% fgseaRes[grep("5991851", pathway), leadingEdge][[1]])

    expect_true(!"11909" %in% fgseaRes[grep("5992314", pathway), leadingEdge][[1]])
    expect_true("69386" %in% fgseaRes[grep("5992314", pathway), leadingEdge][[1]])

    # analyzing one pathway is done in a different way
    fgsea1Res <- fgseaSimple(examplePathways[1237], exampleRanks, nperm=nperm, maxSize=500)
    expect_gt(fgseaRes[1, nMoreExtreme], 50 * nperm / 1000)

    # specifying number of threads
    fgseaRes <- fgseaSimple(examplePathways, exampleRanks, nperm=2000, maxSize=100, nproc=2)

    # all nMoreExtreme being even is a sign of invalid parallelization
    expect_false(all(fgseaRes$nMoreExtreme %% 2 == 0))
})

test_that("fgseaSimple is reproducible independent of bpparam settings", {

    data(examplePathways)
    data(exampleRanks)
    nperm <- 2000

    set.seed(42)
    fr <- fgseaSimple(examplePathways[1:2], exampleRanks, nperm=nperm, maxSize=500, nproc=1)


    set.seed(42)
    fr1 <- fgseaSimple(examplePathways[1:2], exampleRanks, nperm=nperm, maxSize=500)
    expect_equal(fr1$nMoreExtreme, fr$nMoreExtreme)

    set.seed(42)
    fr2 <- fgseaSimple(examplePathways[1:2], exampleRanks, nperm=nperm, maxSize=500, nproc=0)
    expect_equal(fr2$nMoreExtreme, fr$nMoreExtreme)

    set.seed(42)
    fr3 <- fgseaSimple(examplePathways[1:2], exampleRanks, nperm=nperm, maxSize=500, nproc=2)
    expect_equal(fr3$nMoreExtreme, fr$nMoreExtreme)
})

test_that("fgseaSimple works with zero pathways", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    fgseaRes <- fgseaSimple(examplePathways, exampleRanks, nperm=nperm,
                      minSize=50, maxSize=10)
    expect_equal(nrow(fgseaRes), 0)
    fgseaRes1 <- fgseaSimple(examplePathways[1], exampleRanks, nperm=nperm)
    expect_equal(colnames(fgseaRes), colnames(fgseaRes1))
})

test_that("fgseaLabel works", {
    mat <- matrix(rnorm(1000*20),1000,20)
    labels <- rep(1:2, c(10, 10))

    rownames(mat) <- as.character(seq_len(nrow(mat)))
    pathways <- list(sample(rownames(mat), 100),
                     sample(rownames(mat), 200))

    fgseaRes <- fgseaLabel(pathways, mat, labels, nperm = 1000, minSize = 15, maxSize = 500)
    expect_true(!is.null(fgseaRes))
})

test_that("fgseaLabel example works", {
    skip_on_bioc()

    # hack to run the example with devtools
    helpFile <- system.file("../man/fgseaLabel.Rd", package="fgsea")

    if (helpFile == "") {
        example("fgseaLabel", run.donttest = TRUE, local = FALSE)
    } else {
        tcon <- textConnection("exampleCode", open="w", local=TRUE)
        tools::Rd2ex(helpFile, out=tcon)
        close(tcon)

        econ <- textConnection(exampleCode, open="r")
        source(econ)
    }

    expect_true(!is.null(fgseaRes))
})

test_that("Ties detection in ranking works", {
    data(examplePathways)
    data(exampleRanks)
    exampleRanks.ties <- exampleRanks
    exampleRanks.ties[41] <- exampleRanks.ties[42]
    exampleRanks.ties.zero <- exampleRanks.ties
    exampleRanks.ties.zero[41] <- exampleRanks.ties.zero[42] <- 0

    expect_silent(fgseaSimple(examplePathways, exampleRanks, nperm=100,
                              minSize=10, maxSize=50, BPPARAM=SerialParam()))

    expect_warning(fgseaSimple(examplePathways, exampleRanks.ties, nperm=100,
                               minSize=10, maxSize=50, BPPARAM=SerialParam()))

    expect_silent(fgseaSimple(examplePathways, exampleRanks.ties.zero, nperm=100,
                              minSize=10, maxSize=50, BPPARAM=SerialParam()))
})

test_that("fgseaSimple correctly checks gene names", {
    data(examplePathways)
    data(exampleRanks)
    exampleRanks.dupNames <- exampleRanks
    names(exampleRanks.dupNames)[41] <- names(exampleRanks.dupNames)[42]

    expect_error(fgseaSimple(examplePathways, exampleRanks.dupNames, nperm=100, minSize=10, maxSize=50, nproc=1))

    ranks <- exampleRanks
    names(ranks)[41] <- NA
    expect_error(fgseaSimple(examplePathways, ranks, nperm=100, minSize=10, maxSize=50, nproc=1))

    ranks <- exampleRanks
    names(ranks)[41] <- ""
    expect_error(fgseaSimple(examplePathways, ranks, nperm=100, minSize=10, maxSize=50, nproc=1))

    ranks <- unname(exampleRanks)
    expect_error(fgseaSimple(examplePathways, ranks, nperm=100, minSize=10, maxSize=50, nproc=1))

})

test_that("fgseaSimple returns leading edge ordered by decreasing of absolute statistic value", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    fgseaRes <- fgseaSimple(examplePathways, exampleRanks, nperm=nperm, maxSize=50)

    expect_true(abs(exampleRanks[fgseaRes$leadingEdge[[1]][1]]) >
                abs(exampleRanks[fgseaRes$leadingEdge[[1]][2]]))

    expect_true(abs(exampleRanks[fgseaRes$leadingEdge[[2]][1]]) >
                abs(exampleRanks[fgseaRes$leadingEdge[[2]][2]]))
})

test_that("collapsePathways work", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    pathways <- list(p1=examplePathways$`5991851_Mitotic_Prometaphase`)
    pathways <- c(pathways, list(p2=unique(c(pathways$p1, sample(names(exampleRanks), 20)))))
    pathways <- c(pathways, list(p3=sample(pathways$p1, floor(length(pathways$p1) * 0.8))))
    fgseaRes <- fgseaSimple(pathways, exampleRanks, nperm=nperm, maxSize=500)
    collapsedPathways <- collapsePathways(fgseaRes[order(pval)],
                                          pathways, exampleRanks)
    collapsedPathways$mainPathways
    expect_identical("p1", collapsedPathways$mainPathways)
})


test_that("fgseaSimple throws a warning when there are unbalanced gene-level statistic values", {
    data(exampleRanks)
    data(examplePathways)

    ranks <- sort(exampleRanks, decreasing = TRUE)
    firstNegative <- which(ranks < 0)[1]
    ranks <- c(abs(ranks[1:(firstNegative - 1)]) ^ 0.1, ranks[firstNegative:length(ranks)])

    pathway <- list(testPathway = names(ranks)[1:100])
    set.seed(1)
    expect_warning(fgseaSimple(pathway, ranks, nperm = 200, minSize = 15, maxSize = 500))
})

test_that("fgseaSimple and fgseaMultilevel properly handle duplicated genes in gene sets", {
    data(exampleRanks)
    data(examplePathways)


    pathways <- list(p1=examplePathways$`5991851_Mitotic_Prometaphase`,
                     p2=rep(examplePathways$`5991851_Mitotic_Prometaphase`, 2))

    set.seed(42)
    fr1 <- fgseaSimple(pathways, exampleRanks, nperm=1000)
    expect_equal(fr1$size[1], fr1$size[2])

    fr2 <- suppressWarnings(fgseaMultilevel(pathways, exampleRanks, eps = 1e-4))
    expect_equal(fr2$size[1], fr2$size[2])
})


test_that("fgsea works correctly if gene sets have signed ES = 0 and scoreType != std", {
    data("exampleRanks")

    stats <- exampleRanks
    stats <- sort(stats, decreasing=TRUE)

    # positive scoreType with gs: ES+(gs) = 0
    gsInTail <- list("gsInTail" = names(stats)[(length(stats) - 14) : length(stats)])

    set.seed(1)
    expect_silent(fr <- fgseaSimple(gsInTail, stats, nperm = 1000, scoreType = "pos"))

    # negative scoreType with gs: ES-(gs) = 0
    gsInHead <- list("gsInHead" = names(stats)[1:15])

    set.seed(1)
    expect_silent(fr <- fgseaSimple(gsInHead, stats, nperm = 1000, scoreType = "neg"))
})

test_that("fgsea skips pathways containing all the possible genes", {
    data("exampleRanks")
    fr <- fgseaSimple(list(p=names(exampleRanks)), exampleRanks, nperm = 1, maxSize=Inf)
    expect_true(!is.null(fr))
})

test_that("fgseaMultilevel handels superdiscrete cases (like issue #151)", {
    set.seed(42)
    stats <- rep(1, 5000)
    names(stats) <- paste0("g", seq_along(stats))
    system.time(res <- fgseaMultilevel(pathways=list(p=names(stats)[1:10]),
                           stats=stats, scoreType = "pos", eps=0, sampleSize = 21))
    expect_true(is.na(res$log2err))
})

test_that("leadingEdge interacts correctly with scoreType", {
    data(examplePathways)
    data(exampleRanks)
    set.seed(42)
    nperm <- 100
    fr <- fgseaSimple(examplePathways, exampleRanks, nperm=nperm, minSize=15, maxSize=50)
    fr <- fr[order(-abs(NES))]

    posP <- head(fr[ES > 0, pathway], 1)
    negP <- head(fr[ES < 0, pathway], 1)

    pp <- examplePathways[c(posP, negP)]

    frStd <- fgseaSimple(pp, exampleRanks, nperm=nperm, scoreType = "std")
    expect_true(exampleRanks[frStd$leadingEdge[[1]][1] ] > 0)
    expect_true(exampleRanks[frStd$leadingEdge[[2]][1] ] < 0)


    frPos <- fgseaSimple(pp, exampleRanks, nperm=nperm, scoreType = "pos")
    expect_true(exampleRanks[frPos$leadingEdge[[1]][1] ] > 0)
    expect_true(exampleRanks[frPos$leadingEdge[[2]][1] ] > 0)

    frNeg <- fgseaSimple(pp, exampleRanks, nperm=nperm, scoreType = "neg")
    expect_true(exampleRanks[frNeg$leadingEdge[[1]][1] ] < 0)
    expect_true(exampleRanks[frNeg$leadingEdge[[2]][1] ] < 0)
})

test_that("fgseaSimple is reproducible between platforms, issues #170, #80", {
    feats <- c(3, -1, -4.1, 42, 0, 12, 13, -13, 0.01, 0)
    names(feats) <- paste0("gene", 1:length(feats))

    some_sets <- list(
        "pathway1" = c("gene2", "gene3"),
        "pathway4" = paste0("gene", 1:7)
    )
    set.seed(42)
    fr <- fgseaSimple(some_sets, feats, nproc = 1, nperm=1000)
    expect_identical(fr$nMoreExtreme, c(123, 80))
})
