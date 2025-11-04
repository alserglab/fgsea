context("GSEA stat calculation")

test_that("calcGseaStats works", {
    stats <- 10:-10

    sample <- 1:5
    expect_equal(calcGseaStat(stats, selectedStats = sample), 1)
    sample <- 16:21
    expect_equal(calcGseaStat(stats, selectedStats = sample), -1)

    sample <- (1:5)*2
    expect_equal(calcGseaStat(stats, selectedStats = sample), 0.71)
})

test_that("calcGseaStats(returnLeadingEdge=TRUE) works", {
    stats <- 10:-10

    sample <- c(10, 1:5)
    gseaRes <- calcGseaStat(stats, selectedStats = sample,
                            returnLeadingEdge = TRUE,
                            returnAllExtremes = TRUE)
    expect_equal(gseaRes$res, calcGseaStat(stats, selectedStats = sample))
    expect_true(1 %in% gseaRes$leadingEdge)
    expect_false(10 %in% gseaRes$leadingEdge)

    sample <- c(10, 16:21)
    gseaRes <- calcGseaStat(stats, selectedStats = sample,
                            returnLeadingEdge = TRUE,
                            returnAllExtremes = TRUE)
    expect_equal(gseaRes$res, calcGseaStat(stats, selectedStats = sample))
    expect_true(21 %in% gseaRes$leadingEdge)
    expect_false(10 %in% gseaRes$leadingEdge)

    sample <- 10:12
    gseaRes <- calcGseaStat(stats, selectedStats = sample,
                            returnLeadingEdge = TRUE,
                            returnAllExtremes = TRUE)
    expect_equal(gseaRes$res, calcGseaStat(stats, selectedStats = sample))
    expect_false(12 %in% gseaRes$leadingEdge)
    expect_false(10 %in% gseaRes$leadingEdge)
})

test_that("leading edge is consistent", {
    stats <- exampleRanks
    pathway <- examplePathways[[1]]

    set.seed(1)
    gseaRes <- fgseaSimple(list(p=pathway), stats, nperm=2)
    gseaResRev <- fgseaSimple(list(p=pathway), -stats, nperm=2)

    expect_identical(gseaRes$leadingEdge, gseaResRev$leadingEdge)

    gseaResPos <- fgseaSimple(list(p=pathway), stats, nperm=2, scoreType = "pos")
    expect_identical(gseaRes$leadingEdge, gseaResPos$leadingEdge)
})

test_that("calcGseaStats returns zero when both sides are equally enriched", {
    stats <- 10:-10
    sample <- 10:12
    expect_equal(calcGseaStat(stats, selectedStats = sample), 0)

    ess <- calcGseaStatCumulative(stats, selectedStats = sample, gseaParam = 1)
    for (i in seq_along(sample)) {
        expect_equal(ess[i], calcGseaStat(stats, sample[seq_len(i)]))
    }
})

test_that("calcGseaStats* work with zero gene-level stat", {

    expect_equal(
        calcGseaStat(c(10:1, 0, -1:-20), selectedStats = 11, 1),
        calcGseaStat(c(10:1, -1e-9, -1:-20), selectedStats = 11, 1))

    expect_equal(
        calcGseaStat(c(10:1, 0, 0, -1:-29), selectedStats = 11:12, 1),
        calcGseaStat(c(10:1, 0.1, -0.2, -1:-29), selectedStats = 11:12, 1))

    stats <- c(10:1, 0, 0, -1:-29)
    sample <- 11:13

    ess <- calcGseaStatCumulative(stats, selectedStats = sample, gseaParam = 1)
    for (i in seq_along(sample)) {
        expect_equal(ess[i], calcGseaStat(stats, sample[seq_len(i)]))
    }
})

test_that("calcGseaStatsCumulative works", {
    set.seed(42)
    stats <- rnorm(100)
    sample <- sample(seq_along(stats), 10)
    ess <- calcGseaStatCumulative(stats, selectedStats = sample, gseaParam = 1)
    for (i in seq_along(sample)) {
        expect_equal(ess[i], calcGseaStat(stats, sample[seq_len(i)]))
    }
})

test_that("fgsea results are reproducible with set.seed", {
    skip_on_os("mac")
    set.seed(42)
    data(examplePathways)
    data(exampleRanks)
    erm <- sample(exampleRanks, 1000)
    epw <- lapply(examplePathways, function(a) { return(a[a %in% names(erm)]) })
    epw <- epw[sapply(epw, length) >= 5]
    set.seed(42)
    res1 = fgseaSimple(pathways = epw, stats = erm, minSize=15, maxSize=500, nperm=1000, nproc=0)
    set.seed(42)
    res2 = fgseaSimple(pathways = epw, stats = erm, minSize=15, maxSize=500, nperm=1000, nproc=0)
    epsilon <- 1e-5
    for (i in seq_along(length(res1))) {
      expect_lte(abs(res1[i]$pval / res2[i]$pval - 1), epsilon)
    }
})


test_that(paste0("calcGseaStat and calcGseaStatCumulative calculate the same values for different",
                 " scoreType parameters"), {
    set.seed(42)
    data(exampleRanks)
    stats <- sort(exampleRanks, decreasing = TRUE)
    randomGeneSet <- sample(1:length(exampleRanks), size = 15)

    subSets <- lapply(1:length(randomGeneSet), function(x) randomGeneSet[1:x])

    scoresStd1 <- unlist(lapply(subSets, function(x) calcGseaStat(x, stats=stats,
                                                                  gseaParam=1,
                                                                  scoreType="std")))
    scoresStd2 <- fgsea:::calcGseaStatCumulative(stats, randomGeneSet, gseaParam=1, scoreType = "std")
    expect_equal(scoresStd1, scoresStd2)

    scoresPos1 <- unlist(lapply(subSets, function(x) calcGseaStat(x, stats=stats,
                                                                  gseaParam=1,
                                                                  scoreType="pos")))
    scoresPos2 <- fgsea:::calcGseaStatCumulative(stats, randomGeneSet, gseaParam=1, scoreType = "pos")
    expect_equal(scoresPos1, scoresPos2)

    scoresNeg1 <- unlist(lapply(subSets, function(x) calcGseaStat(x, stats=stats,
                                                                  gseaParam=1,
                                                                  scoreType="neg")))
    scoresNeg2 <- fgsea:::calcGseaStatCumulative(stats, randomGeneSet, gseaParam=1, scoreType = "neg")
    expect_equal(scoresNeg1, scoresNeg2)
})

test_that("preparePathwaysAndStats stops if stats contain NA values", {

  set.seed(42)
  data(exampleRanks)
  data(examplePathways)

  expect_silent(preparePathwaysAndStats(examplePathways, exampleRanks, minSize = 1,
                                        maxSize = 500, gseaParam = 1, scoreType = "std"))

  ranks1 <- exampleRanks
  ranks1[sample(seq_along(ranks1), 10)] <- NA

  expect_error(preparePathwaysAndStats(examplePathways, ranks1, minSize = 1,
                                       maxSize = 500, gseaParam = 1, scoreType = "std"))

  ranks2 <- exampleRanks
  ranks2[sample(seq_along(ranks2), 1)] <- NA

  expect_error(preparePathwaysAndStats(examplePathways, ranks2, minSize = 1,
                                       maxSize = 500, gseaParam = 1, scoreType = "std"))
})


test_that("preparePathwaysAndStats stops if stats contain infinite values", {

  set.seed(42)
  data(exampleRanks)
  data(examplePathways)

  originalRanks <- exampleRanks
  expect_silent(preparePathwaysAndStats(examplePathways, originalRanks, minSize = 1,
                                        maxSize = 500, gseaParam = 1, scoreType = "std"))

  ranks1 <- originalRanks

  indx <- sample(seq_along(ranks1), 1)
  ranks1[indx] <- Inf
  expect_error(preparePathwaysAndStats(examplePathways, ranks1, minSize = 1,
                                       maxSize = 500, gseaParam = 1, scoreType = "std"))

  ranks1[indx] <- -Inf
  expect_error(preparePathwaysAndStats(examplePathways, ranks1, minSize = 1,
                                       maxSize = 500, gseaParam = 1, scoreType = "std"))
  ranks2 <- originalRanks
  indxs <- sample(seq_along(ranks2), 10)
  ranks2[indxs] <- Inf

  expect_error(preparePathwaysAndStats(examplePathways, ranks2, minSize = 1,
                                       maxSize = 500, gseaParam = 1, scoreType = "std"))

  ranks2[indxs] <- -Inf
  expect_error(preparePathwaysAndStats(examplePathways, ranks2, minSize = 1,
                                       maxSize = 500, gseaParam = 1, scoreType = "std"))

  ranks2[indxs[1:5]] <- Inf
  ranks2[indxs[6:10]] <- -Inf
  expect_error(preparePathwaysAndStats(examplePathways, ranks2, minSize = 1,
                                       maxSize = 500, gseaParam = 1, scoreType = "std"))
})
