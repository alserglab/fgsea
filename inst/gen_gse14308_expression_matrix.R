library(GEOquery)
library(limma)
# for collapseBy
source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")

gse14308 <- getGEO("GSE14308")[[1]]

pData(gse14308)$condition <- sub("-.*$", "", gse14308$title)

gse14308 <- gse14308[!rownames(gse14308) %in% c("NONE", "NA", "none", "NULL")]
exprs(gse14308) <- log2(exprs(gse14308) + 1)
exprs(gse14308) <- normalizeBetweenArrays(exprs(gse14308), method = "quantile")

es <- collapseBy(gse14308, fData(gse14308)$ENTREZ_GENE_ID, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]
es <- exprs(es)

es <- es[head(order(apply(es, 1, mean), decreasing = TRUE), 10000), ]
write.table(es, "inst/extdata/gse14308_expression_matrix.tsv",
            sep="\t", quote = FALSE)
