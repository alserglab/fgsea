# code to prepare `exampleExpressionMatrix` dataset goes here
expressionMatrixFile <- "./inst/extdata/gse14308_expression_matrix.tsv"
exampleExpressionMatrix <- as.matrix(read.table(expressionMatrixFile, sep = "\t"))

usethis::use_data(exampleExpressionMatrix, overwrite = TRUE)
