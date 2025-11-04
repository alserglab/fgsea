#' @useDynLib fgsea
#' @import Rcpp
NULL

#' Example vector of gene-level statistics obtained for Th1 polarization.
#'
#' The data were obtained by doing differential expression between
#' Naive and Th1-activated states for GEO dataset GSE14308. The
#' exact script is available as system.file("gen_gene_ranks.R", package="fgsea")
#' @docType data
#' @name exampleRanks
NULL

#' Example of expression values obtained for GSE14308.
#'
#' Expression data was obtained by preprocessing the GSE14308 dataset.
#' For the matrix of gene expression value, the following steps were performed:
#' \itemize{
#' \item expression values were log2-scaled
#' \item quantile-type normalization was perfomred between arrays
#' \item rows were collapsed by `ENTREZID`
#' \item rows were sorted in descending order by mean expression value per gene
#' \item finally, top-10_000 genes were taken
#' }
#' The exact script is available as
#' system.file("gen_gse14308_expression_matrix.R", package="fgsea")
#' @docType data
#' @name exampleExpressionMatrix
NULL

#' Example list of mouse Reactome pathways.
#'
#' The list was obtained by selecting all the pathways from `reactome.db`
#' package that contain mouse genes. The exact script is available as
#' system.file("gen_reactome_pathways.R", package="fgsea")
#' @docType data
#' @name examplePathways
NULL
