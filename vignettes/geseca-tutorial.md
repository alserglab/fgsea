An R-package for performing GEne SEt Co-regulation Analysis (GESECA).

## Quick run

Loading necessary libraries

``` r
library(fgsea)
library(data.table)
```

Loading example expression matrix and pathways:

``` r
data(exampleExpressionMatrix)
str(exampleExpressionMatrix)
#>  num [1:10000, 1:12] 16.4 16.3 16.2 16.2 16.3 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:10000] "20090" "22121" "20091" "67671" ...
#>   ..$ : chr [1:12] "GSM357839" "GSM357841" "GSM357842" "GSM357843" ...
```

`exampleExpressionMatrix` - numeric expression matrix of the GSE14308
dataset. Rows correspond to genes (ENTREZID is used as identifiers),
columns correspond to samples. **Note:** `geseca` expects that the
values of expression matrix are log-scaled and the samples are quantile
normalized. To get more information about steps that were performed to
process the GSE14308 expression matrix, use: `?exampleExpressionMatrix`.

``` r
data(examplePathways)
str(head(examplePathways, 3))
#> List of 3
#>  $ 1221633_Meiotic_Synapsis                                     : chr [1:64] "12189" "13006" "15077" "15078" ...
#>  $ 1368092_Rora_activates_gene_expression                       : chr [1:9] "11865" "12753" "12894" "18143" ...
#>  $ 1368110_Bmal1:Clock,Npas2_activates_circadian_gene_expression: chr [1:16] "11865" "11998" "12753" "12952" ...
```

`examplePathways` - list of mouse pathways from `reactome.db` package.
The list of gene sets should use the same identifiers that are used in
`rownames(exampleExpressionMatrix)`. To get more information about this
list, use: `?examplePathways`

Running `geseca`:

``` r
gesecaRes <- geseca(exampleExpressionMatrix, examplePathways, minSize = 15, maxSize = 500)
head(gesecaRes[order(pval), ], 5)
#>                                            pathway   pctEvar         pval
#> 1:                              5990980_Cell_Cycle 1.0525898 2.276767e-39
#> 2:                     5990979_Cell_Cycle,_Mitotic 1.0466758 3.891185e-39
#> 3:                    5991851_Mitotic_Prometaphase 0.7485590 1.409623e-24
#> 4: 5992217_Resolution_of_Sister_Chromatid_Cohesion 0.7075162 1.265783e-22
#> 5:                                 5991454_M_Phase 0.6020886 6.237558e-22
#>            padj  log2err size
#> 1: 2.587638e-37 1.634419  417
#> 2: 2.587638e-37 1.628072  362
#> 3: 6.249329e-23 1.287104   96
#> 4: 4.208729e-21 1.229504   88
#> 5: 1.659190e-20 1.212545  211
```

Optionally, it is possible to scale the values of each gene-row to a
unit variance by passing the `scale=TRUE` to `geseca`.
