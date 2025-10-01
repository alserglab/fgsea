cowplotText <- function(text, style) {
    ggdraw() + do.call(cowplot::draw_text, c(list(text=text), style))
}

cowplotLabel <- function(label, style) {
    ggdraw() + do.call(cowplot::draw_label, c(list(label=label), style))
}

valueToExpExpression <- function(val, log2err=NULL) {
    a <- floor(log10(val))
    b <- val/(10**a)
    x <- sprintf("%.1f\u00B710", b)
    if (is.null(log2err)) {
        y <- as.character(a)
    } else {
        y <- sprintf("%s\u00B1%.2f", a, log2err / log2(10))
    }

    as.expression(substitute(x ^ y, list(y=y,  x=x)))
}


#' Plots table of enrichment graphs using ggplot and gridExtra.
#' @param pathways Pathways to plot table, as in `fgsea` function.
#' @param stats Gene-level stats, as in `fgsea` function.
#' @param fgseaRes Table with fgsea results.
#' @param gseaParam GSEA-like parameter. Adjusts displayed statistic values,
#'      values closer to 0 flatten plots. Default = 1, value of 0.5 is a good
#'      choice too.
#' @param colwidths Vector of five elements corresponding to column width for
#'      grid.arrange. Can be both units and simple numeric vector, in latter case
#'      it defines proportions, not actual sizes. If column width is set to zero, the column is not drawn.
#' @param pathwayLabelStyle list with style parameter adjustments for pathway labels.
#'      For example, `list(size=10, color="red")` set the font size to 10 and color to red.
#'      See `cowplot::draw_text` for possible options.
#' @param headerLabelStyle similar to `pathwayLabelStyle` but for the table header.
#' @param valueStyle similar to `pathwayLabelStyle` but for NES and p-value columns.
#' @param axisLabelStyle list with style parameter adjustments for stats axis labels.
#'      See `ggplot2::element_text` for possible options.
#' @param render (deprecated)
#' @return ggplot object with enrichment barcode plots
#' @import ggplot2
#' @import grid
#' @import cowplot
#' @export
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' fgseaRes <- fgsea(examplePathways, exampleRanks, minSize=15, maxSize=500)
#' topPathways <- fgseaRes[head(order(pval), n=15)][order(NES), pathway]
#' plotGseaTable(examplePathways[topPathways], exampleRanks,
#'               fgseaRes, gseaParam=0.5)
plotGseaTable <- function(pathways, stats, fgseaRes,
                          gseaParam=1,
                          colwidths=c(5, 3, 0.8, 1.2, 1.2),
                          pathwayLabelStyle=NULL,
                          headerLabelStyle=NULL,
                          valueStyle=NULL,
                          axisLabelStyle=NULL,
                          render=NULL) {

    pathwayLabelStyleDefault <- list(size=12, hjust=1, x=0.95, vjust=0)
    pathwayLabelStyle <- modifyList(pathwayLabelStyleDefault, as.list(pathwayLabelStyle))

    headerLabelStyleDefault <- list(size=12)
    headerLabelStyle <- modifyList(headerLabelStyleDefault, as.list(headerLabelStyle))

    valueStyleDefault <- list(size=12, vjust=0)
    valueStyle <- modifyList(valueStyleDefault, as.list(valueStyle))

    if (!is.null(render)) {
        warning("render argument is deprecated, a ggplot object is always returned")
    }

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    pathways <- lapply(pathways, function(p) {
        unname(as.vector(na.omit(match(p, names(statsAdj)))))
    })

    # fixes #40
    pathways <- pathways[sapply(pathways, length) > 0]

    ps <- lapply(names(pathways), function(pn) {
        p <- pathways[[pn]]
        annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
        list(
            cowplotText(pn, pathwayLabelStyle),
            ggplot() +
                geom_segment(aes(x=p, xend=p,
                                 y=0, yend=statsAdj[p]),
                             linewidth=0.2) +
                scale_x_continuous(limits=c(0, length(statsAdj)),
                                   expand=c(0, 0)) +
                scale_y_continuous(limits=c(-1, 1),
                                   expand=c(0, 0)) +
                xlab(NULL) + ylab(NULL) +
                theme(panel.background = element_blank(),
                      plot.background=element_blank(),
                      axis.line=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      panel.grid = element_blank(),
                      axis.title=element_blank(),
                      plot.margin = rep(unit(0,"null"),4),
                      panel.spacing = rep(unit(0,"null"),4)
                ),
            cowplotText(sprintf("%.2f", annotation$NES), valueStyle),
            cowplotLabel(valueToExpExpression(annotation$pval), valueStyle),
            cowplotLabel(valueToExpExpression(annotation$padj), valueStyle)
            )
    })

    rankPlot <-
        ggplot() +
        geom_blank() +
        scale_x_continuous(limits=c(0, length(statsAdj)),
                           expand=c(0, 0)) +
        scale_y_continuous(limits=c(-1, 1),
                           expand=c(0, 0)) +
        xlab(NULL) + ylab(NULL) +
        theme(panel.background = element_blank(),
              plot.background=element_blank(),
              axis.line=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid = element_blank(),
              axis.title=element_blank(),
              axis.text=do.call(element_text, as.list(axisLabelStyle)),
              plot.margin = unit(c(0,0,0.5,0), "npc"),
              panel.spacing = unit(c(0,0,0,0), "npc")
        )

    grobs <- c(
        list(cowplotText("Pathway",
                         modifyList(headerLabelStyle, pathwayLabelStyle[c("hjust", "x")])
                         )),
        lapply(c("Gene ranks", "NES", "pval", "padj"), cowplotText, style=headerLabelStyle),
        unlist(ps, recursive = FALSE),
        list(nullGrob(),
             rankPlot,
             nullGrob(),
             nullGrob(),
             nullGrob()))

    # not drawing column if corresponding colwidth is set to zero
    grobsToDraw <- rep(as.numeric(colwidths) != 0, length(grobs)/length(colwidths))


    p <- plot_grid(plotlist=grobs[grobsToDraw],
                   ncol=sum(as.numeric(colwidths) != 0),
                   rel_widths=colwidths[as.numeric(colwidths) != 0])

    p
}


#' Returns data required for doing an enrichment plot.
#' @param pathway Gene set to plot.
#' @param stats Gene-level statistics.
#' @param gseaParam GSEA parameter.
#' @return returns list with the following data:
#' * `curve` - data.table with the coordinates of the enrichment curve;
#' * `ticks` - data.table with statistic entries for each pathway gene, adjusted with gseaParam;
#' * `stats` - data.table with statistic values for all of the genes, adjusted with gseaParam;
#' * `posES`, `negES`, `spreadES` - values of the positive enrichment score,
#'  negative enrichment score, and difference between them;
#' * `maxAbsStat` - maximal absolute value of statistic entries, adjusted with gseaParam
#'
#' Note: for input `stats=0` and `gseaParam=0`, the returned `ticks` and `stats` will be zero, while
#' values of one are used in calculation of enrichment score.
#' @export
#' @examples
#' library(ggplot2)
#' data(examplePathways)
#' data(exampleRanks)
#'
#' pd <- plotEnrichmentData(
#'     pathway = examplePathways[["5991130_Programmed_Cell_Death"]],
#'     stats = exampleRanks
#' )
#'
#' with(pd,
#'      ggplot(data=curve) +
#'          geom_line(aes(x=rank, y=ES), color="green") +
#'          geom_ribbon(data=stats,
#'                      mapping=aes(x=rank, ymin=0,
#'                                  ymax=stat/maxAbsStat*(spreadES/4)),
#'                      fill="grey") +
#'          geom_segment(data=ticks,
#'                       mapping=aes(x=rank, y=-spreadES/16,
#'                                   xend=rank, yend=spreadES/16),
#'                       size=0.2) +
#'          geom_hline(yintercept=posES, colour="red", linetype="dashed") +
#'          geom_hline(yintercept=negES, colour="red", linetype="dashed") +
#'          geom_hline(yintercept=0, colour="black") +
#'          theme(
#'              panel.background = element_blank(),
#'              panel.grid.major=element_line(color="grey92")
#'          ) +
#'          labs(x="rank", y="enrichment score"))
plotEnrichmentData <- function(pathway, stats,
                              gseaParam=1) {

    if (any(!is.finite(stats))){
        stop("Not all stats values are finite numbers")
    }

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statSigns <- sign(statsAdj)
    statsAdj <- (abs(statsAdj) ^ gseaParam)

    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    pathway <- unique(pathway)

    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                                   returnAllExtremes = TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.table(rank=c(0, xs, n + 1), ES=c(0, ys, 0))
    ticks <- data.table(rank=pathway, stat=statSigns[pathway]*statsAdj[pathway])
    stats <- data.table(rank=seq_along(stats), stat=statSigns*statsAdj)

    res <- list(
        curve=toPlot,
        ticks=ticks,
        stats=stats,
        posES=max(tops),
        negES=min(bottoms),
        spreadES=max(tops)-min(bottoms),
        maxAbsStat=max(abs(statsAdj)))
}

#' Plots GSEA enrichment plot. For more flexibility use `plotEnrichmentData` function.
#' @param pathway Gene set to plot.
#' @param stats Gene-level statistics.
#' @param gseaParam GSEA parameter.
#' @param ticksSize width of vertical line corresponding to a gene (default: 0.2)
#' @return ggplot object with the enrichment plot.
#' @export
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' \dontrun{
#' plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
#'                exampleRanks)
#' }
plotEnrichment <- function(pathway, stats,
                          gseaParam=1,
                          ticksSize=0.2) {

    pd <- plotEnrichmentData(
        pathway = pathway,
        stats = stats,
        gseaParam = gseaParam)

    with(pd,
         ggplot(data=curve) +
             geom_line(aes(x=rank, y=ES), color="green") +
             geom_segment(data=ticks,
                          mapping=aes(x=rank, y=-spreadES/16,
                                      xend=rank, yend=spreadES/16),
                          linewidth=ticksSize) +
             geom_hline(yintercept=posES, colour="red", linetype="dashed") +
             geom_hline(yintercept=negES, colour="red", linetype="dashed") +
             geom_hline(yintercept=0, colour="black") +
             theme(
                 panel.background = element_blank(),
                 panel.grid.major=element_line(color="grey92")
             ) +
             labs(x="rank", y="enrichment score"))
}
