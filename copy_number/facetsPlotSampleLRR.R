plotSampleLRR <- function(x, fit)
{
    mat = x$jointseg
    cncf = fit$cncf
    dipLogR <- fit$dipLogR
    par(mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0))
    chr = mat$chrom
    len = table(chr)
    altcol = rep_len(c("light blue", "gray"), length(len))
    chr.col = rep(altcol, len)
    nmark = cncf$num.mark
    tmp = cumsum(len)
    start = c(1, tmp[-length(len)] + 1)
    end = tmp
    mid = start + len/2
    plot(mat$cnlr, pch = ".", axes = F, cex = 1.5, ylim = c(-max(c(max(abs(mat$cnlr)), 4)),
                                                            max(c(max(abs(mat$cnlr)), 4))), col = c("grey", "lightblue")[1 + rep(cncf$chrom -
                                                                                                                         2 * floor(cncf$chrom/2), cncf$num.mark)],
         ylab = "log-ratio", xlab="Chromosomes")
     points(rep(cncf$cnlr.median, cncf$num.mark), pch = ".", cex = 2,
            col = "brown")
    labs <- names(mid)
    labs <- sub('21', '', labs)
    labs <- sub('23', 'X', labs)
    axis(side = 1, at = mid, labs, cex.axis = 1, las = 2)
    axis(side = 2, cex.axis = 1, las=2)
    abline(h=0, lty=2, col="lightgrey")
    box()
}
