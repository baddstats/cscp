## -------- Base R schematic: clearer colours + emphasis on CSCP modularisation
## Writes: semilog_gminus1_base_colours.pdf / .png

## Parameters (same as before)
phiL <- 0.60; phiS <- 0.06       # large/small ranges
aL <- 0.55;  aS <- 0.95          # LGCP weights (latent)
bL <- 0.40;  bS <- 1.40          # CSCP bumps in g-1

hstar <- (phiS*phiL)/(2*(phiL - phiS)) * log(bS/bL)
r <- seq(1e-4, 1, length.out = 1400)

g_cs_total <- 1 + bL*exp(-2*r/phiL) + bS*exp(-2*r/phiS)
g_cs_small <-      bS*exp(-2*r/phiS)
g_cs_large <-      bL*exp(-2*r/phiL)
g_lg       <- exp(aL*exp(-r/phiL) + aS*exp(-r/phiS))

## Okabe–Ito colourblind-safe palette
col_cscp   <- "#D55E00"  # vermillion (CSCP total)
col_lgcp   <- "#000000"  # black     (LGCP total)
col_small  <- "#0072B2"  # blue      (CSCP small-scale component)
col_large  <- "#009E73"  # green     (CSCP large-scale component)
col_hstar  <- "#555555"  # dark gray

plot_semilog <- function(file, type=c("pdf","png")) {
  type <- match.arg(type)
  if (type=="pdf") pdf(file, width=7.4, height=4.6, useDingbats=FALSE)
  if (type=="png") png(file, width=1700, height=1050, res=220)
  
  op <- par(mar=c(4.6,5,1.2,12), mgp=c(2.8,0.9,0))  # extra right margin for legend
  on.exit({par(op); if (dev.cur()>1) dev.off()}, add=TRUE)
  
  ylo <- 1e-3; yhi <- 10; yt <- 10^seq(-3,1,1)
  
  plot(NA, NA, xlim=c(0,1), ylim=c(ylo,yhi), log="y",
       xlab="Distance r", ylab=expression(g(r) - 1~"(log scale)"),
       axes=FALSE, xaxs="i", yaxs="i")
  abline(v=seq(0,1,0.1), col=gray(0.92), lwd=1)
  abline(h=yt,               col=gray(0.92), lwd=1)
  axis(1, at=seq(0,1,0.1))
  axis(2, at=yt, labels=parse(text=paste0("10^", seq(-3,1))))
  box()
  
  ## Draw totals bold, components dashed & thinner
  lines(r, g_lg - 1,       lwd=3.5, col=col_lgcp)          # LGCP total (black)
  lines(r, g_cs_total - 1, lwd=3.5, col=col_cscp)          # CSCP total (vermillion)
  
  lines(r, g_cs_small,     lwd=2.2, col=col_small, lty=2)  # CSCP small comp (blue)
  lines(r, g_cs_large,     lwd=2.2, col=col_large, lty=4)  # CSCP large comp (green, dotdash)
  
  ## Shoulder marker & label
  abline(v=hstar, col=col_hstar, lty=3, lwd=2)
  text(hstar, 0.018, expression(h^"*"), col=col_hstar, pos=4, cex=0.95)
  
  ## Compact legend in the right margin (no overlap)
  par(xpd=NA)  # allow drawing in outer margin
  legend(1.02, yhi, xpd=NA, bty="n", cex=0.95,
         lwd=c(3.5,3.5,2.2,2.2,2),
         lty=c(1,1,2,4,3),
         col=c(col_cscp, col_lgcp, col_small, col_large, col_hstar),
         legend=c("CSCP (2 scales) — total",
                  "LGCP (2 scales) — total",
                  "CSCP component (small scale)",
                  "CSCP component (large scale)",
                  expression(h^"*~(equal~component~influence)")))
}

plot_semilog("semilog_gminus1_base_colours.pdf", "pdf")
plot_semilog("semilog_gminus1_base_colours.png", "png")
cat("Wrote: semilog_gminus1_base_colours.pdf/png\n")