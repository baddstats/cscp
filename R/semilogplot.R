############################################################
## Two-scale CSCP vs LGCP (schematic) + two-panel figure  ##
## + helper to estimate phiS, phiL from log{g-1} lines    ##
## Base R only                                            ##
############################################################

## ---- Parameters (choose well-separated scales to see the shoulder clearly)
phiL <- 0.60; phiS <- 0.06             # large/small ranges (exponential)
aL   <- 0.55; aS   <- 0.95             # LGCP latent weights (in exponent)
bL   <- 0.40; bS   <- 1.40             # CSCP bumps in g-1
r    <- seq(1e-4, 1, length.out = 1400)

## ---- Analytic curves
g_cs  <- 1 + bL*exp(-2*r/phiL) + bS*exp(-2*r/phiS)         # CSCP total
g_lg  <- exp(aL*exp(-r/phiL) + aS*exp(-r/phiS))            # LGCP total
gS    <-      bS*exp(-2*r/phiS)                            # CSCP small   (g-1)
gL    <-      bL*exp(-2*r/phiL)                            # CSCP large   (g-1)
hstar <- (phiS*phiL)/(2*(phiL - phiS)) * log(bS/bL)        # equal-contribution radius

## ---- Okabe–Ito palette (colourblind-safe)
col_cscp  <- "#D55E00"   # CSCP total (vermillion)
col_lgcp  <- "#000000"   # LGCP total (black)
col_small <- "#0072B2"   # CSCP small component (blue)
col_large <- "#009E73"   # CSCP large component (green)
col_hstar <- "#555555"   # marker

## ---------- (A) TWO-PANEL FIGURE (base R only)
two_panel <- function(file, type=c("pdf","png")) {
  type <- match.arg(type)
  if (type=="pdf") pdf(file, width=9.2, height=4.6, useDingbats=FALSE)
  if (type=="png") png(file, width=2000, height=1000, res=220)
  
  op <- par(mfrow=c(1,2),
            mar=c(4.6,5,1.2,1.2), mgp=c(2.8,0.9,0))
  on.exit({par(op); if (dev.cur()>1) dev.off()}, add=TRUE)
  
  ## ---- Left: linear g(r)
  plot(NA, NA, xlim=c(0,1), ylim=range(c(1,g_cs,g_lg)),
       xlab="Distance r", ylab=expression(g(r)), xaxs="i", yaxs="i")
  abline(v=seq(0,1,0.1), col=gray(0.93), lwd=1)
  abline(h=pretty(range(c(1,g_cs,g_lg))), col=gray(0.93), lwd=1)
  lines(r, g_lg,        lwd=3.5, col=col_lgcp)
  lines(r, g_cs,        lwd=3.5, col=col_cscp)
  abline(v=hstar, col=col_hstar, lty=3, lwd=2)
  text(hstar, 1.02*min(g_cs), expression(h^"*"), col=col_hstar, pos=4, cex=0.95)
  legend("topright", bty="n", cex=0.95,
         lwd=c(3.5,3.5,2), lty=c(1,1,3), col=c(col_cscp,col_lgcp,col_hstar),
         legend=c("CSCP (2 scales)","LGCP (2 scales)", expression(h^"*")))
  
  ## ---- Right: semilog g(r)-1
  ylo <- 1e-3; yhi <- 10; yt <- 10^seq(-3,1,1)
  plot(NA, NA, xlim=c(0,1), ylim=c(ylo,yhi), log="y",
       xlab="Distance r", ylab=expression(g(r) - 1~"(log scale)"),
       axes=FALSE, xaxs="i", yaxs="i")
  abline(v=seq(0,1,0.1), col=gray(0.93), lwd=1)
  abline(h=yt,               col=gray(0.93), lwd=1)
  axis(1, at=seq(0,1,0.1))
  axis(2, at=yt, labels=parse(text=paste0("10^", seq(-3,1))))
  box()
  
  lines(r, g_lg - 1,  lwd=3.5, col=col_lgcp)
  lines(r, g_cs - 1,  lwd=3.5, col=col_cscp)
  lines(r, gS,        lwd=2.2, col=col_small, lty=2)   # small-scale comp
  lines(r, gL,        lwd=2.2, col=col_large, lty=4)   # large-scale comp
  
  abline(v=hstar, col=col_hstar, lty=3, lwd=2)
  text(hstar, 0.018, expression(h^"*"), col=col_hstar, pos=4, cex=0.95)
  
  legend("topright", inset=0.02, bty="n", cex=0.95,
         lwd=c(3.5,3.5,2.2,2.2,2),
         lty=c(1,1,2,4,3),
         col=c(col_cscp,col_lgcp,col_small,col_large,col_hstar),
         legend=c("CSCP total","LGCP total",
                  "CSCP small component","CSCP large component", expression(h^"*")))
}

two_panel("pcf_twopanel_base.pdf", "pdf")
two_panel("pcf_twopanel_base.png", "png")
cat("Wrote: pcf_twopanel_base.pdf/png\n")

## ---------- (B) FIT TWO LINES TO log{g-1} TO ESTIMATE phiS, phiL -------------
## Given vectors r (distances) and y = g(r)-1, and two r-intervals where each
## component dominates, fit log(y) ~ r in each interval.
## Returns slopes m1, m2 (so phi = -2/m), intercepts (log b), and predictions.

fit_two_lines_log_gminus1 <- function(r, y, range_small, range_large) {
  stopifnot(length(r)==length(y))
  keep <- which(is.finite(r) & is.finite(y) & y > 0)
  r <- r[keep]; y <- y[keep]
  
  idxS <- which(r >= range_small[1] & r <= range_small[2] & y > 0)
  idxL <- which(r >= range_large[1] & r <= range_large[2] & y > 0)
  if (length(idxS) < 5 || length(idxL) < 5)
    stop("Not enough points in one of the fitting ranges.")
  
  fitS <- lm(log(y[idxS]) ~ r[idxS])
  fitL <- lm(log(y[idxL]) ~ r[idxL])
  
  mS <- coef(fitS)[2];  mL <- coef(fitL)[2]              # slopes
  bS_hat <- exp(coef(fitS)[1]);  bL_hat <- exp(coef(fitL)[1])
  phiS_hat <- -2 / mS;  phiL_hat <- -2 / mL
  
  list(
    phi_small = phiS_hat,
    phi_large = phiL_hat,
    b_small   = bS_hat,
    b_large   = bL_hat,
    fit_small = fitS,
    fit_large = fitL,
    predict   = function(rnew) bS_hat*exp(-2*rnew/phiS_hat) + bL_hat*exp(-2*rnew/phiL_hat)
  )
}

## ---- Demo on the analytic CSCP curve (pretend it's empirical)
y_cs <- g_cs - 1
est <- fit_two_lines_log_gminus1(
  r, y_cs,
  range_small = c(0.015, 0.05),   # region dominated by small-scale term
  range_large = c(0.25, 0.8)      # region dominated by large-scale term
)

print(est[c("phi_small","phi_large","b_small","b_large")])

## ---- (Optional) overlay the fitted two-line model on the semilog panel
overlay_fit <- function(file, type=c("pdf","png")) {
  type <- match.arg(type)
  if (type=="pdf") pdf(file, width=6.4, height=4.6, useDingbats=FALSE)
  if (type=="png") png(file, width=1500, height=1000, res=220)
  
  op <- par(mar=c(4.6,5,1.2,1.2), mgp=c(2.8,0.9,0))
  on.exit({par(op); if (dev.cur()>1) dev.off()}, add=TRUE)
  
  ylo <- 1e-3; yhi <- 10; yt <- 10^seq(-3,1,1)
  plot(NA, NA, xlim=c(0,1), ylim=c(ylo,yhi), log="y",
       xlab="Distance r", ylab=expression(g(r) - 1~"(log scale)"),
       axes=FALSE, xaxs="i", yaxs="i")
  abline(v=seq(0,1,0.1), col=gray(0.93), lwd=1)
  abline(h=yt,               col=gray(0.93), lwd=1)
  axis(1, at=seq(0,1,0.1))
  axis(2, at=yt, labels=parse(text=paste0("10^", seq(-3,1))))
  box()
  
  lines(r, y_cs,                 lwd=3.2, col=col_cscp)                         # true CSCP total
  lines(r, est$predict(r),       lwd=2.4, col=col_cscp, lty=3)                   # two-line fit
  lines(r, bS*exp(-2*r/phiS),    lwd=2.0, col=col_small, lty=2)                  # true small comp
  lines(r, bL*exp(-2*r/phiL),    lwd=2.0, col=col_large, lty=4)                  # true large comp
  
  legend("topright", bty="n", cex=0.95,
         lwd=c(3.2,2.4,2,2), lty=c(1,3,2,4),
         col=c(col_cscp, col_cscp, col_small, col_large),
         legend=c("CSCP total (truth)","Two-line fit",
                  "Small component (truth)","Large component (truth)"))
}

overlay_fit("pcf_semilog_fit_overlay.pdf", "pdf")
overlay_fit("pcf_semilog_fit_overlay.png", "png")
cat("Wrote: pcf_semilog_fit_overlay.pdf/png\n")






## ---------- Parameters (same spirit as before; well-separated scales)
phiL <- 0.60; phiS <- 0.06          # ranges
aL <- 0.55;  aS <- 0.95             # LGCP weights (in exponent)
bL <- 0.40;  bS <- 1.40             # CSCP bumps in g-1
r  <- seq(1e-4, 1, length.out = 1400)

## Curves
g_cs_total <- 1 + bL*exp(-2*r/phiL) + bS*exp(-2*r/phiS)     # CSCP total
g_cs_small <-      bS*exp(-2*r/phiS)                        # (g-1) small
g_cs_large <-      bL*exp(-2*r/phiL)                        # (g-1) large
g_lg       <- exp(aL*exp(-r/phiL) + aS*exp(-r/phiS))        # LGCP total
hstar <- (phiS*phiL)/(2*(phiL - phiS)) * log(bS/bL)

## Two-line fit to log{g-1} for CSCP (small- & large-dominated bands)
fit_two_lines_log_gminus1 <- function(r, y, range_small, range_large){
  keep <- which(is.finite(r) & is.finite(y) & y > 0); r <- r[keep]; y <- y[keep]
  idxS <- which(r >= range_small[1] & r <= range_small[2] & y > 0)
  idxL <- which(r >= range_large[1] & r <= range_large[2] & y > 0)
  fitS <- lm(log(y[idxS]) ~ r[idxS]); fitL <- lm(log(y[idxL]) ~ r[idxL])
  mS <- coef(fitS)[2]; mL <- coef(fitL)[2]
  bS.h <- exp(coef(fitS)[1]); bL.h <- exp(coef(fitL)[1])
  phiS.h <- -2/mS; phiL.h <- -2/mL
  list(phi_small=phiS.h, phi_large=phiL.h, b_small=bS.h, b_large=bL.h,
       predict=function(rr) bS.h*exp(-2*rr/phiS.h) + bL.h*exp(-2*rr/phiL.h))
}
y_cs <- g_cs_total - 1
est  <- fit_two_lines_log_gminus1(r, y_cs,
                                  range_small=c(0.015,0.05),
                                  range_large=c(0.25,0.8))

## Colours (Okabe–Ito palette)
col_cscp  <- "#D55E00"  # CSCP total – orange
col_lgcp  <- "#000000"  # LGCP total – black
col_small <- "#0072B2"  # CSCP small comp – blue
col_large <- "#009E73"  # CSCP large comp – green
col_fit   <- "#D55E00"  # same hue, dotted
col_hstar <- "#555555"

plot_semilog_overlay <- function(file, type=c("pdf","png")){
  type <- match.arg(type)
  if (type=="pdf") pdf(file, width=7.4, height=4.6, useDingbats=FALSE)
  if (type=="png") png(file, width=1700, height=1050, res=220)
  
  op <- par(mar=c(4.6,5,1.2,12), mgp=c(2.8,0.9,0))
  on.exit({par(op); if (dev.cur()>1) dev.off()}, add=TRUE)
  
  ylo <- 1e-3; yhi <- 10; yt <- 10^seq(-3,1,1)
  plot(NA, NA, xlim=c(0,1), ylim=c(ylo,yhi), log="y",
       xlab="Distance r", ylab=expression(g(r) - 1~"(log scale)"),
       axes=FALSE, xaxs="i", yaxs="i")
  abline(v=seq(0,1,0.1), col=gray(0.93), lwd=1)
  abline(h=yt,               col=gray(0.93), lwd=1)
  axis(1, at=seq(0,1,0.1))
  axis(2, at=yt, labels=parse(text=paste0("10^", seq(-3,1))))
  box()
  
  ## Totals first (bold), then components, then fit
  lines(r, g_lg - 1,          lwd=3.6, col=col_lgcp)     # LGCP total
  lines(r, g_cs_total - 1,    lwd=3.6, col=col_cscp)     # CSCP total
  lines(r, g_cs_small,        lwd=2.3, col=col_small, lty=2)   # CSCP small comp
  lines(r, g_cs_large,        lwd=2.3, col=col_large, lty=4)   # CSCP large comp
  lines(r, est$predict(r),    lwd=2.6, col=col_fit,   lty=3)   # two-line fit
  
  abline(v=hstar, col=col_hstar, lty=3, lwd=2)
  text(hstar, 0.018, expression(h^"*"), col=col_hstar, pos=4, cex=0.95)
  
  par(xpd=NA)
  legend(1.02, yhi, xpd=NA, bty="n", cex=0.95,
         lwd=c(3.6,3.6,2.6,2.3,2.3,2),
         lty=c(1,1,3,2,4,3),
         col=c(col_cscp,col_lgcp,col_fit,col_small,col_large,col_hstar),
         legend=c("CSCP total (2 scales)",
                  "LGCP total (2 scales)",
                  "Two-line fit to CSCP",
                  "CSCP small component",
                  "CSCP large component",
                  expression(h^"*~(equal~components)")))
}

plot_semilog_overlay("semilog_gminus1_with_LGCP.pdf", "pdf")
plot_semilog_overlay("semilog_gminus1_with_LGCP.png", "png")
cat("Wrote: semilog_gminus1_with_LGCP.pdf/png\n")