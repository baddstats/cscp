
DATADIR <- "./data-auto"
FIGDIR <- "./pix-auto"
if(!dir.exists(DATADIR)) dir.create(DATADIR)
if(!dir.exists(FIGDIR)) dir.create(FIGDIR)

plotit <- function(tag, expr, w=6, h=6, dir=FIGDIR) {
  fnpdf <- paste0(dir, "/", tag, ".pdf")
  pdf(fnpdf, width=w, height=h)
  eval(substitute(expr))
  dev.off()
  fneps <- paste0(dir, "/", tag, ".eps")
  postscript(fneps, width=w, height=h, horizontal=FALSE)
  eval(substitute(expr))
  dev.off()
  return(c(fnpdf, fneps))
}

plotit("chisq", {
  curve(dchisq(x, df=1), to=10, ylim=c(0, 0.5),
        ylab="probability density", main="Chi-squared")
  curve(dchisq(x, df=2), add=TRUE, col=2)
  curve(dchisq(x, df=3), add=TRUE, col=3)
  curve(dchisq(x, df=4), add=TRUE, col=4)
  legend("topright", lty=1, col=1:4,
         legend=paste("df =", 1:4))
})

plotit("lognorm", {
  curve(dlnorm(x), to=3, ylim=c(0, 1.5),
        ylab="probability density", main="Lognormal")
  curve(dlnorm(x, meanlog=1), add=TRUE, col=2)
  curve(dlnorm(x, sdlog=2), add=TRUE, col=3)
  legend("topright", lty=1, col=1:4,
         legend=paste("mu =", c(0, 1, 0), " sd = ", c(1, 1, 2)))
})


