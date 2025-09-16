# X: object of class ppp
# k: number of CSCP components (presently only k=1 implemented)
# ...: additional arguments to be passed to pcf.ppp

fitCSCP <- function(X,k=1,plot=FALSE,fitalpha=TRUE,...){
  ghat <- pcf(X,...)
  rr <- ghat$r
  suppressWarnings(gg <- log(ghat$iso-1))
  if(k==1){
    if(fitalpha){
      fit <- lm(gg~rr)
      pars <- coef(fit)
    } else {
      fit <- lm(gg~0+rr+offset(rep(log(2),length(rr))))
      pars <- c(log(2),coef(fit))
    }
  } else {
    stop("k>1 not yet implemented")
  }
  
  if(plot){
    plot(rr,gg,type="l",xlab="r",ylab="log(ghat-1)",main="semilog pcf")
    abline(fit,lty=2,col=2)
  }
  
  return(list(par=c(phi=-2/as.numeric(pars[2]),alpha=sqrt(exp(as.numeric(pars[1]))/2)),lm=fit,ghat=ghat))
}

effrange <- function(phi) -log(0.05)*phi
