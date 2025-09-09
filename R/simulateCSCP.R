#'
#'  simulateCSCP.R
#'

require(spatstat)
require(spatstat.utils)

simulateCSCP <- function(p, W = owin(), nsim=1, drop=TRUE) {
  if(!is.list(p) || !all(sapply(p, is.list)))
    stop("Argument p should be a list of lists of parameter values",
         call.=FALSE)
  W <- as.owin(W)
  check.1.integer(nsim)
  drop <- isTRUE(drop)
  argh <- list(W=W, nsim=nsim, drop=FALSE) # sic
  ## 
  k <- length(p)
  Zlist <- vector(mode="list", length=k)
  for(i in seq_len(k))
    Zlist[[i]] <- do.call(simulateGRF, append(p[[i]], argh))
  Lamlist <- vector(mode="list", length=nsim)
  Ylist <- vector(mode="list", length=nsim)
  for(j in seq_len(nsim)) {
    Lamj <- 0
    for(i in seq_len(k)) {
      Zij <- Zlist[[i]][[j]]
      Lamj <- Lamj + Zij^2
    }
    Ylist[[j]] <- rpoispp(Lamj)
  }
  result <- simulationresult(Ylist, nsim, drop)
  attr(result, "Lambda") <- simulationresult(Lamlist, nsim, drop)
  return(result)
}

simulateGRF <- function(model=c("exponential", "gauss", "stable",
                                "gencauchy", "matern"),
                        ..., 
                        W=owin(), nsim=1, drop=FALSE) {
  model <- match.arg(model)
  switch(model,
         exponential = {
           rGRFexpo(W=W, ..., nsim=nsim, drop=drop)
         }, gauss = {
           rGRFgauss(W=W, ..., nsim=nsim, drop=drop)
         }, stable = {
           rGRFstable(W=W, ..., nsim=nsim, drop=drop)
         }, gencauchy = {
           rGRFgencauchy(W=W, ..., nsim=nsim, drop=drop)
         }, matern = {
           rGRFmatern(W=W, ..., nsim=nsim, drop=drop)
         }, stop(paste("Model", sQuote(model), "not recognised"),
                 call.=FALSE))
}
  
